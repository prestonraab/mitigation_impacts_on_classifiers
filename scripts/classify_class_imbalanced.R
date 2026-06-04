#!/usr/bin/env Rscript

# classify_class_imbalanced.R - Classification on systematically class-imbalanced TB data
# Tests Combat supervised vs unsupervised vs unadjusted across multiple imbalance levels

suppressMessages(suppressWarnings({
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", 
                        "ROCR", "ggplot2", "gridExtra", "reshape2", 
                        "dplyr", "purrr", "nnls", "batchelor",
                        "argparse", "class", "xgboost", "sda", "reticulate", "COCONUT")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# Configure reticulate to use the pixi environment Python
if (requireNamespace("reticulate", quietly = TRUE)) {
  # Explicitly disable automatic virtualenv/conda env creation
  Sys.setenv(RETICULATE_AUTOCREATE_VENV = "FALSE")
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
  
  # Find the python in the current pixi environment
  pixi_python <- file.path(getwd(), ".pixi/envs/default/bin/python")
  if (file.exists(pixi_python)) {
    reticulate::use_python(pixi_python, required = TRUE)
  }
}


# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Execute class-imbalanced TB classification")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, combat, or combat_sup")
parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type: logistic, elnet, elasticnet, svm, rf, nnet, knn, xgboost, or shrinkageLDA")
parser$add_argument("--scenario-name", type = "character", required = TRUE,
                   help = "Scenario name (e.g., GSE37250_SA-USA-imbal20-testIndia-rep1)")
parser$add_argument("--data", type = "character", required = TRUE,
                   help = "Path to TB_real_data.RData")
parser$add_argument("--seed", type = "integer", default = 123,
                   help = "Base random seed for reproducible subsetting (default: 123)")
parser$add_argument("--imbalance-levels", type = "character", required = TRUE,
                   help = "Comma-separated imbalance proportions used to build scenarios, e.g. '0.1,0.2,0.3,0.4,0.5'")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output CSV file path")

args <- parser$parse_args()

# Parameter validation
valid_adjusters <- c("unadjusted", "combat", "combat_sup", "rankin", "coconut", "naive", "rank_twice", "recombat", "recombat_sup")
valid_classifiers <- c("rda", "logistic", "elnet", "elasticnet", "svm", "rf", "nnet", "knn", "xgboost", "shrinkageLDA")

if (!args$adjuster %in% valid_adjusters) {
  cat(sprintf("Error: Invalid adjuster '%s'. Must be one of: %s\n",
              args$adjuster, paste(valid_adjusters, collapse=", ")))
  quit(status=1)
}

if (!args$classifier %in% valid_classifiers) {
  cat(sprintf("Error: Invalid classifier '%s'. Must be one of: %s\n", 
              args$classifier, paste(valid_classifiers, collapse=", ")))
  quit(status=1)
}

adjuster <- args$adjuster
classifier <- args$classifier
scenario_name <- args$scenario_name
data_file <- args$data
base_seed <- args$seed
imbalance_levels <- as.numeric(strsplit(args$imbalance_levels, ",")[[1]])
output_file <- args$output

# ====================================================================
# HELPER FUNCTIONS
# ====================================================================

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")

# Parse scenario name — derives training pair by matching against all combn pairs
parse_scenario_name <- function(scenario_name) {
  parts <- strsplit(scenario_name, "-imbal|-test|-rep")[[1]]
  if (length(parts) < 3) {
    stop(sprintf("Cannot parse scenario name: %s", scenario_name))
  }

  training_pair_str <- parts[1]
  imbalance_pct     <- as.numeric(parts[2]) / 100
  test_dataset      <- parts[3]
  replicate         <- if (length(parts) >= 4) as.numeric(parts[4]) else 1

  training_pairs_list <- combn(ALL_STUDIES, 2, simplify = FALSE)
  pair_names <- sapply(training_pairs_list, paste, collapse = "-")
  pair_idx   <- which(pair_names == training_pair_str)
  if (length(pair_idx) == 0) {
    stop(sprintf("Cannot find training pair '%s' among all combinations", training_pair_str))
  }
  train_datasets <- training_pairs_list[[pair_idx]]

  return(list(
    train_datasets = train_datasets,
    test_dataset   = test_dataset,
    imbalance_pct  = imbalance_pct,
    training_pair  = training_pair_str,
    replicate      = replicate,
    pair_idx       = pair_idx
  ))
}

# ── Subsetting helpers (mirror create_class_imbalanced_data.R exactly) ──────

calculate_optimal_samples <- function(n_active, n_latent, target_active_pct, target_total) {
  target_active <- round(target_total * target_active_pct)
  target_latent <- target_total - target_active
  if (target_active > n_active || target_latent > n_latent ||
      target_active < 3 || target_latent < 3) {
    return(NULL)
  }
  list(total = target_total, active = target_active, latent = target_latent,
       actual_active_pct = target_active / target_total)
}

determine_optimal_arrangement_consistent <- function(ds1_stats, ds2_stats, all_imbalance_pcts) {
  max_totals <- function(hi, lo) {
    sapply(all_imbalance_pcts, function(p)
      min(floor(hi$n_active / p), floor(hi$n_latent / (1 - p)),
          floor(lo$n_active / (1 - p)), floor(lo$n_latent / p)))
  }
  t1 <- min(max_totals(ds1_stats, ds2_stats))
  t2 <- min(max_totals(ds2_stats, ds1_stats))
  if (t1 >= t2 && t1 >= 30)
    list(high_active_dataset = ds1_stats$dataset, low_active_dataset = ds2_stats$dataset,
         total_samples = t1, feasible = TRUE)
  else if (t2 >= 30)
    list(high_active_dataset = ds2_stats$dataset, low_active_dataset = ds1_stats$dataset,
         total_samples = t2, feasible = TRUE)
  else
    list(feasible = FALSE)
}

create_imbalanced_subset <- function(data, labels, target_active_pct, target_total,
                                     seed_offset = 0, base_seed = 123) {
  set.seed(base_seed + seed_offset)
  active_idx <- which(labels == 1)
  latent_idx <- which(labels == 0)
  optimal <- calculate_optimal_samples(length(active_idx), length(latent_idx),
                                       target_active_pct, target_total)
  if (is.null(optimal)) {
    stop(sprintf("Cannot create subset: infeasible target %.0f%% with total %.0f",
                 target_active_pct * 100, target_total))
  }
  keep <- sort(c(sample(active_idx, optimal$active), sample(latent_idx, optimal$latent)))
  list(data = data[, keep, drop = FALSE], labels = labels[keep])
}

# Define the train_and_evaluate_classifier function
train_and_evaluate_classifier <- function(classifier_type, train_data, train_labels, test_data, test_labels) {
  
  source("scripts/helper.R")
  source("scripts/adjusters.R")
  
  # Initialize variables
  trained_model <- NULL
  test_predictions <- NULL
  
  if (classifier_type == "shrinkageLDA") {
    cat("Training RDA (sda)...\n")
    
    # Transpose data: R (features x samples) -> sda expects (samples x features)
    X_train <- t(train_data)
    X_test <- t(test_data)
    y_train <- as.factor(train_labels)
    
    # Ensure data is in matrix format
    if (!is.matrix(X_train)) X_train <- as.matrix(X_train)
    if (!is.matrix(X_test)) X_test <- as.matrix(X_test)
    
    # Train RDA model
    lda_fit <- sda(Xtrain = X_train, L = y_train, diagonal = FALSE)
    
    # Generate predictions on test set
    pred <- predict(lda_fit, X_test)
    
    # For binary classification, extract probability of positive class
    if (nlevels(y_train) == 2) {
      test_predictions <- pred$posterior[, 2]
    } else {
      test_predictions <- apply(pred$posterior, 1, max)
    }
    
    trained_model <- list(mod = lda_fit)
    
  } else {
    # Get prediction function for classifier type
    learner_fit <- getPredFunctions(classifier_type)
    
    # Validate data before training
    if (!is.matrix(train_data)) train_data <- as.matrix(train_data)
    if (!is.matrix(test_data)) test_data <- as.matrix(test_data)
    
    # Train model
    cat(sprintf("Training %s classifier...\n", classifier_type))
    trained_model <- trainPipe(train_set = train_data, train_label = train_labels, lfit = learner_fit)
    
    # Generate predictions on test set
    cat(sprintf("Generating predictions...\n"))
    test_predictions <- predWrapper(trained_model$mod, test_data, classifier_type)
  }
  
  # Calculate performance metrics
  perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc")
  
  # Create predictions list in format expected by perf_wrapper
  predictions_list <- list()
  predictions_list[[adjuster]] <- test_predictions
  
  # Calculate performance using original perf_wrapper function
  perf_results <- perf_wrapper(perf_measures, predictions_list, test_labels)
  
  # Calculate confusion matrix elements and derived metrics
  confusion_results <- confusion_matrix_wrapper(predictions_list, test_labels)
  
  # Combine performance metrics and confusion matrix metrics
  combined_results <- rbind(perf_results, confusion_results)
  
  # Extract performance values for this method
  perf_values <- combined_results[, adjuster]
  names(perf_values) <- rownames(combined_results)
  
  return(list(
    model = trained_model,
    predictions = test_predictions,
    performance = perf_values
  ))
}

# Adjuster helper functions
adjust_coconut <- function(matrix_, batch, group, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Starting COCONUT harmonization.\n")
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
    cat("DEBUG: Unique batches: ", length(unique(batch)), "\n")
  }
  
  if (!requireNamespace("COCONUT", quietly = TRUE)) {
    stop("Package 'COCONUT' is required but not installed.")
  }
  
  gse_list <- list()
  for (b in unique(batch)) {
    idx <- which(batch == b)
    disease_vec <- as.numeric(group[idx])
    
    pheno_df <- data.frame(
      disease_state = disease_vec,
      dummy = 1, # Add dummy column to prevent dimension dropping
      row.names = colnames(matrix_[, idx])
    )
    
    gse_list[[as.character(b)]] <- list(
      pheno = pheno_df,
      genes = matrix_[, idx, drop = FALSE]
    )
  }
  
  res <- COCONUT::COCONUT(GSEs = gse_list, control.0.col = "disease_state")
  
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  
  for (b in names(res$COCONUTList)) {
    disease_cols <- colnames(res$COCONUTList[[b]]$genes)
    result_matrix[, disease_cols] <- as.matrix(res$COCONUTList[[b]]$genes)
  }
  for (b in names(res$controlList$GSEs)) {
    control_cols <- colnames(res$controlList$GSEs[[b]]$genes)
    result_matrix[, control_cols] <- as.matrix(res$controlList$GSEs[[b]]$genes)
  }
  
  return(result_matrix)
}

adjust_rankin <- function(matrix_, n_svd = 1L, train_indices = NULL, debug = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("LIBRARY MISSING: 'reticulate' is required for Rank-In.")
  }
  
  # Locate rankin.py
  rankin_script <- "scripts/rankin.py"
  if (!file.exists(rankin_script)) {
    rankin_script <- "rankin.py"
  }
  if (!file.exists(rankin_script)) {
    stop("CODE MISSING: Rank-In script (rankin.py) not found in scripts/ or current directory.")
  }
  
  reticulate::source_python(rankin_script)
  
  # Pass matrix and train_indices to Python (reticulate)
  result_list <- rank_in_from_r(
    unname(as.list(as.data.frame(t(matrix_)))), 
    train_indices = train_indices,
    n_svd = as.integer(n_svd)
  )
  
  # Reconstruct matrix (genes x samples) row-wise as it is returned as a list of gene vectors
  result_matrix <- matrix(unlist(result_list), nrow = nrow(matrix_), ncol = ncol(matrix_), byrow = TRUE)
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  
  return(result_matrix)
}


# ====================================================================
# MAIN ANALYSIS FUNCTION
# ====================================================================

main_class_imbalanced_analysis <- function() {
  
  cat("=== CLASS-IMBALANCED TB CLASSIFICATION ===\n")
  cat(sprintf("Adjuster: %s\n", adjuster))
  cat(sprintf("Classifier: %s\n", classifier))
  cat(sprintf("Scenario: %s\n", scenario_name))
  cat("==========================================\n\n")
  
  # Load full TB data and subset inline (no intermediate files)
  if (!file.exists(data_file)) {
    stop(sprintf("Data file not found: %s", data_file))
  }
  load(data_file)  # loads dat_lst and label_lst (all 6 datasets)
  dat_lst_original   <- dat_lst
  label_lst_original <- label_lst

  scenario_info <- parse_scenario_name(scenario_name)

  # Compute arrangement and seed offset (mirrors create_class_imbalanced_data.R)
  train_datasets <- scenario_info$train_datasets
  test_dataset   <- scenario_info$test_dataset
  imbalance_pct  <- scenario_info$imbalance_pct
  replicate_idx  <- scenario_info$replicate

  original_stats <- lapply(ALL_STUDIES, function(ds) {
    lbls <- label_lst_original[[ds]]
    list(dataset = ds, n_active = sum(lbls == 1), n_latent = sum(lbls == 0))
  })
  names(original_stats) <- ALL_STUDIES

  arrangement <- determine_optimal_arrangement_consistent(
    original_stats[[train_datasets[1]]], original_stats[[train_datasets[2]]], imbalance_levels
  )
  if (!arrangement$feasible) {
    stop(sprintf("Infeasible arrangement for pair %s", scenario_info$training_pair))
  }

  high_ds <- arrangement$high_active_dataset
  low_ds  <- arrangement$low_active_dataset
  total_n <- arrangement$total_samples

  test_datasets_for_pair <- setdiff(ALL_STUDIES, train_datasets)
  seed_offset <- scenario_info$pair_idx * 100000 +
                 which(imbalance_levels == imbalance_pct) * 10000 +
                 which(test_datasets_for_pair == test_dataset) * 100 +
                 replicate_idx

  high_sub <- create_imbalanced_subset(dat_lst_original[[high_ds]],
                                       label_lst_original[[high_ds]],
                                       imbalance_pct, total_n,
                                       seed_offset = seed_offset, base_seed = base_seed)
  low_sub  <- create_imbalanced_subset(dat_lst_original[[low_ds]],
                                       label_lst_original[[low_ds]],
                                       1 - imbalance_pct, total_n,
                                       seed_offset = seed_offset + 50000, base_seed = base_seed)

  dat_lst   <- list()
  label_lst <- list()
  dat_lst[[high_ds]]   <- high_sub$data;   label_lst[[high_ds]]   <- high_sub$labels
  dat_lst[[low_ds]]    <- low_sub$data;    label_lst[[low_ds]]    <- low_sub$labels
  dat_lst[[test_dataset]]   <- dat_lst_original[[test_dataset]]
  label_lst[[test_dataset]] <- label_lst_original[[test_dataset]]
  
  cat(sprintf("Parsed scenario:\n"))
  cat(sprintf("  Training pair: %s\n", scenario_info$training_pair))
  cat(sprintf("  Train datasets: %s\n", paste(scenario_info$train_datasets, collapse = ", ")))
  cat(sprintf("  Test dataset: %s\n", scenario_info$test_dataset))
  cat(sprintf("  Imbalance level: %.0f%%\n", scenario_info$imbalance_pct * 100))
  cat(sprintf("  Replicate: %d\n", scenario_info$replicate))
  
  # Validate that all required datasets are present
  required_datasets <- c(scenario_info$train_datasets, scenario_info$test_dataset)
  missing_datasets <- setdiff(required_datasets, names(dat_lst))
  if (length(missing_datasets) > 0) {
    stop(sprintf("Missing datasets in loaded data: %s", paste(missing_datasets, collapse = ", ")))
  }
  
  # Prepare training and test data
  train_datasets <- scenario_info$train_datasets
  test_dataset <- scenario_info$test_dataset
  
  # Ensure all datasets have the same genes by taking intersection
  all_datasets <- c(train_datasets, test_dataset)
  common_genes <- Reduce(intersect, lapply(dat_lst[all_datasets], rownames))
  
  cat(sprintf("Gene intersection: %d common genes across all datasets\n", length(common_genes)))
  
  # Prepare training data (combine training datasets)
  dat_lst_train_subset <- lapply(dat_lst[train_datasets], function(x) x[common_genes, , drop=FALSE])
  dat <- do.call(cbind, dat_lst_train_subset)
  batch <- rep(1:length(train_datasets), times=sapply(dat_lst_train_subset, ncol))
  group <- do.call(c, label_lst[train_datasets])
  
  # Prepare test data
  dat_test <- dat_lst[[test_dataset]][common_genes, , drop=FALSE]
  group_test <- label_lst[[test_dataset]]
  
  # Feature reduction (top 1000 most variable genes)
  n_highvar_genes <- min(1000, nrow(dat))
  genes_sel_names <- order(rowVars(dat), decreasing=TRUE)[1:n_highvar_genes]
  dat <- dat[genes_sel_names, ]
  dat_test <- dat_test[genes_sel_names, ]
  
  # Log transformation
  needs_log_transform_train <- max(dat) > 100 || mean(dat) > 20 || (max(dat) / median(dat)) > 50
  needs_log_transform_test <- max(dat_test) > 100 || mean(dat_test) > 20 || (max(dat_test) / median(dat_test)) > 50
  
  if (needs_log_transform_train) {
    train_min <- min(dat)
    if (train_min < 0) dat <- dat - train_min
    dat <- log2(dat + 1)
  }
  
  if (needs_log_transform_test) {
    test_min <- min(dat_test)
    if (test_min < 0) dat_test <- dat_test - test_min
    dat_test <- log2(dat_test + 1)
  }
  
  cat(sprintf("Data preparation completed:\n"))
  cat(sprintf("  Training samples: %d\n", ncol(dat)))
  cat(sprintf("  Test samples: %d\n", ncol(dat_test)))
  cat(sprintf("  Features (genes): %d\n", nrow(dat)))
  cat(sprintf("  Training batches: %d\n", length(unique(batch))))
  
  # Print class distribution to show the imbalance effect
  train_active <- sum(group == 1)
  train_latent <- sum(group == 0)
  test_active <- sum(group_test == 1)
  test_latent <- sum(group_test == 0)
  
  cat(sprintf("TRAINING class distribution: LTBI=%d (%.1f%%), Active TB=%d (%.1f%%)\n", 
              train_latent, 100 * train_latent / (train_active + train_latent),
              train_active, 100 * train_active / (train_active + train_latent)))
  cat(sprintf("TEST class distribution: LTBI=%d (%.1f%%), Active TB=%d (%.1f%%)\n", 
              test_latent, 100 * test_latent / (test_active + test_latent),
              test_active, 100 * test_active / (test_active + test_latent)))
  
  # Print batch-specific class distributions
  cat("Training batch-specific class distributions:\n")
  for (i in 1:length(train_datasets)) {
    batch_mask <- batch == i
    batch_group <- group[batch_mask]
    batch_active <- sum(batch_group == 1)
    batch_latent <- sum(batch_group == 0)
    batch_total <- batch_active + batch_latent

    cat(sprintf("  %s (batch %d): LTBI=%d (%.1f%%), Active TB=%d (%.1f%%), Total=%d\n",
                train_datasets[i], i, batch_latent, 100 * batch_latent / batch_total,
                batch_active, 100 * batch_active / batch_total, batch_total))
  }

  # ====================================================================
  # BATCH CORRECTION
  # ====================================================================

  cat(sprintf("Applying batch correction method: %s\n", adjuster))
  handler <- BATCH_CORRECTION_METHODS[[adjuster]]
  if (is.null(handler)) stop(sprintf("Unknown batch correction method: '%s'", adjuster))
  batch_result       <- handler(dat = dat, dat_test = dat_test, batch = batch, group = group)
  dat_corrected      <- batch_result$dat_corrected
  dat_test_corrected <- batch_result$dat_test_corrected

  # ====================================================================
  # CLASSIFICATION
  # ====================================================================
  
  cat(sprintf("Training %s classifier...\n", classifier))
  
  # Train and evaluate classifier
  result <- train_and_evaluate_classifier(
    classifier_type = classifier,
    train_data = dat_corrected,
    train_labels = group,
    test_data = dat_test_corrected,
    test_labels = group_test
  )
  
  # ====================================================================
  # SAVE RESULTS
  # ====================================================================
  
  # Extract performance metrics from the result
  performance <- result$performance
  
  # Calculate batch-class confounding metrics
  batch_tb_ratios <- sapply(1:length(train_datasets), function(i) {
    batch_mask <- batch == i
    batch_group <- group[batch_mask]
    sum(batch_group == 1) / length(batch_group)
  })
  confounding_strength <- sd(batch_tb_ratios)
  
  # Create results dataframe
  results_df <- data.frame(
    scenario_name = scenario_name,
    training_pair = scenario_info$training_pair,
    test_dataset = scenario_info$test_dataset,
    train_dataset_1 = train_datasets[1],
    train_dataset_2 = train_datasets[2],
    imbalance_pct = scenario_info$imbalance_pct,
    replicate = scenario_info$replicate,
    adjuster = adjuster,
    classifier = classifier,
    mcc = performance["mcc"],
    accuracy = performance["acc"],
    sensitivity = performance["sensitivity"],
    specificity = performance["specificity"],
    auc = performance["auc"],
    f1_score = performance["f"],
    train_samples = ncol(dat),
    test_samples = ncol(dat_test),
    train_active = train_active,
    train_latent = train_latent,
    test_active = test_active,
    test_latent = test_latent,
    train_active_pct = train_active / (train_active + train_latent),
    test_active_pct = test_active / (test_active + test_latent),
    batch_confounding_strength = confounding_strength,
    analysis_type = "class_imbalanced",
    stringsAsFactors = FALSE
  )
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(results_df, output_file, row.names = FALSE)
  
  cat(sprintf("Results saved to: %s\n", output_file))
  cat(sprintf("MCC: %.4f, Accuracy: %.4f, AUC: %.4f\n", 
              performance["mcc"], performance["acc"], performance["auc"]))
  
  return(results_df)
}

# ====================================================================
# EXECUTE WITH ERROR HANDLING
# ====================================================================

tryCatch({
  result <- main_class_imbalanced_analysis()
}, error = function(e) {
  cat(sprintf("[ERROR] %s\n", e$message), file = stderr())
  quit(status = 1)
})