#!/usr/bin/env Rscript

# classify_adjusters.R - Single job adjuster comparison script
# Executes single adjuster × classifier × dataset × seed combination

# Suppress warnings and messages for cleaner output
options(warn = -1)
options(repos = c(CRAN = "https://cloud.r-project.org"))
# Do NOT clear workspace if testing_mode exists
if (!exists("testing_mode")) {
  rm(list=ls())
}

# Load required libraries
suppressMessages(suppressWarnings({
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", 
                        "ROCR", "ggplot2", "gridExtra", "reshape2", 
                        "dplyr", "purrr", "nnls", "batchelor",
                        "argparse", "class", "xgboost", "sda", "klaR")
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

parser <- ArgumentParser(description = "Execute single adjuster comparison job for batch correction analysis")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, naive, rank_samples, rank_twice, npn, combat, combat_mean, combat_sup, mnn, fast_mnn, ruvg, yugene, cublock, angel, tdm, rnabc, shambhala2, coconut, rankin, recombat, or recombat_sup")
parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type: rda, elnet, elasticnet, svm, rf, nnet, knn, xgboost, or shrinkageLDA")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets to include: 3, 4, 5, or 6")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Test study name (e.g., GSE37250_SA, USA, India, etc.)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output CSV file path")

# Parse arguments if not sourced
if (!exists("testing_mode")) {
  args <- parser$parse_args()
} else {
  # Default empty args for testing mode to avoid errors later
  args <- list(adjuster="naive", classifier="svm", num_datasets=3, test_study="test", output="test.csv")
}

# Arguments are automatically validated as required by argparse

# Parameter validation
valid_adjusters <- c("unadjusted", "naive", "rank_samples", "rank_twice", "npn", "combat", "combat_mean", "combat_sup", "mnn", "fast_mnn", "ruvg", "ruvr", "yugene", "cublock", "angel", "tdm", "rnabc", "shambhala2", "coconut", "rankin", "recombat", "recombat_sup")
valid_classifiers <- c("rda", "elnet", "elasticnet", "svm", "rf", "nnet", "knn", "xgboost", "shrinkageLDA")
valid_num_datasets <- c(2, 3, 4, 5)

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

if (!args$num_datasets %in% valid_num_datasets) {
  cat(sprintf("Error: Invalid num-datasets '%d'. Must be one of: %s\n", 
              args$num_datasets, paste(valid_num_datasets, collapse=", ")))
  quit(status=1)
}

# Extract validated parameters
adjuster <- args$adjuster
classifier <- args$classifier
num_datasets <- args$num_datasets
test_study <- args$test_study
output_file <- args$output

# Validate output directory exists
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  cat(sprintf("Error: Output directory does not exist: %s\n", output_dir))
  quit(status=1)
}

# ====================================================================
# ERROR HANDLING AND LOGGING WRAPPER
# ====================================================================

# Create job ID for logging
job_id <- sprintf("adjuster_%s_%s_%d_%s", adjuster, classifier, num_datasets, test_study)

# Main job wrapper with comprehensive error handling
main_job_wrapper <- function() {
  tryCatch({
    # Print job parameters for logging
    cat("=== ADJUSTER COMPARISON JOB ===\n")
    cat(sprintf("Job ID: %s\n", job_id))
    cat(sprintf("Adjuster: %s\n", adjuster))
    cat(sprintf("Classifier: %s\n", classifier))
    cat(sprintf("Num datasets: %d\n", num_datasets))
    cat(sprintf("Test study: %s\n", test_study))
    cat(sprintf("Output: %s\n", output_file))
    cat(sprintf("Start time: %s\n", Sys.time()))
    cat("===============================\n\n")
    
    # Execute main analysis
    result <- main_analysis_function()
    
    # Success logging
    cat(sprintf("\n[SUCCESS] Job %s completed at %s\n", job_id, Sys.time()))
    
    return(result)
    
  }, error = function(e) {
    # Detailed error logging
    cat(sprintf("[ERROR] Job %s failed at %s\n", job_id, Sys.time()), file = stderr())
    cat(sprintf("[ERROR] Error: %s\n", e$message), file = stderr())
    cat(sprintf("[ERROR] Parameters: adjuster=%s, classifier=%s, num_datasets=%d, test_study=%s\n", 
                adjuster, classifier, num_datasets, test_study), file = stderr())
    
    # Additional debugging information
    cat(sprintf("[ERROR] Working directory: %s\n", getwd()), file = stderr())
    cat(sprintf("[ERROR] Output file: %s\n", output_file), file = stderr())
    
    # Check input files
    data_path <- "data/TB_real_data.RData"
    helper_path <- "scripts/helper.R"
    
    cat(sprintf("[ERROR] Data file exists: %s\n", file.exists(data_path)), file = stderr())
    cat(sprintf("[ERROR] Helper file exists: %s\n", file.exists(helper_path)), file = stderr())
    
    # Memory usage
    gc_info <- capture.output(gc())
    cat(sprintf("[ERROR] Memory usage: %s\n", paste(gc_info, collapse="; ")), file = stderr())
    
    # Session info
    session_info <- capture.output(sessionInfo())
    cat("[ERROR] Session info:\n", file = stderr())
    cat(paste(session_info[1:5], collapse="\n"), file = stderr())
    cat("\n", file = stderr())
    
    # Exit with error code for Snakemake
    quit(status = 1)
  })
}


# ====================================================================
# MAIN ANALYSIS FUNCTION
# ====================================================================

main_analysis_function <- function() {
  # Load data and dependencies
  data_path <- "data/TB_real_data.RData"
  if (!file.exists(data_path)) {
    stop(sprintf("Data file not found: %s", data_path))
  }
  
  load(data_path)
  source("scripts/helper.R")
  source("scripts/adjusters.R")

  # ====================================================================
  # REAL DATA PREPARATION LOGIC
  # ====================================================================
  
  filter_studies <- function(dat_lst, label_lst, n_studies, test_study) {
    all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
    train_studies <- all_studies[all_studies != test_study][seq_len(n_studies)]
    selected_studies <- c(train_studies, test_study)

    dat_lst <- dat_lst[selected_studies]
    label_lst <- label_lst[selected_studies]
    study_names <- names(dat_lst)
    cat(sprintf("Running %d-study analysis with train: %s, test: %s\n",
                n_studies, paste(train_studies, collapse=", "), test_study))

    list(dat_lst = dat_lst, label_lst = label_lst, study_names = study_names)
  }
  
  #' Prepare training and test data
  prepare_datasets <- function(dat_lst, label_lst, test_name, study_names) {
    train_name <- setdiff(study_names, test_name)
    
    # Ensure all datasets have the same genes by taking intersection
    all_datasets <- c(train_name, test_name)
    common_genes <- Reduce(intersect, lapply(dat_lst[all_datasets], rownames))
    
    cat(sprintf("  Gene intersection: %d common genes across all datasets\n", length(common_genes)))
    
    # Subset all datasets to common genes
    dat_lst_subset <- lapply(dat_lst[all_datasets], function(x) x[common_genes, , drop=FALSE])
    
    # Combine training datasets
    dat <- do.call(cbind, dat_lst_subset[train_name])
    batch <- rep(1:length(train_name), times=sapply(dat_lst_subset[train_name], ncol))
    batches_ind <- lapply(1:length(train_name), function(x) which(batch == x))
    batch_names <- levels(factor(batch))
    group <- unlist(lapply(label_lst[train_name], as.character))
    # Convert back to numeric/factor if needed for classifiers, but COCONUT needs consistent labels
    group <- ifelse(group == "Control" | group == "0", 0, 1)
    y_sgbatch_train <- lapply(batch_names, function(x) group[batch == x])
    
    dat_test <- dat_lst_subset[[test_name]]
    group_test <- as.character(label_lst[[test_name]])
    group_test <- ifelse(group_test == "Control" | group_test == "0", 0, 1)
    
    list(dat=dat, batch=batch, batches_ind=batches_ind, batch_names=batch_names, 
         group=group, y_sgbatch_train=y_sgbatch_train, 
         dat_test=dat_test, group_test=group_test)
  }
  
  #' Select highly variable genes and reduce feature space
  reduce_features <- function(dat, dat_test, n_genes=1000) {
    genes_sel_names <- order(rowVars(dat), decreasing=TRUE)[1:n_genes]
    list(dat=dat[genes_sel_names, ], 
         dat_test=dat_test[genes_sel_names, ])
  }
  
  # ====================================================================
  # EXECUTE DATA PREPARATION
  # ====================================================================
  # Execute main logic
  cat("Starting data preparation...\n")
  
  # Filter studies based on num_datasets parameter
  filtered_data <- filter_studies(dat_lst, label_lst, num_datasets, test_study)
  dat_lst_filtered <- filtered_data$dat_lst
  label_lst_filtered <- filtered_data$label_lst
  study_names <- filtered_data$study_names
  
  # Validate test study is in the filtered list
  if (!test_study %in% study_names) {
    stop(sprintf("Test study '%s' not found in selected studies: %s", 
                 test_study, paste(study_names, collapse=", ")))
  }
  
  test_name <- test_study
  cat(sprintf("Using test study: %s\n", test_name))
  
  # Prepare datasets
  datasets <- prepare_datasets(dat_lst_filtered, label_lst_filtered, test_name, study_names)

  if(is.null(datasets$dat_test)) {
    stop(sprintf("Test dataset '%s' is NULL or missing from data", test_name))
  }
  if(ncol(datasets$dat_test) == 0) {
    stop(sprintf("Test dataset '%s' has no samples", test_name))
  }

  dat      <- datasets$dat
  dat_test <- datasets$dat_test

  dat      <- maybe_log2_transform(dat)
  dat_test <- maybe_log2_transform(dat_test)

  reduced  <- reduce_features(dat, dat_test, n_genes = 1000)
  dat      <- reduced$dat
  dat_test <- reduced$dat_test

  cat(sprintf("Data preparation completed: %d training / %d test samples, %d genes, %d batches\n",
              ncol(dat), ncol(dat_test), nrow(dat), length(unique(datasets$batch))))
  
  # ====================================================================
  # BATCH CORRECTION 
  # ====================================================================
  
  apply_batch_corrections <- function(dat, batch, group, dat_test, method, group_test) {
    handler <- BATCH_CORRECTION_METHODS[[method]]
    if (is.null(handler)) stop(sprintf("Unknown batch correction method: '%s'", method))
    handler(dat = dat, dat_test = dat_test, batch = batch, group = group)
  }
  
  #' Global scaling: scale entire dataset to have overall variance = 1
  #' Preserves relative gene importance while putting data on consistent scale
  #' @param dat_train Training data matrix
  #' @param dat_test Test data matrix
  #' @return List with scaled training and test data
  global_scale <- function(dat_train, dat_test) {
    # Compute global mean and SD from training data
    train_mean <- mean(dat_train, na.rm = TRUE)
    train_sd <- sd(as.vector(dat_train), na.rm = TRUE)
    
    # Check for invalid scaling parameters
    if (is.na(train_mean) || is.na(train_sd) || train_sd == 0 || !is.finite(train_sd) || train_sd < 1e-10) {
      cat(sprintf("Global scaling: mean=%s, sd=%s\n", 
                  ifelse(is.na(train_mean), "NA", sprintf("%.4f", train_mean)),
                  ifelse(is.na(train_sd), "NA", sprintf("%.4f", train_sd))))
      
      # If standard deviation is too small, use unit scaling (subtract mean only)
      if (!is.na(train_mean) && is.finite(train_mean)) {
        cat("[WARNING] Using unit scaling (mean centering only) due to low variance\n")
        dat_train_scaled <- dat_train - train_mean
        dat_test_scaled <- dat_test - train_mean
        
        return(list(
          dat_train = dat_train_scaled,
          dat_test = dat_test_scaled
        ))
      } else {
        stop("Scaling produced invalid values (NA, 0, or infinite standard deviation)")
      }
    }
    
    # Apply same transformation to both train and test
    dat_train_scaled <- (dat_train - train_mean) / train_sd
    dat_test_scaled <- (dat_test - train_mean) / train_sd
    
    # Check for NaN/Inf in results
    if (any(!is.finite(dat_train_scaled)) || any(!is.finite(dat_test_scaled))) {
      stop("Scaling produced non-finite values in output data")
    }
    
    cat(sprintf("Global scaling: mean=%.4f, sd=%.4f\n", train_mean, train_sd))
    
    list(
      dat_train = dat_train_scaled,
      dat_test = dat_test_scaled
    )
  }
  
  # ====================================================================
  # EXECUTE BATCH CORRECTION
  # ====================================================================
  
  cat(sprintf("Applying batch correction method: %s\n", adjuster))
  
  # Apply batch correction
  batch_corr_result <- apply_batch_corrections(dat, datasets$batch, datasets$group, dat_test, adjuster, datasets$group_test)
  dat_corrected <- batch_corr_result$dat_corrected
  dat_test_corrected <- batch_corr_result$dat_test_corrected
  
  # Global scaling (not per-gene normalization)
  cat("Applying global scaling to training and test data...\n")
  
  if(is.null(dat_test_corrected)) {
    stop("Test data is NULL after batch correction")
  }
  
  # Check if batch correction produced valid data
  is_constant_train <- all(dat_corrected == dat_corrected[1,1], na.rm = TRUE)
  is_constant_test <- all(dat_test_corrected == dat_test_corrected[1,1], na.rm = TRUE)
  
  if((!is.na(is_constant_train) && is_constant_train) || (!is.na(is_constant_test) && is_constant_test)) {
    cat("[WARNING] Batch correction produced constant data - using original data instead\n")
    dat_corrected <- dat
    dat_test_corrected <- dat_test
  }
  
  scaled_data <- global_scale(dat_corrected, dat_test_corrected)
  dat_train_norm <- scaled_data$dat_train
  dat_test_norm <- scaled_data$dat_test
  
  # Validation
  if (any(is.na(dat_train_norm)) || any(is.na(dat_test_norm))) {
    stop("Scaling produced NA values")
  }
  if (is.null(dat_test_norm)) {
    stop("Test data became NULL during scaling")
  }
  
  cat(sprintf("Scaled data dimensions - Train: %d x %d, Test: %d x %d\n", 
              nrow(dat_train_norm), ncol(dat_train_norm),
              nrow(dat_test_norm), ncol(dat_test_norm)))
  
  cat(sprintf("Batch correction completed successfully\n"))
  cat(sprintf("  Method: %s\n", adjuster))
  cat(sprintf("  Training data shape: %d x %d\n", nrow(dat_train_norm), ncol(dat_train_norm)))
  cat(sprintf("  Test data shape: %d x %d\n", nrow(dat_test_norm), ncol(dat_test_norm)))
  
  # ====================================================================
  # SINGLE CLASSIFIER TRAINING AND EVALUATION
  # ====================================================================
  
  #' Train single classifier and evaluate performance
  train_and_evaluate_classifier <- function(classifier_type, train_data, train_labels, test_data, test_labels) {
    
    # Initialize variables
    trained_model <- NULL
    test_predictions <- NULL
    
    if (classifier_type == "shrinkageLDA") {
      cat("Training RDA (sda)...\n")
      
      # Transpose data: R (features x samples) -> sda expects (samples x features)
      X_train <- t(train_data)
      X_test <- t(test_data)
      y_train <- as.factor(train_labels)
      
      # Ensure data is in matrix format (sda requires matrix, not data.frame)
      if (!is.matrix(X_train)) {
        X_train <- as.matrix(X_train)
      }
      if (!is.matrix(X_test)) {
        X_test <- as.matrix(X_test)
      }
      
      cat(sprintf("  X_train dimensions: %d x %d (class: %s)\n", nrow(X_train), ncol(X_train), class(X_train)[1]))
      cat(sprintf("  X_test dimensions: %d x %d (class: %s)\n", nrow(X_test), ncol(X_test), class(X_test)[1]))
      
      # Train RDA model
      lda_fit <- sda(
        Xtrain = X_train,
        L = y_train,
        diagonal = FALSE      # Allow correlations between features
      )
      
      # Generate predictions on test set
      pred <- predict(lda_fit, X_test)
      
      # For binary classification, extract probability of positive class
      if (nlevels(y_train) == 2) {
        # Probability of positive class (second column)
        test_predictions <- pred$posterior[, 2]
      } else {
        # Multiclass: pick max probability
        test_predictions <- apply(pred$posterior, 1, max)
      }
      
      trained_model <- list(mod = lda_fit)
      
    } else {
      # Get prediction function for classifier type
      learner_fit <- getPredFunctions(classifier_type)
      
      # Validate data before training
      if (!is.matrix(train_data)) {
        cat(sprintf("Converting train_data to matrix (was %s)\n", class(train_data)))
        train_data <- as.matrix(train_data)
      }
      if (!is.matrix(test_data)) {
        cat(sprintf("Converting test_data to matrix (was %s)\n", class(test_data)))
        test_data <- as.matrix(test_data)
      }
      
      # Train model
      cat(sprintf("Training %s classifier...\n", classifier_type))
      trained_model <- trainPipe(train_set = train_data, train_label = train_labels, 
                                lfit = learner_fit)
      
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
  
  # ====================================================================
  # EXECUTE CLASSIFIER TRAINING AND EVALUATION
  # ====================================================================
  
  cat(sprintf("Training and evaluating classifier: %s\n", classifier))

  if(is.null(dat_test_norm)) {
    stop("Test data is NULL after normalization - this should not happen")
  }
  if(ncol(dat_test_norm) == 0) {
    stop("Test data has zero columns after normalization")
  }
  if(nrow(dat_test_norm) == 0) {
    stop("Test data has zero rows after normalization")
  }
  
  # Classifier-specific early validation
  n_train_samples <- ncol(dat_train_norm)
  n_test_samples <- ncol(dat_test_norm)
  n_features <- nrow(dat_train_norm)
  
  cat(sprintf("Dataset summary before classifier training:\n"))
  cat(sprintf("  Training samples: %d\n", n_train_samples))
  cat(sprintf("  Test samples: %d\n", n_test_samples))
  cat(sprintf("  Features: %d\n", n_features))
  cat(sprintf("  Training labels: %d unique values\n", length(unique(datasets$group))))
  
  # ====================================================================
  # FEATURE REDUCTION FOR HIGH-DIMENSIONAL DATA
  # ====================================================================
  # Classifiers like nnet can't handle 10k+ features. 
  # Reduce to top 1000 most variable genes.
  
  max_features_for_classifier <- 1000
  if (n_features > max_features_for_classifier) {
    cat(sprintf("[FEATURE REDUCTION] Reducing from %d to %d features for classifier\n", 
                n_features, max_features_for_classifier))
    
    # Select top variable genes from training data
    gene_vars <- rowVars(dat_train_norm)
    top_gene_indices <- order(gene_vars, decreasing = TRUE)[1:max_features_for_classifier]
    
    # Apply to both train and test
    dat_train_norm <- dat_train_norm[top_gene_indices, ]
    dat_test_norm <- dat_test_norm[top_gene_indices, ]
    
    n_features <- nrow(dat_train_norm)
    cat(sprintf("[FEATURE REDUCTION] New feature count: %d\n", n_features))
  }
  
  # Early validation for problematic cases
  if(n_train_samples < 10) {
    warning(sprintf("Very small training set (%d samples) - results may be unreliable", n_train_samples))
  }
  

  
  # Check for class imbalance
  class_counts <- table(datasets$group)
  if(min(class_counts) < 3) {
    warning(sprintf("Severe class imbalance detected: %s", paste(names(class_counts), class_counts, sep="=", collapse=", ")))
  }
  
  # Train classifier and evaluate performance
  result <- train_and_evaluate_classifier(
    classifier_type = classifier,
    train_data = dat_train_norm,
    train_labels = datasets$group,
    test_data = dat_test_norm,
    test_labels = datasets$group_test
  )
  
  cat("Classification completed successfully\n")
  cat(sprintf("Performance metrics for %s + %s:\n", adjuster, classifier))
  for (metric in names(result$performance)) {
    cat(sprintf("  %s: %.6f\n", metric, result$performance[metric]))
  }
  
  # ====================================================================
  # CSV OUTPUT FORMAT IMPLEMENTATION
  # ====================================================================
  
  #' Create output data frame with required columns
  create_output_dataframe <- function(adjuster, classifier, n_datasets, test_study, performance_metrics) {
    # Create one row per metric
    output_rows <- lapply(names(performance_metrics), function(metric) {
      data.frame(
        adjuster = adjuster,
        classifier = classifier,
        n_datasets = n_datasets,
        test_study = test_study,
        metric = metric,
        value = performance_metrics[metric],
        stringsAsFactors = FALSE
      )
    })
    
    # Combine all rows
    output_df <- do.call(rbind, output_rows)
    return(output_df)
  }
  
  # ====================================================================
  # GENERATE AND WRITE OUTPUT
  # ====================================================================
  
  cat("Generating output CSV...\n")
  
  # Create output data frame
  output_df <- create_output_dataframe(
    adjuster = adjuster,
    classifier = classifier,
    n_datasets = num_datasets,
    test_study = test_study,
    performance_metrics = result$performance
  )
  
  # Write to CSV file with error handling
  tryCatch({
    write.csv(output_df, file = output_file, row.names = FALSE)
    
    # Verify file was created
    if (!file.exists(output_file)) {
      stop(sprintf("File was not created: %s", output_file))
    }
    
    # Verify file has content
    file_size <- file.info(output_file)$size
    if (is.na(file_size) || file_size == 0) {
      stop(sprintf("File was created but is empty: %s", output_file))
    }
    
    cat(sprintf("Results written to: %s\n", output_file))
    cat(sprintf("Output contains %d rows (one per metric)\n", nrow(output_df)))
    cat(sprintf("File size: %d bytes\n", file_size))
    
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to write output file: %s\n", e$message), file = stderr())
    cat(sprintf("[ERROR] Output file path: %s\n", output_file), file = stderr())
    cat(sprintf("[ERROR] Output directory exists: %s\n", dir.exists(dirname(output_file))), file = stderr())
    cat(sprintf("[ERROR] Output directory writable: %s\n", file.access(dirname(output_file), 2) == 0), file = stderr())
    stop(sprintf("Failed to write output file: %s", e$message))
  })
  
  # Display output for verification
  cat("\nOutput preview:\n")
  print(output_df)
  
  cat(sprintf("\n=== JOB COMPLETED SUCCESSFULLY ===\n"))
  cat(sprintf("Adjuster: %s\n", adjuster))
  cat(sprintf("Classifier: %s\n", classifier))
  cat(sprintf("Datasets: %d\n", num_datasets))
  cat(sprintf("Test study: %s\n", test_study))
  cat(sprintf("Output: %s\n", output_file))
  cat("==================================\n")
  
  return(output_df)
}

# ====================================================================
# EXECUTE MAIN JOB
# ====================================================================

# Execute main logic if not in testing mode
if (!exists("testing_mode")) {
  res <- main_job_wrapper()
  quit(save = "no", status = 0)
}
