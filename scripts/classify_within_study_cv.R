#!/usr/bin/env Rscript

# classify_within_study_cv.R - Within-study cross-validation baseline
# Performs stratified k-fold CV within a single study to establish upper bound performance

options(warn = -1)
suppressMessages(suppressWarnings({
  rm(list=ls())
}))

# Load required libraries
suppressMessages(suppressWarnings({
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", 
                        "ROCR", "ggplot2", "gridExtra", "reshape2", 
                        "dplyr", "purrr", "nnls", "batchelor",
                        "argparse", "class", "xgboost", "caret", "sda")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Within-study cross-validation baseline")

parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets (for consistency with other rules)")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Study name to perform CV on")
parser$add_argument("--n-folds", type = "integer", default = 3,
                   help = "Number of CV folds (default: 3)")
parser$add_argument("--n-features", type = "integer", default = 0,
                   help = "Number of top variable genes to select (0 = use all genes, default: 0)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output CSV file path")

args <- parser$parse_args()

classifier <- args$classifier
num_datasets <- args$num_datasets
test_study <- args$test_study
n_folds <- args$n_folds
n_features <- args$n_features
output_file <- args$output

valid_classifiers <- c("rda", "logistic", "elnet", "elasticnet", "svm", "rf", "nnet", "knn", "xgboost", "shrinkageLDA")

if (!classifier %in% valid_classifiers) {
  cat(sprintf("Error: Invalid classifier '%s'\n", classifier))
  quit(status=1)
}

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  cat(sprintf("Error: Output directory does not exist: %s\n", output_dir))
  quit(status=1)
}

# ====================================================================
# MAIN ANALYSIS
# ====================================================================

main_analysis <- function() {
  cat("=== WITHIN-STUDY CV BASELINE ===\n")
  cat(sprintf("Classifier: %s\n", classifier))
  cat(sprintf("Test study: %s\n", test_study))
  cat(sprintf("N folds: %d\n", n_folds))
  cat(sprintf("N features: %s\n", if(n_features == 0) "all" else n_features))
  cat(sprintf("Start time: %s\n", Sys.time()))
  cat("================================\n\n")
  
  # Load data
  data_path <- "data/TB_real_data.RData"
  if (!file.exists(data_path)) {
    stop(sprintf("Data file not found: %s", data_path))
  }
  
  load(data_path)
  source("scripts/helper.R")
  
  # Extract the test study data
  if (!test_study %in% names(dat_lst)) {
    stop(sprintf("Test study '%s' not found in data", test_study))
  }
  
  dat <- dat_lst[[test_study]]
  labels <- label_lst[[test_study]]
  
  cat(sprintf("Loaded study: %s\n", test_study))
  cat(sprintf("  Samples: %d\n", ncol(dat)))
  cat(sprintf("  Genes: %d\n", nrow(dat)))
  cat(sprintf("  TB cases: %d\n", sum(labels == 1)))
  cat(sprintf("  Controls: %d\n", sum(labels == 0)))
  
  # Optional feature reduction
  if (n_features > 0 && n_features < nrow(dat)) {
    n_highvar_genes <- min(n_features, nrow(dat))
    genes_sel <- order(rowVars(dat), decreasing=TRUE)[1:n_highvar_genes]
    dat <- dat[genes_sel, ]
    cat(sprintf("Reduced to %d most variable genes\n", n_highvar_genes))
  } else {
    cat(sprintf("Using all %d genes (no feature selection)\n", nrow(dat)))
  }
  
  # Global scaling
  dat_mean <- mean(dat)
  dat_sd <- sd(as.vector(dat))
  dat_scaled <- (dat - dat_mean) / dat_sd
  
  cat(sprintf("Global scaling: mean=%.4f, sd=%.4f\n", dat_mean, dat_sd))
  
  # Create stratified folds
  set.seed(42)  # For reproducibility
  folds <- createFolds(labels, k = n_folds, list = TRUE, returnTrain = FALSE)
  
  cat(sprintf("\nCreated %d stratified folds:\n", n_folds))
  for (i in 1:n_folds) {
    fold_labels <- labels[folds[[i]]]
    cat(sprintf("  Fold %d: %d samples (%d TB, %d Control)\n", 
                i, length(folds[[i]]), sum(fold_labels == 1), sum(fold_labels == 0)))
  }
  
  # Perform CV
  all_predictions <- numeric(ncol(dat))
  all_labels <- labels
  
  cat("\nPerforming cross-validation...\n")
  
  for (fold_idx in 1:n_folds) {
    cat(sprintf("\n--- Fold %d/%d ---\n", fold_idx, n_folds))
    
    test_idx <- folds[[fold_idx]]
    train_idx <- setdiff(1:ncol(dat), test_idx)
    
    dat_train <- dat_scaled[, train_idx]
    dat_test <- dat_scaled[, test_idx]
    labels_train <- labels[train_idx]
    labels_test <- labels[test_idx]
    
    cat(sprintf("Train: %d samples, Test: %d samples\n", 
                length(train_idx), length(test_idx)))
    
    # Train classifier
    if (classifier == "shrinkageLDA") {
      # RDA using sda package
      X_train <- t(dat_train)  # Transpose to samples x features
      X_test <- t(dat_test)
      y_train <- as.factor(labels_train)
      
      # Ensure data is in matrix format (sda requires matrix, not data.frame)
      if (!is.matrix(X_train)) {
        X_train <- as.matrix(X_train)
      }
      if (!is.matrix(X_test)) {
        X_test <- as.matrix(X_test)
      }
      
      # Train RDA model
      lda_fit <- sda(
        Xtrain = X_train,
        L = y_train,
        diagonal = FALSE
      )
      
      # Generate predictions
      pred <- predict(lda_fit, X_test)
      
      # Extract probabilities for binary classification
      if (nlevels(y_train) == 2) {
        predictions <- pred$posterior[, 2]  # Probability of positive class
      } else {
        predictions <- apply(pred$posterior, 1, max)
      }
      
    } else if (classifier == "nnet" || classifier == "nn") {
      # Neural network with increased MaxNWts for all genes
      trn_transposed <- t(dat_train)
      rownames(trn_transposed) <- NULL
      data <- data.frame(y=as.factor(labels_train), trn_transposed)
      
      n_samples <- nrow(data)
      network_size <- min(10, max(2, n_samples %/% 3))
      
      # Calculate required weights: (n_features + 1) * hidden_size + (hidden_size + 1) * n_outputs
      n_features <- ncol(dat_train)
      max_weights <- (n_features + 1) * network_size + (network_size + 1) * 2
      max_weights <- max(max_weights + 1000, 150000)  # Add buffer, minimum 150k
      
      mod <- nnet::nnet(y ~ ., data = data, size = network_size, MaxNWts = max_weights, 
                  decay = 0.01, linout = FALSE, trace = FALSE, maxit = 200)
      
      tst_transposed <- t(dat_test)
      rownames(tst_transposed) <- NULL
      newdata <- data.frame(tst_transposed)
      predictions <- as.vector(predict(mod, newdata = newdata))
      
    } else {
      # Standard classifiers
      learner_fit <- getPredFunctions(classifier)
      trained_model <- trainPipe(train_set = dat_train, train_label = labels_train, 
                                 lfit = learner_fit)
      predictions <- predWrapper(trained_model$mod, dat_test, classifier)
    }
    
    # Store predictions
    all_predictions[test_idx] <- predictions
    
    cat(sprintf("Fold %d complete\n", fold_idx))
  }
  
  cat("\n=== Cross-validation complete ===\n")
  
  # Calculate performance metrics
  perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc")
  
  predictions_list <- list()
  predictions_list[["within_study_cv"]] <- all_predictions
  
  perf_results <- perf_wrapper(perf_measures, predictions_list, all_labels)
  confusion_results <- confusion_matrix_wrapper(predictions_list, all_labels)
  
  combined_results <- rbind(perf_results, confusion_results)
  perf_values <- combined_results[, "within_study_cv"]
  names(perf_values) <- rownames(combined_results)
  
  cat("\nPerformance metrics:\n")
  for (metric in names(perf_values)) {
    cat(sprintf("  %s: %.6f\n", metric, perf_values[metric]))
  }
  
  # Create output dataframe
  output_rows <- lapply(names(perf_values), function(metric) {
    data.frame(
      adjuster = "within_study_cv",
      classifier = classifier,
      n_datasets = num_datasets,
      test_study = test_study,
      metric = metric,
      value = perf_values[metric],
      stringsAsFactors = FALSE
    )
  })
  
  output_df <- do.call(rbind, output_rows)
  
  # Write output
  write.csv(output_df, file = output_file, row.names = FALSE)
  
  if (!file.exists(output_file)) {
    stop(sprintf("File was not created: %s", output_file))
  }
  
  cat(sprintf("\nResults written to: %s\n", output_file))
  cat(sprintf("Output contains %d rows\n", nrow(output_df)))
  
  cat("\n=== JOB COMPLETED SUCCESSFULLY ===\n")
  
  return(output_df)
}

# Execute with error handling
tryCatch({
  result <- main_analysis()
  quit(save = "no", status = 0)
}, error = function(e) {
  cat(sprintf("[ERROR] Job failed: %s\n", e$message), file = stderr())
  cat(sprintf("[ERROR] Classifier: %s, Test study: %s\n", classifier, test_study), file = stderr())
  quit(status = 1)
})
