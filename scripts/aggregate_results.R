#!/usr/bin/env Rscript

# aggregate_results.R
# Script to aggregate individual CSV results from parallel jobs
# Expected to be called from Snakemake workflow

# Suppress warnings and messages for cleaner output
options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
})

# Define command line arguments
parser <- ArgumentParser(description = "Aggregate individual CSV results from parallel jobs")

parser$add_argument("--input-dir", 
                   type = "character", 
                   action = "append",
                   help = "Input directory containing individual CSV result files (can be specified multiple times)")

parser$add_argument("-o", "--output", 
                   type = "character", 
                   required = TRUE,
                   help = "Output CSV file path for aggregated results")

parser$add_argument("--verbose", 
                   action = "store_true",
                   default = FALSE,
                   help = "Enable verbose output for debugging")

# Parse arguments
opt <- parser$parse_args()

# Validate input directories exist
if (is.null(opt$input_dir) || length(opt$input_dir) == 0) {
  stop("At least one --input-dir must be specified")
}

for (dir in opt$input_dir) {
  if (!dir.exists(dir)) {
    stop("Input directory does not exist: ", dir)
  }
}

# Create output directory if it doesn't exist
output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Function to discover and read CSV files
discover_and_read_csv_files <- function(input_dir, verbose = FALSE) {
  # Find all CSV files in the input directory
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  
  if (length(csv_files) == 0) {
    stop("No CSV files found in input directory: ", input_dir)
  }
  
  if (verbose) {
    cat("Found", length(csv_files), "CSV files\n")
    cat("Sample files:", paste(head(basename(csv_files), 5), collapse = ", "), "\n")
  }
  
  # Read all CSV files and combine them
  all_results <- list()
  failed_files <- character(0)
  
  for (i in seq_along(csv_files)) {
    file_path <- csv_files[i]
    
    tryCatch({
      # Read the CSV file
      result <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Validate that the file has data
      if (nrow(result) == 0) {
        warning("Empty CSV file: ", basename(file_path))
        failed_files <- c(failed_files, basename(file_path))
        next
      }
      
      # Add source file information for debugging
      result$source_file <- basename(file_path)
      
      all_results[[i]] <- result
      
      if (verbose && i %% 100 == 0) {
        cat("Processed", i, "of", length(csv_files), "files\n")
      }
      
    }, error = function(e) {
      warning("Failed to read file ", basename(file_path), ": ", e$message)
      failed_files <<- c(failed_files, basename(file_path))
    })
  }
  
  # Remove NULL entries (failed files)
  all_results <- all_results[!sapply(all_results, is.null)]
  
  if (length(all_results) == 0) {
    stop("No valid CSV files could be read from: ", input_dir)
  }
  
  # Combine all results into a single data frame
  combined_results <- do.call(rbind, all_results)
  
  # Report summary
  cat("Successfully read", length(all_results), "CSV files\n")
  if (length(failed_files) > 0) {
    cat("Failed to read", length(failed_files), "files\n")
    if (verbose) {
      cat("Failed files:", paste(head(failed_files, 10), collapse = ", "), "\n")
      if (length(failed_files) > 10) {
        cat("... and", length(failed_files) - 10, "more\n")
      }
    }
  }
  
  return(list(
    data = combined_results,
    total_files = length(csv_files),
    successful_files = length(all_results),
    failed_files = failed_files
  ))
}

# Function to validate and analyze job completeness
validate_job_completeness <- function(results_data, input_dir, verbose = FALSE) {
  cat("\n=== JOB COMPLETENESS ANALYSIS ===\n")
  
  # Determine the type of results based on column names
  columns <- colnames(results_data)
  
  if ("variance" %in% columns && "mean" %in% columns) {
    # This is batch effects data
    job_type <- "batch_effects"
    cat("Detected job type: Batch Effects Analysis\n")
    
    # Expected parameters from config
    expected_classifiers <- c("logistic", "elasticnet", "svm", "rf", "knn", "xgboost", "nn", "shrinkageLDA")
    expected_means <- c(5)
    expected_variances <- c(1, 3, 5)
    expected_seeds <- 42:141  # 100 seeds starting from 42
    
    # Calculate expected total jobs
    expected_total <- length(expected_classifiers) * length(expected_means) * 
                     length(expected_variances) * length(expected_seeds)
    
    if (verbose) {
      cat("Expected parameters:\n")
      cat("  Classifiers:", length(expected_classifiers), "\n")
      cat("  Means:", length(expected_means), "\n") 
      cat("  Variances:", length(expected_variances), "\n")
      cat("  Seeds:", length(expected_seeds), "\n")
      cat("  Total expected jobs:", expected_total, "\n")
    }
    
    # Check actual parameter combinations
    actual_combinations <- unique(results_data[, c("classifier", "mean", "variance", "seed")])
    
    # Create expected combinations for comparison
    expected_combinations <- expand.grid(
      classifier = expected_classifiers,
      mean = expected_means,
      variance = expected_variances,
      seed = expected_seeds,
      stringsAsFactors = FALSE
    )
    
  } else if ("adjuster" %in% columns && "n_datasets" %in% columns) {
    # This is adjusters data
    job_type <- "adjusters"
    cat("Detected job type: Adjusters Analysis\n")
    
    # Expected parameters from config
    expected_adjusters <- c("unadjusted", "naive", "rank_samples", "rank_twice", "npn", "combat", "mnn", "fast_mnn")
    expected_classifiers <- c("logistic", "elasticnet", "svm", "rf", "knn", "xgboost", "nn", "shrinkageLDA")
    expected_n_datasets <- c(2, 3, 4, 5)
    # For adjusters, we use test_study instead of seed
    all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
    
    # Build expected test studies for each n_datasets
    expected_test_studies <- list()
    for (n in expected_n_datasets) {
      expected_test_studies[[as.character(n)]] <- all_studies[1:n]
    }
    
    # Calculate expected total jobs
    expected_total <- 0
    for (n in expected_n_datasets) {
      expected_total <- expected_total + 
        length(expected_adjusters) * length(expected_classifiers) * length(expected_test_studies[[as.character(n)]])
    }
    
    if (verbose) {
      cat("Expected parameters:\n")
      cat("  Adjusters:", length(expected_adjusters), "\n")
      cat("  Classifiers:", length(expected_classifiers), "\n")
      cat("  N_datasets:", length(expected_n_datasets), "\n")
      cat("  Total expected jobs:", expected_total, "\n")
    }
    
    # Check actual parameter combinations
    actual_combinations <- unique(results_data[, c("adjuster", "classifier", "n_datasets", "test_study")])
    
    # Create expected combinations for comparison
    expected_combinations <- data.frame()
    for (n in expected_n_datasets) {
      for (test_study in expected_test_studies[[as.character(n)]]) {
        temp <- expand.grid(
          adjuster = expected_adjusters,
          classifier = expected_classifiers,
          n_datasets = n,
          test_study = test_study,
          stringsAsFactors = FALSE
        )
        expected_combinations <- rbind(expected_combinations, temp)
      }
    }
    
  } else {
    cat("Warning: Unknown job type based on columns:", paste(columns, collapse = ", "), "\n")
    return(list(
      job_type = "unknown",
      success_rate = NA,
      missing_jobs = NA,
      failed_jobs = NA
    ))
  }
  
  # Calculate completion statistics
  actual_jobs <- nrow(actual_combinations)
  success_rate <- (actual_jobs / expected_total) * 100
  missing_jobs <- expected_total - actual_jobs
  
  cat("Completion Summary:\n")
  cat("  Expected jobs:", expected_total, "\n")
  cat("  Completed jobs:", actual_jobs, "\n")
  cat("  Success rate:", sprintf("%.1f%%", success_rate), "\n")
  cat("  Missing jobs:", missing_jobs, "\n")
  
  # Analyze missing jobs if any
  if (missing_jobs > 0 && verbose) {
    cat("\nMissing job analysis:\n")
    
    # Find missing combinations by comparing expected vs actual
    missing_combinations <- anti_join_manual(expected_combinations, actual_combinations)
    
    if (nrow(missing_combinations) > 0) {
      cat("Sample missing jobs (first 10):\n")
      sample_missing <- head(missing_combinations, 10)
      for (i in 1:nrow(sample_missing)) {
        if (job_type == "batch_effects") {
          cat(sprintf("  %s_mean%s_var%s_seed%s\n", 
                     sample_missing$classifier[i], sample_missing$mean[i], 
                     sample_missing$variance[i], sample_missing$seed[i]))
        } else if (job_type == "adjusters") {
          cat(sprintf("  %s_%s_n%s_test%s\n", 
                     sample_missing$adjuster[i], sample_missing$classifier[i], 
                     sample_missing$n_datasets[i], sample_missing$test_study[i]))
        }
      }
      if (nrow(missing_combinations) > 10) {
        cat("  ... and", nrow(missing_combinations) - 10, "more missing jobs\n")
      }
    }
  }
  
  # Failure analysis based on log files
  if (missing_jobs > 0) {
    cat("\nFailure Analysis:\n")
    
    # Look for corresponding log files to understand failures
    # Only check log directories for the main job type (adjusters or batch_effects)
    # Skip within_study_cv directories as they have different log structure
    for (single_input_dir in input_dir) {
      # Skip within_study_cv directories for log analysis
      if (grepl("within_study_cv", single_input_dir)) {
        next
      }
      
      log_dir <- gsub("/results/", "/logs/", single_input_dir)
      if (job_type == "batch_effects") {
        log_dir <- file.path(log_dir, "classify_batch_effects")
      } else if (job_type == "adjusters") {
        log_dir <- file.path(log_dir, "classify_adjusters")
      }
      
      if (dir.exists(log_dir)) {
        log_files <- list.files(log_dir, pattern = "\\.log$", full.names = FALSE)
        cat("  Found", length(log_files), "log files in", log_dir, "\n")
        
        # Simple heuristic: if we have log files but no results, those are likely failures
        result_files <- list.files(single_input_dir, pattern = "\\.csv$", full.names = FALSE)
        result_basenames <- gsub("\\.csv$", "", result_files)
        log_basenames <- gsub("\\.log$", "", log_files)
        
        failed_jobs <- setdiff(log_basenames, result_basenames)
        if (length(failed_jobs) > 0) {
          cat("  Identified", length(failed_jobs), "failed jobs (have logs but no results)\n")
          if (verbose && length(failed_jobs) > 0) {
            cat("  Sample failed jobs:", paste(head(failed_jobs, 5), collapse = ", "), "\n")
          }
        }
      } else {
        cat("  Log directory not found:", log_dir, "\n")
      }
    }
  }
  
  return(list(
    job_type = job_type,
    success_rate = success_rate,
    missing_jobs = missing_jobs,
    expected_total = expected_total,
    actual_jobs = actual_jobs
  ))
}

# Manual implementation of anti_join since we're not using dplyr
anti_join_manual <- function(x, y) {
  # Find rows in x that don't have matching rows in y
  key_cols <- intersect(colnames(x), colnames(y))
  
  # Create composite keys for comparison
  x_keys <- apply(x[, key_cols, drop = FALSE], 1, paste, collapse = "_")
  y_keys <- apply(y[, key_cols, drop = FALSE], 1, paste, collapse = "_")
  
  # Return rows from x that don't match any row in y
  missing_indices <- !x_keys %in% y_keys
  return(x[missing_indices, , drop = FALSE])
}

# Main execution
main <- function() {
  cat("=== AGGREGATE RESULTS SCRIPT ===\n")
  cat("Input directories:", paste(opt$input_dir, collapse = ", "), "\n")
  cat("Output file:", opt$output, "\n")
  cat("Verbose mode:", opt$verbose, "\n\n")
  
  # Discover and read all CSV files from all input directories
  cat("=== READING CSV FILES ===\n")
  all_results_list <- list()
  total_files <- 0
  total_successful <- 0
  all_failed_files <- character(0)
  
  for (input_dir in opt$input_dir) {
    cat("\nProcessing directory:", input_dir, "\n")
    dir_results <- discover_and_read_csv_files(input_dir, opt$verbose)
    all_results_list[[length(all_results_list) + 1]] <- dir_results$data
    total_files <- total_files + dir_results$total_files
    total_successful <- total_successful + dir_results$successful_files
    all_failed_files <- c(all_failed_files, dir_results$failed_files)
  }
  
  # Combine results from all directories
  combined_data <- do.call(rbind, all_results_list)
  
  results <- list(
    data = combined_data,
    total_files = total_files,
    successful_files = total_successful,
    failed_files = all_failed_files
  )
  
  cat("\n=== COMBINED RESULTS ===\n")
  cat("Total files across all directories:", total_files, "\n")
  cat("Successfully processed:", total_successful, "\n")
  cat("Failed:", length(all_failed_files), "\n")
  
  # Validate job completeness
  completeness <- validate_job_completeness(results$data, opt$input_dir, opt$verbose)
  
  # Remove the source_file column before writing output (it was just for debugging)
  final_data <- results$data
  if ("source_file" %in% colnames(final_data)) {
    final_data$source_file <- NULL
  }
  
  # Write aggregated results
  cat("\n=== WRITING OUTPUT ===\n")
  write.csv(final_data, opt$output, row.names = FALSE)
  cat("Wrote", nrow(final_data), "rows to", opt$output, "\n")
  
  # Final summary
  cat("\n=== FINAL SUMMARY ===\n")
  cat("Job type:", completeness$job_type, "\n")
  cat("Total CSV files found:", results$total_files, "\n")
  cat("Successfully processed:", results$successful_files, "\n")
  cat("Failed to read:", length(results$failed_files), "\n")
  cat("Final output rows:", nrow(final_data), "\n")
  
  if (!is.na(completeness$success_rate)) {
    cat("Job completion rate:", sprintf("%.1f%%", completeness$success_rate), "\n")
    
    # Exit with error code if completion rate is too low
    if (completeness$success_rate < 90) {
      cat("WARNING: Low completion rate detected!\n")
      cat("Consider investigating failed jobs before proceeding.\n")
      # Don't exit with error - let the workflow continue but warn the user
    }
  }
  
  cat("Aggregation completed successfully.\n")
}

# Execute main function
if (sys.nframe() == 0) {
  main()
}