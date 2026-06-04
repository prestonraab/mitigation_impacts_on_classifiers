#!/usr/bin/env Rscript

# aggregate_class_imbalanced_results.R - Aggregate class imbalanced classification results

suppressMessages(suppressWarnings({
  required_packages <- c("argparse", "dplyr")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Aggregate class imbalanced classification results")

parser$add_argument("--input-dir", type = "character", required = TRUE,
                   help = "Directory containing individual result CSV files")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output aggregated CSV file path")

args <- parser$parse_args()

# ====================================================================
# MAIN AGGREGATION FUNCTION
# ====================================================================

aggregate_class_imbalanced_results <- function() {
  
  cat("=== CLASS IMBALANCED RESULTS AGGREGATION ===\n")
  cat(sprintf("Input directory: %s\n", args$input_dir))
  cat(sprintf("Output file: %s\n", args$output))
  cat("============================================\n\n")
  
  # Check if input directory exists
  if (!dir.exists(args$input_dir)) {
    stop(sprintf("Input directory not found: %s", args$input_dir))
  }
  
  # Find all CSV files in the input directory
  csv_files <- list.files(args$input_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  
  if (length(csv_files) == 0) {
    stop(sprintf("No CSV files found in directory: %s", args$input_dir))
  }
  
  cat(sprintf("Found %d CSV files to aggregate\n", length(csv_files)))
  
  # Initialize list to store dataframes
  all_results <- list()
  
  # Read and combine all CSV files
  for (i in seq_along(csv_files)) {
    file_path <- csv_files[i]
    file_name <- basename(file_path)
    
    tryCatch({
      # Read the CSV file
      df <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Validate that it has the expected columns
      required_cols <- c("scenario_name", "adjuster", "classifier", "mcc", "imbalance_pct", "training_pair")
      missing_cols <- setdiff(required_cols, colnames(df))
      
      if (length(missing_cols) > 0) {
        cat(sprintf("Warning: File %s missing columns: %s\n", file_name, paste(missing_cols, collapse = ", ")))
        next
      }
      
      # Add file source information
      df$source_file <- file_name
      
      # Add to results list
      all_results[[i]] <- df
      
      if (i %% 50 == 0) {
        cat(sprintf("Processed %d/%d files\n", i, length(csv_files)))
      }
      
    }, error = function(e) {
      cat(sprintf("Error reading file %s: %s\n", file_name, e$message))
    })
  }
  
  # Remove NULL entries (failed reads)
  all_results <- all_results[!sapply(all_results, is.null)]
  
  if (length(all_results) == 0) {
    stop("No valid CSV files could be read")
  }
  
  cat(sprintf("Successfully read %d files\n", length(all_results)))
  
  # Combine all dataframes
  combined_results <- do.call(rbind, all_results)
  
  cat(sprintf("Combined results: %d rows, %d columns\n", nrow(combined_results), ncol(combined_results)))
  
  # Add some derived columns for analysis
  combined_results <- combined_results %>%
    mutate(
      # Convert imbalance percentage to factor for easier plotting
      imbalance_level = factor(sprintf("%.0f%%", imbalance_pct * 100), 
                              levels = c("20%", "30%", "40%", "50%")),
      
      # Create a unique trial identifier
      trial_id = paste(scenario_name, adjuster, classifier, sep = "_"),
      
      # Extract training pair info if not already present
      training_pair = ifelse("training_pair" %in% colnames(combined_results), 
                            training_pair, 
                            gsub("_imbal.*", "", scenario_name)),
      
      # Create grouping variables for analysis
      classifier_group = case_when(
        classifier %in% c("logistic", "elasticnet") ~ "Linear",
        classifier %in% c("rf", "xgboost") ~ "Tree-based", 
        classifier %in% c("svm", "nnet") ~ "Non-linear",
        classifier %in% c("knn", "shrinkageLDA") ~ "Distance-based",
        TRUE ~ "Other"
      )
    )
  
  # Print summary statistics
  cat("\nSummary statistics:\n")
  cat(sprintf("  Unique scenarios: %d\n", length(unique(combined_results$scenario_name))))
  cat(sprintf("  Adjusters: %s\n", paste(unique(combined_results$adjuster), collapse = ", ")))
  cat(sprintf("  Classifiers: %s\n", paste(unique(combined_results$classifier), collapse = ", ")))
  cat(sprintf("  Imbalance levels: %s\n", paste(unique(combined_results$imbalance_level), collapse = ", ")))
  cat(sprintf("  Dataset combinations: %d\n", length(unique(combined_results$training_pair))))
  
  # Check for missing values in key metrics
  missing_mcc <- sum(is.na(combined_results$mcc))
  if (missing_mcc > 0) {
    cat(sprintf("  Warning: %d rows have missing MCC values\n", missing_mcc))
  }
  
  # Print performance summary by adjuster
  cat("\nPerformance summary by adjuster (mean MCC ± SD):\n")
  perf_summary <- combined_results %>%
    group_by(adjuster) %>%
    summarise(
      n_trials = n(),
      mean_mcc = mean(mcc, na.rm = TRUE),
      sd_mcc = sd(mcc, na.rm = TRUE),
      median_mcc = median(mcc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_mcc))
  
  print(perf_summary)
  
  # Print performance summary by imbalance level
  cat("\nPerformance summary by imbalance level (mean MCC ± SD):\n")
  imbal_summary <- combined_results %>%
    group_by(imbalance_level) %>%
    summarise(
      n_trials = n(),
      mean_mcc = mean(mcc, na.rm = TRUE),
      sd_mcc = sd(mcc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(imbalance_level)
  
  print(imbal_summary)
  
  # Create output directory if needed
  output_dir <- dirname(args$output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save aggregated results
  write.csv(combined_results, args$output, row.names = FALSE)
  
  cat(sprintf("\n✅ Aggregated results saved to: %s\n", args$output))
  
  return(combined_results)
}

# ====================================================================
# EXECUTE WITH ERROR HANDLING
# ====================================================================

tryCatch({
  result <- aggregate_class_imbalanced_results()
}, error = function(e) {
  cat(sprintf("[ERROR] %s\n", e$message), file = stderr())
  quit(status = 1)
})