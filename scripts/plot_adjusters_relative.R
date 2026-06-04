#!/usr/bin/env Rscript
# Standalone script to generate aggregated relative performance plot

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Generate aggregated relative performance plot")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                   help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output PNG file for aggregated relative plot")
parser$add_argument("--adjusters", type = "character", default = NULL,
                   help = "Comma-separated list of adjusters to include")

opt <- parser$parse_args()

# Source utility functions
source("scripts/adjuster_plot_utils.R")
source("scripts/generate_relative_plot_aggregated.R")

# Read data
cat("Reading data from:", opt$input, "\n")
data <- read.csv(opt$input, stringsAsFactors = FALSE)

# Filter adjusters if specified
if (!is.null(opt$adjusters)) {
  selected_adjusters <- trimws(strsplit(opt$adjusters, ",")[[1]])
  # Always include within_study_cv as it's needed as baseline reference
  selected_adjusters <- unique(c(selected_adjusters, "within_study_cv"))
  data <- data[data$adjuster %in% selected_adjusters, ]
  cat("Filtered to", length(selected_adjusters), "adjusters\n")
}

# Filter to MCC metric and add classifier labels
mxe_data <- data[data$metric == "mcc" & !is.na(data$n_datasets), ]
mxe_data$classifier_label <- factor(mxe_data$classifier,
  levels = c("logistic", "elasticnet", "svm", "rf", "knn", "xgboost", "nnet", "shrinkageLDA"),
  labels = c("Logistic", "ElasticNet", "SVM", "Random Forest", "KNN", "XGBoost", "Neural Net", "RDA"))

# Generate aggregated relative plot
generate_relative_plot_aggregated(
  mxe_data = mxe_data,
  top_adjuster = NULL,  # Not used
  top_adjuster_label = NULL,  # Not used
  output_file = opt$output,
  width = 20,
  height = 16,
  dpi = 300
)

cat("Done!\n")
