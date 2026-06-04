#!/usr/bin/env Rscript

# plot_adjusters_main_only.R
# Generate only the main adjuster effectiveness plot

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(scales)
})

source("scripts/adjuster_plot_utils.R")
source("scripts/generate_main_plot.R")

parser <- ArgumentParser(description = "Create main adjuster effectiveness plot")
parser$add_argument("-i", "--input", type = "character", required = TRUE)
parser$add_argument("-o", "--output", type = "character", required = TRUE)
parser$add_argument("--width", type = "double", default = 20)
parser$add_argument("--height", type = "double", default = 16)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL)

args <- parser$parse_args()

# Load data
data <- read.csv(args$input, stringsAsFactors = FALSE)
if (!is.null(args$adjusters)) {
  adjusters_filter <- trimws(strsplit(args$adjusters, ",")[[1]])
  data <- data[data$adjuster %in% adjusters_filter, ]
}

mxe_data <- data[data$metric == "mcc" & !is.na(data$n_datasets), ]
mxe_data$classifier_label <- factor(mxe_data$classifier,
  levels = c("logistic", "elasticnet", "svm", "rf", "knn", "xgboost", "nnet", "shrinkageLDA"),
  labels = c("Logistic", "ElasticNet", "SVM", "Random Forest", "KNN", "XGBoost", "Neural Net", "RDA"))

# Calculate ordering
adjuster_performance <- mxe_data %>%
  group_by(adjuster) %>%
  summarise(overall_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean))

unique_adjusters <- adjuster_performance$adjuster
adjuster_labels <- sapply(unique_adjusters, format_adjuster_label)

mxe_data$adjuster_label <- factor(mxe_data$adjuster,
  levels = unique_adjusters,
  labels = adjuster_labels)

# Generate plot
generate_main_plot(mxe_data, args$output, args$width, args$height, args$dpi, TRUE)
