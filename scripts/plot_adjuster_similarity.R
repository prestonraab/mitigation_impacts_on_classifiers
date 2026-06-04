#!/usr/bin/env Rscript

# plot_adjuster_similarity.R
# Generate an Adjuster x Adjuster similarity heatmap based on Euclidean distance of rank_diff embeddings.

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(seriation)
  library(ragg)
  library(circlize)
})

# Load utils if present
if (file.exists("scripts/adjuster_plot_utils.R")) {
  source("scripts/adjuster_plot_utils.R")
}

# Hodges-Lehmann estimator
calc_hodges_lehmann <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA)
  return(median(outer(x, x, "+") / 2))
}

parser <- ArgumentParser(description = "Generate Adjuster x Adjuster similarity heatmap")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")

args <- parser$parse_args()

# ==============================================================================
# Constants & Helpers
# ==============================================================================

target_adjusters <- if (!is.null(args$adjusters)) trimws(strsplit(args$adjusters, ",")[[1]]) else NULL

# ==============================================================================
# Data Processing (Same as GAPR heatmap)
# ==============================================================================
# Load and clean raw data
raw_data <- read.csv(args$input, stringsAsFactors = FALSE) %>%
  mutate(across(c(adjuster, classifier), trimws))

# Apply adjuster filter if specified
if (!is.null(target_adjusters)) {
  raw_data <- raw_data %>% filter(adjuster %in% target_adjusters)
}

raw_data <- raw_data %>%
  filter(
    metric == "mcc",
    !is.na(value),
    adjuster != "within_study_cv",
    !classifier %in% c("logistic", "rda")
  )

ranked_data <- raw_data %>%
  mutate(
    rank = rank(-value, ties.method = "average"),
    .by = c(classifier, n_datasets, test_study)
  ) %>%
  mutate(
    adjuster_label = if (exists("format_adjuster_label")) {
      sapply(adjuster, format_adjuster_label)
    } else {
      adjuster
    }
  )

classifier_stats <- ranked_data %>%
  summarise(
    avg_rank = calc_hodges_lehmann(rank),
    .by = c(classifier, adjuster_label)
  )

group_centers <- ranked_data %>%
  summarise(
    center = calc_hodges_lehmann(rank),
    .by = adjuster_label
  )

heatmap_data <- classifier_stats %>%
  left_join(group_centers, by = "adjuster_label") %>%
  mutate(rank_diff = avg_rank - center) %>%
  mutate(rank_diff = ifelse(is.na(rank_diff), 0, rank_diff)) %>%
  select(adjuster_label, classifier, rank_diff) %>%
  pivot_wider(names_from = classifier, values_from = rank_diff)

# Embedding matrix: Adjusters (rows) x Classifiers (columns)
mat <- as.matrix(heatmap_data[,-1])
rownames(mat) <- heatmap_data$adjuster_label
mat[is.na(mat)] <- 0

# ==============================================================================
# Visualization using GAPR::GAP()
# ==============================================================================

# For the Adjuster Similarity heatmap, we provide the proximity (distance) matrix directly
# and set isProximityMatrix = TRUE.

# GAPR expects the distance matrix to be converted to a matrix if it's a dist object
dist_mat <- as.matrix(dist(mat, method = "euclidean"))

# Workaround for potential internal scoping issues in GAPR palette handling
if (!exists("num_colorspectrum")) {
  assign("num_colorspectrum", 100, envir = .GlobalEnv)
}

# The user wants to cluster/arrange adjusters and plot similarities
# Using R2E order and flip on the proximity matrix
similarity_result <- GAPR::GAP(
  data = dist_mat,
  isProximityMatrix = TRUE,
  row.order = 'r2e',
  row.flip = 'r2e',
  original.color = 'Greys', # Sequential palette for strictly non-negative distances
  row.color = 'Greys',
  col.color = 'Greys',
  colorbar.margin = 0.5,
  row.label.size = 12, # Increased for better visibility
  col.label.size = 12, # Increased for better visibility
  border = TRUE,
  border.width = 1,
  exp.row_order = TRUE,
  exp.row_names = TRUE,
  exp.originalmatrix = TRUE,
  exp.row_prox = TRUE,
  PNGfilename = args$output,
  PNGwidth = 3600,
  PNGheight = 2400,
  PNGres = 300,
  show.plot = FALSE
)

cat("Saved adjuster similarity heatmap to:", args$output, "\n")

cat("Saved similarity heatmap to:", args$output, "\n")
