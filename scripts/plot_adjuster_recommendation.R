#!/usr/bin/env Rscript

# plot_adjuster_recommendation.R
# Plot the best adjusters for each specific classifier to answer:
# "I've chosen my classifier, which adjuster should I use?"

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(ggbeeswarm)
})

# Load utils if present
if (file.exists("scripts/adjuster_plot_utils.R")) {
  source("scripts/adjuster_plot_utils.R")
}

parser <- ArgumentParser(description = "Plot best adjusters per classifier")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 15)
parser$add_argument("--height", type = "double", default = 7.5)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")
parser$add_argument("--top-n", type = "integer", default = 10,
                    help = "Number of top adjusters to show per classifier")

args <- parser$parse_args()

# ==============================================================================
# Constants & Helpers
# ==============================================================================

classifier_labels_map <- c(
  "rda" = "RDA", "logistic" = "Logistic Regression", "elasticnet" = "ElasticNet",
  "svm" = "SVM", "rf" = "Random Forest", "knn" = "KNN",
  "xgboost" = "XGBoost", "nnet" = "Neural Net", "shrinkageLDA" = "RDA"
)

# Calculate Hodges-Lehmann estimator (Pseudomedian)
calc_hodges_lehmann <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA)
  return(median(outer(x, x, "+") / 2))
}

# ==============================================================================
# Data Processing
# ==============================================================================

target_adjusters <- if (!is.null(args$adjusters)) trimws(strsplit(args$adjusters, ",")[[1]]) else NULL

raw_data <- read.csv(args$input, stringsAsFactors = FALSE) %>%
  mutate(across(c(adjuster, classifier), trimws)) %>%
  filter(
    metric == "mcc",
    !is.na(value),
    !classifier %in% c("logistic", "rda"),
    is.null(target_adjusters) | adjuster %in% target_adjusters | adjuster == "within_study_cv"
  )

# Calculate MCC Delta: MCC_adjuster - MCC_within_study_cv
cv_baseline <- raw_data %>%
  filter(adjuster == "within_study_cv") %>%
  select(classifier, n_datasets, test_study, cv_mcc = value)

delta_data <- raw_data %>%
  filter(adjuster != "within_study_cv") %>%
  left_join(cv_baseline, by = c("classifier", "n_datasets", "test_study")) %>%
  mutate(mcc_delta = value - cv_mcc) %>%
  filter(!is.na(mcc_delta)) %>%
  mutate(
    classifier_label = recode(classifier, !!!classifier_labels_map),
    adjuster_label = if (exists("format_adjuster_label")) {
      sapply(adjuster, format_adjuster_label)
    } else {
      adjuster
    }
  )

# 1. Calculate average performance for all adjuster-classifier pairs
all_stats <- delta_data %>%
  summarise(
    avg_delta = calc_hodges_lehmann(mcc_delta),
    .by = c(classifier_label, adjuster_label)
  )

# 2. Identify the Top 4 adjusters for each classifier
top4_stats <- all_stats %>%
  group_by(classifier_label) %>%
  slice_max(avg_delta, n = 4, with_ties = FALSE) %>%
  mutate(selection_type = "Top 4") %>%
  ungroup()

# 3. Create a global pool of adjusters that are "Top 4" for at least one classifier
global_pool <- unique(top4_stats$adjuster_label)

# 4. For each classifier, find the worst-performing adjuster from this global pool
bottom_of_pool_stats <- all_stats %>%
  filter(adjuster_label %in% global_pool) %>%
  group_by(classifier_label) %>%
  slice_min(avg_delta, n = 1, with_ties = FALSE) %>%
  mutate(selection_type = "Bottom 1 of Top Pool") %>%
  ungroup()

# 5. Combine Top 4 and the "worst of pool" for the final selection
selected_stats <- bind_rows(top4_stats, bottom_of_pool_stats) %>%
  distinct(classifier_label, adjuster_label, .keep_all = TRUE)

# Filter main data to include only selected adjusters
plot_data <- delta_data %>%
  inner_join(select(selected_stats, classifier_label, adjuster_label, avg_delta, selection_type), 
             by = c("classifier_label", "adjuster_label"))

# Create composite factor for per-facet ordering
plot_data <- plot_data %>%
  mutate(
    facet_item = paste0(classifier_label, "___", adjuster_label),
    facet_item = reorder(facet_item, avg_delta) # Higher delta (better) is higher on Y axis
  )

# ==============================================================================
# Visualization
# ==============================================================================

# Shared beeswarm position logic
pos_swarm_h <- position_beeswarm(method = "hex", cex = 1.0, groupOnX = FALSE)

p <- ggplot(plot_data, aes(x = mcc_delta, y = facet_item)) +

  # Background distribution
  geom_violin(
    aes(fill = adjuster_label),
    width = 0.7,
    alpha = 0.1,
    color = "grey80",
    trim = TRUE,
    show.legend = FALSE
  ) +

  # Reference line at 0 (Baseline)
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", size = 0.5) +

  # HL center lines
  stat_summary(
    fun = calc_hodges_lehmann,
    fun.min = calc_hodges_lehmann,
    fun.max = calc_hodges_lehmann,
    geom = "crossbar",
    width = 0.5,
    size = 0.4, 
    color = "grey40"
  ) +

  # HL ranks per study
  geom_point(
    aes(color = adjuster_label),
    position = pos_swarm_h,
    size = 2.0,
    alpha = 0.6,
    show.legend = FALSE
  ) +

  # Faceting
  facet_wrap(~classifier_label, scales = "free_y", ncol = 3) +

  # Scales and theme
  scale_y_discrete(labels = function(x) gsub(".*___", "", x)) +
  scale_x_reverse() +
  
  labs(
    title = NULL,
    subtitle = NULL,
    y = NULL,
    x = "MCC Adjuster - MCC Within-Study CV"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.major.x = element_line(color = "gray92", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5, "mm"),
    panel.spacing = unit(1, "lines")
  )

# Save output
ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved recommendation plot to:", args$output, "\n")
