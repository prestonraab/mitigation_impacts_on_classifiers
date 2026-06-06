#!/usr/bin/env Rscript

# plot_average_rank_by_classifier.R
# Merged version: Uses Snippet 1 color logic + Snippet 2 manual boxplot layout.

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
})

# Load utils if present (Crucial for correct label formatting)
if (file.exists("scripts/adjuster_plot_utils.R")) {
  source("scripts/adjuster_plot_utils.R")
}

parser <- ArgumentParser(description = "Create MCC violin plot by classifier with adjuster means")
parser$add_argument("-i", "--input", type = "character", required = TRUE, help = "Input CSV")
parser$add_argument("-o", "--output", type = "character", required = TRUE, help = "Output PNG")
parser$add_argument("--width", type = "double", default = 8)
parser$add_argument("--height", type = "double", default = 14)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL, help = "Filter adjusters")
parser$add_argument("--n-datasets", type = "character", default = "4", help = "Filter n_datasets")

args <- parser$parse_args()

# Load Data
data <- read.csv(args$input, stringsAsFactors = FALSE)
data <- data %>%
  filter(metric == "mcc", !is.na(value), !is.na(n_datasets), !is.na(classifier))

# Clean whitespace to prevent matching errors
data$adjuster <- trimws(data$adjuster)
data$classifier <- trimws(data$classifier)

# Exclude within-study CV (not a cross-study method)
data <- data %>% filter(adjuster != "within_study_cv")

# Filter N-Datasets
n_datasets_label <- args$n_datasets
if (tolower(args$n_datasets) != "all") {
  n_val <- as.integer(args$n_datasets)
  data <- data %>% filter(n_datasets == n_val)
  n_datasets_label <- sprintf("%d-Study", n_val)
} else {
  n_datasets_label <- "All Studies"
}

# Filter Adjusters
if (!is.null(args$adjusters)) {
  adjusters_filter <- trimws(strsplit(args$adjusters, ",")[[1]])
  data <- data %>% filter(adjuster %in% adjusters_filter)
}

# Label Classifiers
classifier_labels <- c(
  "logistic" = "Logistic\nRegression", "elasticnet" = "ElasticNet", "svm" = "SVM",
  "rf" = "Random\nForest", "knn" = "KNN", "xgboost" = "XGBoost",
  "nnet" = "Neural\nNet", "shrinkageLDA" = "RDA"
)
data$classifier_label <- classifier_labels[data$classifier]

# Format Adjuster Labels (Uses utils logic if available)
if (exists("format_adjuster_label")) {
    data$adjuster_label <- sapply(data$adjuster, format_adjuster_label)
} else {
    data$adjuster_label <- data$adjuster
}

# Order Adjusters (Best Median on Top)
adjuster_order <- data %>%
  group_by(adjuster_label) %>%
  summarise(median_mcc = median(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_mcc))

data$adjuster_label <- factor(
  data$adjuster_label,
  levels = adjuster_order$adjuster_label,
  ordered = TRUE
)

# Order Classifiers
classifier_order <- data %>%
  group_by(classifier, classifier_label) %>%
  summarise(median_mcc = median(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_mcc))

data$classifier_label <- factor(
  data$classifier_label,
  levels = classifier_order$classifier_label,
  ordered = TRUE
)
data <- data %>% filter(!is.na(classifier_label))

# ==============================================================================
# Compute Adjuster Colors via K-Means on Interaction Embeddings
# (Mirrors plot_gapr_heatmap.R preprocessing)
# ==============================================================================

compute_interaction_embeddings <- function(data_file) {
  classifier_labels_map <- c(
    "rda" = "RDA", "logistic" = "LR", "elasticnet" = "ENet",
    "svm" = "SVM", "rf" = "RF", "knn" = "KNN",
    "xgboost" = "XGB", "nnet" = "NN", "shrinkageLDA" = "RDA"
  )
  calc_hl <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    median(outer(x, x, "+") / 2)
  }
  raw <- read.csv(data_file, stringsAsFactors = FALSE) %>%
    mutate(across(c(adjuster, classifier), trimws)) %>%
    filter(metric == "mcc", !is.na(value),
           adjuster != "within_study_cv",
           !classifier %in% "logistic")
  raw$adjuster_label <- if (exists("format_adjuster_label")) {
    sapply(raw$adjuster, format_adjuster_label)
  } else {
    raw$adjuster
  }
  raw$classifier_label <- ifelse(
    raw$classifier %in% names(classifier_labels_map),
    classifier_labels_map[raw$classifier], raw$classifier
  )
  centers <- raw %>%
    group_by(adjuster_label) %>%
    summarise(center = calc_hl(value), .groups = "drop")
  wide <- raw %>%
    group_by(adjuster_label, classifier_label) %>%
    summarise(avg_mcc = mean(value), .groups = "drop") %>%
    left_join(centers, by = "adjuster_label") %>%
    mutate(mcc_diff = avg_mcc - center) %>%
    select(adjuster_label, classifier_label, mcc_diff) %>%
    pivot_wider(names_from = classifier_label, values_from = mcc_diff, values_fill = 0)
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$adjuster_label
  mat[is.na(mat)] <- 0
  mat
}

mat_embed <- compute_interaction_embeddings(args$input)

# Choose k via average silhouette width (k = 2..6)
best_k <- 4L
if (nrow(mat_embed) >= 4L && requireNamespace("cluster", quietly = TRUE)) {
  k_range <- seq(2L, min(6L, nrow(mat_embed) - 1L))
  sil_scores <- vapply(k_range, function(k) {
    set.seed(42L)
    km  <- kmeans(mat_embed, centers = k, nstart = 25L)
    sil <- cluster::silhouette(km$cluster, dist(mat_embed))
    mean(sil[, 3])
  }, numeric(1))
  best_k <- k_range[which.max(sil_scores)]
}
set.seed(42L)
km_result      <- kmeans(mat_embed, centers = best_k, nstart = 25L)
adjuster_cluster <- km_result$cluster  # named integer vector, names = adjuster_label

cluster_colors <- RColorBrewer::brewer.pal(max(3L, best_k), "Dark2")[seq_len(best_k)]

all_adjusters <- levels(data$adjuster_label)
palette_map <- setNames(
  vapply(all_adjusters, function(lbl) {
    if (lbl %in% names(adjuster_cluster)) cluster_colors[adjuster_cluster[[lbl]]] else "#F0E442"
  }, character(1)),
  all_adjusters
)

# ==============================================================================
# Manual Boxplot Calculation (Logic from Snippet 2)
# ==============================================================================

box_stats <- data %>%
  group_by(classifier_label) %>%
  summarise(
    q1 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    lower_whisker = max(min(value), q1 - 1.5 * iqr),
    upper_whisker = min(max(value), q3 + 1.5 * iqr),
    
    # Calculate scale factor for positioning below stack
    n_adjusters = n_distinct(adjuster_label),
    base_density = max(density(value, na.rm = TRUE)$y),
    scale_factor = base_density * n_adjusters
  ) %>%
  mutate(
    box_top    = -0.05 * scale_factor,
    box_bottom = -0.20 * scale_factor,
    mid_y      = (box_top + box_bottom) / 2
  )

# Align factors for join
box_stats$classifier_label <- factor(box_stats$classifier_label, levels = levels(data$classifier_label))

# Identify Outliers
outliers <- data %>%
  left_join(box_stats, by = "classifier_label") %>%
  filter(value < lower_whisker | value > upper_whisker) %>%
  mutate(y_pos = mid_y)

# ==============================================================================
# Plot
# ==============================================================================

p <- ggplot(data, aes(x = value)) +
  
  # 1. Stacked Density
  geom_density(aes(fill = adjuster_label),
               alpha = 1.0,
               position = "stack", 
               color = "white",
               size = 0.1) +  

  # 2. Manual Boxplot (Rect)
  geom_rect(data = box_stats,
            aes(xmin = q1, xmax = q3, ymin = box_bottom, ymax = box_top),
            fill = NA,            
            color = "black",        
            size = 0.4,
            inherit.aes = FALSE) + 
  
  # 3. Median Line
  geom_segment(data = box_stats,
               aes(x = median, xend = median, y = box_bottom, yend = box_top),
               color = "black",
               size = 1.2,
               inherit.aes = FALSE) +
  
  # 4. Whiskers
  geom_segment(data = box_stats,
               aes(x = lower_whisker, xend = q1, y = mid_y, yend = mid_y),
               color = "black", inherit.aes = FALSE) +
  geom_segment(data = box_stats,
               aes(x = q3, xend = upper_whisker, y = mid_y, yend = mid_y),
               color = "black", inherit.aes = FALSE) +
  
  # 5. Whisker Caps
  geom_segment(data = box_stats,
               aes(x = lower_whisker, xend = lower_whisker, 
                   y = box_bottom + (box_top-box_bottom)*0.25, 
                   yend = box_top - (box_top-box_bottom)*0.25),
               color = "black", inherit.aes = FALSE) +
  geom_segment(data = box_stats,
               aes(x = upper_whisker, xend = upper_whisker, 
                   y = box_bottom + (box_top-box_bottom)*0.25, 
                   yend = box_top - (box_top-box_bottom)*0.25),
               color = "black", inherit.aes = FALSE) +

  # 6. Outliers
  geom_point(data = outliers,
             aes(x = value, y = y_pos),
             shape = 16, 
             size = 0.8, 
             alpha = 0.5,
             inherit.aes = FALSE) +

  facet_wrap(~ classifier_label, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_fill_manual(values = palette_map, name = "Adjuster") +
  
  # Expand bottom to fit boxplot
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) + 
  
  labs(
    x = "MCC Distribution",
    y = NULL
  ) +
  
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
    strip.placement = "outside",
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved plot to:", args$output, "\n")