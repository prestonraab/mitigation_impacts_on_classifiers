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
# Define Color Palette (Logic from Snippet 1)
# ==============================================================================

group_highlight_1 <- c("Within-Study CV")
group_highlight_2 <- c("ComBat Sup.")
group_highlight_3 <- c("Coconut")
group_aggregate_1 <- c("ReComBat", "CuBlock", "ComBat", "ComBat Mean", "NPN", "Rank Twice", "Naive")
group_aggregate_2 <- c()

# Individual colors for MNN family + Rank Features
group_mnn      <- c("MNN")
group_fastmnn  <- c("FastMNN")
group_rank_feat <- c("Rank Features", "YuGene", "Angel")
group_log_only <- c("Log Only", "Shambhala2, RUVg", "RNAbc")

all_adjusters <- levels(data$adjuster_label)
palette_map <- character(length(all_adjusters))
names(palette_map) <- all_adjusters

# ==============================================================================
# Define Color Palette (Modern Scientific Aesthetic)
# ==============================================================================

# Custom professional palette
color_vermilion <- "#D55E00" # Strong Highlight
color_blue      <- "#7c3bb4ff" # Secondary Highlight
color_purple    <- "#bfa6d5ff" # Tertiary Highlight
color_slate     <- "#acb1b7ff" # Muted Gray 
color_mnn       <- "#4292C6" # Deeper Blue (MNN)
color_fastmnn   <- "#6ab7e0ff" # Light-Medium Blue (FastMNN)
color_rank_feat <- "#BDD7E7" # Light Blue (Rank Features)
color_log_only  <- "#636363" # Dark Grey (Log Only)
color_fallback   <- "#F0E442" # Warning

all_adjusters <- levels(data$adjuster_label)
palette_map <- character(length(all_adjusters))
names(palette_map) <- all_adjusters

for (lbl in all_adjusters) {
  if (lbl %in% group_highlight_1) {
    palette_map[lbl] <- color_vermilion
  } else if (lbl %in% group_highlight_2) {
    palette_map[lbl] <- color_blue
  } else if (lbl %in% group_aggregate_1) {
    palette_map[lbl] <- color_slate
  } else if (lbl %in% group_aggregate_2) {
    palette_map[lbl] <- color_sky
  } else if (lbl %in% group_mnn) {
    palette_map[lbl] <- color_mnn
  } else if (lbl %in% group_fastmnn) {
    palette_map[lbl] <- color_fastmnn
  } else if (lbl %in% group_rank_feat) {
    palette_map[lbl] <- color_rank_feat
  } else if (lbl %in% group_log_only) {
    palette_map[lbl] <- color_log_only
  } else {
    palette_map[lbl] <- color_fallback
  }
}

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