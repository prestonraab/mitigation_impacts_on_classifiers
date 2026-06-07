#!/usr/bin/env Rscript

# plot_class_imbalance_trend.R - Plot adjuster performance trends across class imbalance levels
# Aggregation strategy:
# 1. Average MCC over training pairs and replicates -> (classifier, adjuster, test_set, imbalance)
# 2. Rank adjusters within (classifier, test_set, imbalance)
# 3. Average ranks over classifiers -> (adjuster, test_set, imbalance)
# 4. Plot mean rank over test sets with min/max uncertainty zone

library(argparse)
library(dplyr)
library(ggplot2)
library(tidyr)

# Parse command line arguments
parser <- ArgumentParser(description = "Create class imbalance trend plot with aggregated rankings")
parser$add_argument("--input-data", required = TRUE, help = "Path to class imbalanced results CSV")
parser$add_argument("-o", "--output", required = TRUE, help = "Output PNG file path")

args <- parser$parse_args()

# Load data
cat("Loading class imbalanced results from:", args$input_data, "\n")
data <- read.csv(args$input_data, stringsAsFactors = FALSE)

cat("Data dimensions:", nrow(data), "rows\n")
cat("Adjusters:", paste(unique(data$adjuster), collapse = ", "), "\n")
cat("Classifiers:", paste(unique(data$classifier), collapse = ", "), "\n")
cat("Imbalance levels:", paste(sort(unique(data$imbalance_pct)), collapse = ", "), "\n")

cat("Columns in data:", paste(names(data), collapse = ", "), "\n")

# Check training pairs
cat("\nTraining pairs in data:\n")
training_pairs_summary <- data %>%
  select(training_pair, train_dataset_1, train_dataset_2, test_dataset) %>%
  distinct() %>%
  arrange(training_pair, test_dataset)
print(head(training_pairs_summary, 20))
cat("Number of unique training pairs:", n_distinct(data$training_pair), "\n")
cat("Number of unique test datasets:", n_distinct(data$test_dataset), "\n")

# Filter to target adjusters (all adjusters present in data)
target_adjusters <- unique(data$adjuster)
data_filtered <- data %>%
  filter(adjuster %in% target_adjusters, classifier != "knn")

cat("Filtered to", nrow(data_filtered), "rows with target adjusters\n")

# Fix training_pair column - create proper identifier from training datasets
cat("\nFixing training_pair column...\n")
data_filtered <- data_filtered %>%
  mutate(
    # Sort the two training datasets alphabetically and combine
    training_pair_fixed = paste(
      pmin(train_dataset_1, train_dataset_2),
      pmax(train_dataset_1, train_dataset_2),
      sep = "-"
    )
  )

cat("Original training_pair values:", n_distinct(data_filtered$training_pair), "unique\n")
cat("Fixed training_pair values:", n_distinct(data_filtered$training_pair_fixed), "unique\n")
cat("Sample of fixed training pairs:\n")
print(head(unique(data_filtered$training_pair_fixed), 10))

# Replace the old column
data_filtered$training_pair <- data_filtered$training_pair_fixed
data_filtered$training_pair_fixed <- NULL

# Step 1: Average MCC over replicates
# Result: one value per (classifier, adjuster, test_set, training_pair, imbalance)
cat("\nStep 1: Averaging MCC over replicates...\n")
mcc_averaged <- data_filtered %>%
  group_by(classifier, adjuster, test_dataset, training_pair, imbalance_pct) %>%
  summarise(
    mean_mcc = mean(mcc, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

cat("  Result:", nrow(mcc_averaged), "unique combinations\n")
cat("    Classifiers:", n_distinct(mcc_averaged$classifier), "\n")
cat("    Adjusters:", n_distinct(mcc_averaged$adjuster), "\n")
cat("    Test sets:", n_distinct(mcc_averaged$test_dataset), "\n")
cat("    Training pairs:", n_distinct(mcc_averaged$training_pair), "\n")
cat("    Imbalance levels:", n_distinct(mcc_averaged$imbalance_pct), "\n")



# Step 2: Rank adjusters within each (classifier, test_set, training_pair) group
cat("\nStep 2: Ranking adjusters within (classifier, test_set, training_pair)...\n")
ranked_by_classifier <- mcc_averaged %>%
  group_by(classifier, test_dataset, training_pair) %>%
  arrange(desc(mean_mcc)) %>%
  mutate(rank = rank(-mean_mcc, ties.method = "average")) %>%
  ungroup()

cat("  Ranked", nrow(ranked_by_classifier), "combinations\n")

# Step 3: Average ranks over classifiers and training pairs
# Result: one value per (adjuster, test_set, imbalance)
cat("\nStep 3: Averaging ranks over classifiers and training pairs...\n")
ranks_averaged_over_classifiers <- ranked_by_classifier %>%
  group_by(adjuster, test_dataset, imbalance_pct) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

cat("  Result:", nrow(ranks_averaged_over_classifiers), "unique (adjuster, test_set, imbalance) combinations\n")
cat("  Adjusters:", paste(unique(ranks_averaged_over_classifiers$adjuster), collapse = ", "), "\n")
cat("  Test sets:", n_distinct(ranks_averaged_over_classifiers$test_dataset), "\n")
cat("  Imbalance levels:", paste(sort(unique(ranks_averaged_over_classifiers$imbalance_pct)), collapse = ", "), "\n")


# Step 4: Calculate mean, min, max over test sets for plotting
cat("\nStep 4: Computing statistics over test sets...\n")

plot_data <- ranks_averaged_over_classifiers %>%
  group_by(adjuster, imbalance_pct) %>%
  summarise(
    avg_rank = mean(mean_rank, na.rm = TRUE),
    min_rank = min(mean_rank, na.rm = TRUE),
    max_rank = max(mean_rank, na.rm = TRUE),
    n_test_sets = n(),
    .groups = "drop"
  ) %>%
  rename(mean_rank = avg_rank)  # Rename back for compatibility

cat("  Final plot data:", nrow(plot_data), "points\n")

# Create adjuster labels for plotting
adjuster_label_map <- c(
  "unadjusted"   = "Unadjusted",
  "naive"        = "Naive",
  "rank_twice"   = "RankTwice",
  "combat"       = "ComBat (Unsupervised)",
  "combat_sup"   = "ComBat (Supervised)",
  "coconut"      = "COCONUT (Supervised)",
  "recombat"     = "reComBat (Unsupervised)",
  "rankin"       = "Rank-In"
)

adjuster_level_order <- c(
  "Unadjusted", "Naive", "RankTwice",
  "ComBat (Unsupervised)", "ComBat (Supervised)",
  "COCONUT (Supervised)",
  "reComBat (Unsupervised)",
  "Rank-In"
)

adjuster_colors <- c(
  "Unadjusted"             = "#E31A1C",
  "Naive"                  = "#999999",
  "RankTwice"              = "#6A3D9A",
  "ComBat (Unsupervised)"  = "#1F78B4",
  "ComBat (Supervised)"    = "#33A02C",
  "COCONUT (Supervised)"   = "#FF7F00",
  "reComBat (Unsupervised)"= "#A6CEE3",
  "Rank-In"                = "#B15928"
)

plot_data <- plot_data %>%
  mutate(
    adjuster_label = adjuster_label_map[adjuster],
    adjuster_label = factor(adjuster_label,
                            levels = adjuster_level_order[adjuster_level_order %in% adjuster_label]),
    imbalance_pct_num = imbalance_pct * 100
  )

# Print summary
cat("\nSummary of plot data:\n")
print(plot_data %>% select(adjuster_label, imbalance_pct_num, mean_rank, min_rank, max_rank))

# Create the plot
p <- ggplot(plot_data, aes(x = imbalance_pct_num, y = mean_rank, color = adjuster_label, fill = adjuster_label)) +
  geom_ribbon(aes(ymin = min_rank, ymax = max_rank), alpha = 0.2, color = NA) +
  # Mean line
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  # Reverse y-axis so rank 1 is at top
  scale_y_reverse(
    breaks = seq(1, 45, by = 5)
  ) +
  scale_x_reverse(
    breaks = sort(unique(plot_data$imbalance_pct_num), decreasing = TRUE),
    labels = paste0(sort(unique(plot_data$imbalance_pct_num), decreasing = TRUE), "%")
  ) +
  scale_color_manual(values = adjuster_colors, name = "Batch Adjuster") +
  scale_fill_manual(values = adjuster_colors, name = "Batch Adjuster") +
  labs(
    x = "Class Imbalance Level (% Active TB in High-Imbalance Training Set)",
    y = "Average Performance Rank"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Save the plot
ggsave(args$output, p, width = 10, height = 7, dpi = 300)

cat("\nPlot saved to:", args$output, "\n")

# Calculate and print trend statistics
cat("\nTrend analysis (linear regression of rank vs imbalance):\n")
trend_stats <- plot_data %>%
  group_by(adjuster_label) %>%
  do({
    model <- lm(mean_rank ~ imbalance_pct_num, data = .)
    data.frame(
      slope = coef(model)[2],
      intercept = coef(model)[1],
      p_value = summary(model)$coefficients[2, 4],
      r_squared = summary(model)$r.squared
    )
  }) %>%
  ungroup()

print(trend_stats)
