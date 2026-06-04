#!/usr/bin/env Rscript

# plot_gapr_heatmap.R
# Generate a hybrid GAPR/ComplexHeatmap for adjuster-classifier performance differences.

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyr)
  library(GAPR)
  library(ComplexHeatmap)
  library(circlize)
  library(dendextend)
})

# Load utils if present
if (file.exists("scripts/adjuster_plot_utils.R")) {
  source("scripts/adjuster_plot_utils.R")
}

# Estimate Hodges-Lehmann pseudomedian
calc_hodges_lehmann <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA)
  return(median(outer(x, x, "+") / 2))
}

parser <- ArgumentParser(description = "Generate hybrid GAPR/ComplexHeatmap")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 8)
parser$add_argument("--height", type = "double", default = 4)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")

args <- parser$parse_args()

# ==============================================================================
# Constants & Helpers
# ==============================================================================

classifier_labels_map <- c(
  "rda" = "RDA", "logistic" = "LR", "elasticnet" = "ENet",
  "svm" = "SVM", "rf" = "RF", "knn" = "KNN",
  "xgboost" = "XGB", "nnet" = "NN", "shrinkageLDA" = "RDA"
)

target_adjusters <- if (!is.null(args$adjusters)) trimws(strsplit(args$adjusters, ",")[[1]]) else NULL

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
    !classifier %in% c("logistic")
  )

# Prepare performance data per study
mcc_data <- raw_data %>%
  mutate(
    mcc = value, 
    .by = c(classifier, n_datasets, test_study)
  ) %>%
  mutate(
    adjuster_label = if (exists("format_adjuster_label")) {
      sapply(adjuster, format_adjuster_label)
    } else {
      adjuster
    }
  )

# Aggregate stats per classifier-adjuster pair
classifier_stats <- mcc_data %>%
  mutate(classifier_label = ifelse(classifier %in% names(classifier_labels_map), 
                                   classifier_labels_map[classifier], 
                                   classifier)) %>%
  summarise(
    avg_mcc = mean(mcc),
    .by = c(classifier_label, adjuster_label)
  )

# Calculate group medians per adjuster
group_medians <- mcc_data %>%
  summarise(
    center = calc_hodges_lehmann(mcc),
    .by = adjuster_label
  )

# Calculate difference: avg_mcc - center
heatmap_data <- classifier_stats %>%
  left_join(group_medians, by = "adjuster_label") %>%
  mutate(mcc_diff = avg_mcc - center) %>%
  mutate(mcc_diff = ifelse(is.na(mcc_diff), 0, mcc_diff)) %>%
  select(adjuster_label, classifier_label, mcc_diff) %>%
  pivot_wider(names_from = classifier_label, values_from = mcc_diff)

# Convert to matrix
mat <- as.matrix(heatmap_data[,-1])
rownames(mat) <- heatmap_data$adjuster_label

# Handle NAs
mat[is.na(mat)] <- 0

# ==============================================================================
# Visualization using Hybrid GAPR + ComplexHeatmap Approach
# ==============================================================================

mat_df <- as.data.frame(mat)

# Run GAP silently to get analytical order and proximity
tmp_png <- tempfile(fileext = ".png")
gap_result <- GAPR::GAP(
  data = mat_df, 
  row.name = rownames(mat),
  row.prox = 'euclidean', 
  col.prox = 'euclidean',
  row.order = 'average', 
  col.order = 'average',
  row.flip = 'r2e', 
  col.flip = 'r2e',
  exp.row_order = TRUE, 
  exp.column_order = TRUE,
  exp.row_names = TRUE, 
  exp.column_names = TRUE,
  exp.row_prox = TRUE, 
  PNGfilename = tmp_png,
  show.plot = FALSE
)
unlink(tmp_png)

# Reconstruct unordered proximity matrix
inv_row_order <- order(gap_result$row_order)
unordered_row_prox <- gap_result$row_prox[inv_row_order, inv_row_order]
rownames(unordered_row_prox) <- rownames(mat)
colnames(unordered_row_prox) <- rownames(mat)

# Recreate tree and rotate to match R2E output
row_dist <- as.dist(unordered_row_prox)
row_hc <- hclust(row_dist, method = "average")
row_dend <- as.dendrogram(row_hc)
row_dend <- dendextend::rotate(row_dend, gap_result$row_names)

# Build color scales
max_val <- max(abs(mat), na.rm = TRUE)
col_fun_main <- circlize::colorRamp2(
  breaks = c(-max_val, 0, max_val),
  colors = c("#D73027", "#ffffbf70", "#4575B4")
)

col_fun_row_prox <- circlize::colorRamp2(
  breaks = seq(min(unordered_row_prox, na.rm=TRUE), max(unordered_row_prox, na.rm=TRUE), length.out=9),
  colors = RColorBrewer::brewer.pal(9, "YlGnBu")
)

# Draw row similarity matrix (1st on left)
ht_row_prox <- Heatmap(
  unordered_row_prox,
  name = "Adjuster Distance",
  col = col_fun_row_prox,
  cluster_rows = row_dend,       
  cluster_columns = row_dend,    
  row_order = gap_result$row_order, 
  column_order = gap_result$row_order, 
  show_row_dend = FALSE,         
  show_column_dend = TRUE,       # Put dendrogram on top
  column_dend_side = "top",
  column_dend_height = unit(2, "cm"),
  show_row_names = TRUE,         # Put row names on far left
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 11.2, fontfamily = "sans"),
  show_column_names = FALSE,
  border = FALSE,
  width = unit(2, "null"),         # 2:1 ratio width
  height = unit(1, "null"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 11.2, fontfamily = "sans"),
    labels_gp = gpar(fontsize = 11.2, fontfamily = "sans")
  )
)

# Draw main matrix (2nd on right)
ht_main <- Heatmap(
  mat,
  name = "MCC Deviation",
  col = col_fun_main,
  cluster_rows = row_dend,       # Maintain alignment with left matrix
  cluster_columns = FALSE,                
  column_order = gap_result$column_order, 
  show_row_names = FALSE,        
  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 11.2, fontfamily = "sans"),
  border = FALSE,
  width = unit(1, "null"),         # Skinny width (half of similarity matrix)
  height = unit(1, "null"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 11.2, fontfamily = "sans"),
    labels_gp = gpar(fontsize = 11.2, fontfamily = "sans")
  )
)

# Render plot to file
png(filename = args$output, width = args$width, height = args$height, units = "in", res = args$dpi, bg = "white")
draw(ht_row_prox + ht_main, ht_gap = unit(5, "mm"), padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

cat("Saved Custom Heatmap PNG to:", args$output, "\n")