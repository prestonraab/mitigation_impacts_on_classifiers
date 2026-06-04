#!/usr/bin/env Rscript

# plot_pca_demonstration.R
# Generate side-by-side PCA plots showing data before and after batch correction.

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(sva)
  library(cowplot)
})

parser <- ArgumentParser(description = "Generate side-by-side PCA demonstration")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input RData file (e.g. data/TB_real_data.RData)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 12)
parser$add_argument("--height", type = "double", default = 6)
parser$add_argument("--dpi", type = "integer", default = 300)

args <- parser$parse_args()

# ==============================================================================
# Data Preparation
# ==============================================================================

load(args$input) # Loads dat_lst, label_lst

# Use two different studies
study1_expr <- dat_lst[[1]]
study1_y <- label_lst[[1]]

study2_expr <- dat_lst[[2]]
study2_y <- label_lst[[2]]

# Intersect genes
common_genes <- intersect(rownames(study1_expr), rownames(study2_expr))
study1_expr <- study1_expr[common_genes, ]
study2_expr <- study2_expr[common_genes, ]

# Combine data
raw_expr <- cbind(study1_expr, study2_expr)
metadata <- data.frame(
  sample_id = colnames(raw_expr),
  study = c(rep("Study 1", ncol(study1_expr)), rep("Study 2", ncol(study2_expr))),
  disease = factor(c(study1_y, study2_y), labels = c("Latent TB", "Active TB"))
)

# Batch Correction (ComBat)
batch <- metadata$study
mod <- model.matrix(~disease, data = metadata)
corrected_expr <- ComBat(dat = as.matrix(raw_expr), batch = batch, mod = mod, par.prior = TRUE)

# ==============================================================================
# PCA Calculations
# ==============================================================================

# Raw data PCA
pca_raw <- prcomp(t(raw_expr), scale. = TRUE)
df_raw <- as.data.frame(pca_raw$x[, 1:2])
df_raw <- cbind(df_raw, metadata)
var_raw <- round(100 * pca_raw$sdev^2 / sum(pca_raw$sdev^2), 1)

# Corrected data PCA
pca_corr <- prcomp(t(corrected_expr), scale. = TRUE)
df_corr <- as.data.frame(pca_corr$x[, 1:2])
df_corr <- cbind(df_corr, metadata)
var_corr <- round(100 * pca_corr$sdev^2 / sum(pca_corr$sdev^2), 1)

# Shared Styling Options
study_colors <- c("Study 1" = "#E03A3A", "Study 2" = "#2F69E0") # Bold Red and Bold Blue

# Custom 4-color palette for Right Plot (Study x Disease Saturation)
df_corr$group <- factor(
  paste(df_corr$study, df_corr$disease, sep = " - "),
  levels = c(
    "Study 1 - Latent TB",
    "Study 1 - Active TB",
    "Study 2 - Latent TB",
    "Study 2 - Active TB"
  )
)

group_colors <- c(
  "Study 1 - Latent TB" = "#C49090",  # Red-grey (desaturated)
  "Study 1 - Active TB" = "#E03A3A",  # Bold Red (saturated)
  "Study 2 - Latent TB" = "#94ACD4",  # Blue-grey (desaturated)
  "Study 2 - Active TB" = "#2F69E0"   # Bold Blue (saturated)
)

group_shapes <- c(
  "Study 1 - Latent TB" = 17,  # Triangle
  "Study 1 - Active TB" = 16,  # Circle
  "Study 2 - Latent TB" = 17,  # Triangle
  "Study 2 - Active TB" = 16   # Circle
)

# ==============================================================================
# Visualization
# ==============================================================================

# Left Plot: Raw Data, colored by Study (Red/Blue), shaped by Disease
p1 <- ggplot(df_raw, aes(x = PC1, y = PC2, color = study, shape = disease)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = study_colors) +
  scale_shape_manual(values = c("Latent TB" = 17, "Active TB" = 16)) +
  labs(
    title = "Before Batch Correction",
    x = paste0("PC1 (", var_raw[1], "%)"),
    y = paste0("PC2 (", var_raw[2], "%)"),
    color = "Dataset (Study)",
    shape = "Disease State"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    panel.border = element_rect(color = "gray90", fill = NA, size = 0.8),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 5, r = 0, b = 0, l = 0),
    legend.text = element_text(margin = margin(r = 25, l = 3)),
    legend.spacing.x = unit(0.2, "cm"),
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 10))
  ) +
  guides(
    color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, order = 1),
    shape = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, order = 2)
  )

# Right Plot: Corrected Data, colored & shaped by Study/Disease interaction (auto-merged legend)
p2 <- ggplot(df_corr, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(
    title = "After Batch Correction (ComBat)",
    x = paste0("PC1 (", var_corr[1], "%)"),
    y = paste0("PC2 (", var_corr[2], "%)"),
    color = "Dataset & Disease State",
    shape = "Dataset & Disease State"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95"),
    panel.border = element_rect(color = "gray90", fill = NA, size = 0.8),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 5, r = 0, b = 0, l = 0),
    legend.text = element_text(margin = margin(r = 25, l = 3)),
    legend.spacing.x = unit(0.2, "cm"),
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 10))
  ) +
  guides(
    color = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 2, byrow = TRUE),
    shape = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 2, byrow = TRUE)
  )

# Combine with cowplot
p_combined <- plot_grid(p1, p2, ncol = 2, labels = "AUTO")

# Save output
save_plot(args$output, p_combined, base_width = args$width, base_height = args$height, dpi = args$dpi)
cat("Saved PCA demonstration to:", args$output, "\n")
