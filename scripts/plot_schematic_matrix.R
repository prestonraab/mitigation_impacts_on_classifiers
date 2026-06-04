#!/usr/bin/env Rscript

# plot_schematic_matrix.R
# Generate a schematic diagram of a data matrix based on real data.

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ComplexHeatmap)
  library(circlize)
})

parser <- ArgumentParser(description = "Generate schematic data matrix diagram")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input RData file (e.g. data/TB_real_data.RData)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 8)
parser$add_argument("--height", type = "double", default = 6)
parser$add_argument("--dpi", type = "integer", default = 300)

args <- parser$parse_args()

# ==============================================================================
# Data Preparation
# ==============================================================================

# Load real data
load(args$input) # Loads dat_lst, label_lst

# Use the first dataset in the list
# Rows are genes, columns are samples in this dataset typically. 
# We need 6 samples (Patients A-F) and some genes.
dataset <- dat_lst[[1]]
mat <- as.matrix(dataset[1:100, 1:6]) # 100 genes, 6 samples

# Scale rows for better visualization
mat <- t(scale(t(mat)))

# Transpose so samples are rows, genes are columns (as requested)
mat_plot <- t(mat)

# Label Rows: Patient A through Patient F
rownames(mat_plot) <- paste("Patient", LETTERS[1:6])

# Label Columns: "Gene 1" ... "Gene 20,000"
# We only show a subset but label the ends to represent the full range
colnames(mat_plot) <- rep("", ncol(mat_plot))
colnames(mat_plot)[1] <- "Gene 1"
colnames(mat_plot)[ncol(mat_plot)] <- "Gene 20,000"

# ==============================================================================
# Visualization
# ==============================================================================

# Define color gradient: Blue (Low) to Red (High)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

png(args$output, width = args$width, height = args$height, units = "in", res = args$dpi)

# Draw Heatmap
# Two row groups separated by a gap: Study 1 and Study 2
ht <- Heatmap(mat_plot,
              name = "Expression",
              col = col_fun,
              
              # Row settings
              row_split = factor(c(rep("Study 1", 3), rep("Study 2", 3)), levels = c("Study 1", "Study 2")),
              row_gap = unit(8, "mm"),
              row_title_rot = 0,
              row_title_gp = gpar(fontsize = 14, fontface = "bold"),
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 12),
              
              # Column settings
              cluster_columns = FALSE,
              column_title = "High-Dimensional Features",
              column_title_side = "bottom",
              column_names_side = "top",
              column_names_rot = 0,
              column_names_centered = TRUE,
              column_names_gp = gpar(fontsize = 12, fontface = "italic"),
              
              # Aesthetics
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              border = TRUE,
              
              # Heatmap legend
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(
                title = "Relative\nExpression",
                at = c(-2, 0, 2),
                labels = c("Low", "Mid", "High")
              )
)

draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

cat("Saved schematic matrix to:", args$output, "\n")
