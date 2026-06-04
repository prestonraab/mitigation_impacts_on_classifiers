#!/usr/bin/env Rscript

# create_visualization_grid.R
# Create grid plots combining multiple visualizations for comparison

suppressMessages(suppressWarnings({
  library(argparse)
  library(ggplot2)
  library(gridExtra)
  library(png)
  library(grid)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Create grid visualizations comparing adjusters")

parser$add_argument("--input-dir", type = "character", required = TRUE,
                   help = "Input directory containing individual plots")
parser$add_argument("--output-file", type = "character", required = TRUE,
                   help = "Output file path for grid plot")
parser$add_argument("--method", type = "character", required = TRUE,
                   help = "Visualization method: pca, lda, or umap")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets to visualize")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Test study name")
parser$add_argument("--adjusters", type = "character", required = TRUE,
                   help = "Comma-separated list of adjusters to include")

args <- parser$parse_args()

# ====================================================================
# GRID CREATION
# ====================================================================

#' Load PNG as grob
#' @param filepath Path to PNG file
#' @return Grid object
load_png_as_grob <- function(filepath) {
  if (!file.exists(filepath)) {
    warning(sprintf("File not found: %s", filepath))
    return(textGrob(sprintf("Missing:\n%s", basename(filepath))))
  }
  
  img <- readPNG(filepath)
  rasterGrob(img, interpolate = TRUE)
}

#' Create comparison grid
#' @param input_dir Input directory
#' @param method Visualization method
#' @param num_datasets Number of datasets
#' @param test_study Test study name
#' @param adjusters Vector of adjuster names
#' @return Grid arrangement
create_comparison_grid <- function(input_dir, method, num_datasets, test_study, adjusters) {
  
  # Load all plots
  grobs <- lapply(adjusters, function(adj) {
    filename <- sprintf("%s_n%d_test%s.png", adj, num_datasets, test_study)
    filepath <- file.path(input_dir, method, filename)
    load_png_as_grob(filepath)
  })
  
  # Create grid with labels
  title <- sprintf("%s Comparison (n=%d datasets, test=%s)", 
                  toupper(method), num_datasets, test_study)
  
  # Arrange horizontally (all adjusters in one row)
  do.call(grid.arrange, c(
    grobs,
    list(
      ncol = length(grobs),
      top = textGrob(title, gp = gpar(fontsize = 16, fontface = "bold"))
    )
  ))
}

# ====================================================================
# MAIN EXECUTION
# ====================================================================

# Parse adjusters list
adjusters <- strsplit(args$adjusters, ",")[[1]]
adjusters <- trimws(adjusters)  # Remove any whitespace

cat(sprintf("Creating %s grid for n=%d, test=%s\n", 
            args$method, args$num_datasets, args$test_study))
cat(sprintf("Adjusters: %s\n", paste(adjusters, collapse = ", ")))

# Create output directory
dir.create(dirname(args$output_file), recursive = TRUE, showWarnings = FALSE)

# Create and save grid
png(args$output_file, width = 800 * length(adjusters), height = 800, res = 150, bg = "white")
create_comparison_grid(args$input_dir, args$method, args$num_datasets, args$test_study, adjusters)
dev.off()

cat(sprintf("Saved grid: %s\n", args$output_file))
