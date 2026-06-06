#!/usr/bin/env Rscript

# visualize_batch_adjustment.R
# Generate PCA, LDA, and UMAP visualizations for batch adjustment methods
# Follows single-responsibility principle with modular functions

# Suppress warnings and messages
options(warn = -1)
suppressMessages(suppressWarnings({
  required_packages <- c("argparse", "ggplot2", "dplyr", "purrr", "umap", "MASS", 
                        "genefilter", "sva", "batchelor", "SummarizedExperiment")
  sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Visualize batch adjustment effects using PCA, LDA, and UMAP")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, naive, rank_samples, rank_twice, npn, combat, combat_mean, combat_sup, mnn, or fast_mnn")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets to include: 3, 4, 5, or 6")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Study to use as test set (e.g., 'India', 'USA', 'Africa')")
parser$add_argument("--output-dir", type = "character", required = TRUE,
                   help = "Output directory for visualization files")
parser$add_argument("--reduce", type = "integer", default = 0,
                   help = "Number of dimensions to reduce to (default: 0, no reduction)")

args <- parser$parse_args()

# Validate arguments
valid_adjusters <- c("unadjusted", "naive", "rank_samples", "rank_twice", "npn", "combat", "combat_mean", "combat_sup", "mnn", "fast_mnn", "ruvr", "gmm", "pace_default", "pace_aggressive", "pace_focused", "pace_conservative", "pace_ultra_aggressive", "pace_extreme_aggressive", "pace_iterative_aggressive")
valid_num_datasets <- c(2, 3, 4, 5)

if (!args$adjuster %in% valid_adjusters) {
  stop(sprintf("Invalid adjuster '%s'. Must be one of: %s", 
               args$adjuster, paste(valid_adjusters, collapse = ", ")))
}

if (!args$num_datasets %in% valid_num_datasets) {
  stop(sprintf("Invalid num-datasets '%d'. Must be one of: %s", 
               args$num_datasets, paste(valid_num_datasets, collapse = ", ")))
}

# Extract parameters
adjuster <- args$adjuster
num_datasets <- args$num_datasets
test_study <- args$test_study
output_dir <- args$output_dir
reduce <- args$reduce

# Create output directory structure
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "pca"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "lda"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "umap"), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Starting visualization: adjuster=%s, num_datasets=%d, test_study=%s\n", 
            adjuster, num_datasets, test_study))

# Load helper functions
source("scripts/helper.R")


# ====================================================================
# DATA LOADING AND PREPARATION (Single Responsibility: Data I/O)
# ====================================================================

#' Load and filter TB real data
#' @param n_studies Number of studies to include
#' @return List with filtered data and labels
load_and_filter_data <- function(n_studies) {
  data_path <- "data/TB_real_data.RData"
  if (!file.exists(data_path)) {
    stop(sprintf("Data file not found: %s", data_path))
  }
  
  load(data_path)
  
  all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
  selected_studies <- all_studies[1:n_studies]
  
  dat_lst <- dat_lst[selected_studies]
  label_lst <- label_lst[selected_studies]
  
  cat(sprintf("Loaded %d studies: %s\n", 
              n_studies, paste(selected_studies, collapse = ", ")))
  
  list(
    dat_lst = dat_lst,
    label_lst = label_lst,
    study_names = selected_studies
  )
}

#' Validate test study selection
#' @param study_names Vector of study names
#' @param test_study Name of test study
#' @return Name of test study (validated)
validate_test_study <- function(study_names, test_study) {
  if (!test_study %in% study_names) {
    stop(sprintf("Test study '%s' not found in available studies: %s",
                 test_study, paste(study_names, collapse = ", ")))
  }
  
  cat(sprintf("Using test study: %s\n", test_study))
  cat(sprintf("Training studies: %s\n", 
              paste(setdiff(study_names, test_study), collapse = ", ")))
  
  test_study
}

#' Prepare training and test datasets
#' @param dat_lst List of data matrices
#' @param label_lst List of label vectors
#' @param test_name Name of test study
#' @param study_names All study names
#' @return List with prepared datasets and metadata
prepare_train_test_split <- function(dat_lst, label_lst, test_name, study_names) {
  train_names <- setdiff(study_names, test_name)
  
  # Combine training data
  dat <- do.call(cbind, dat_lst[train_names])
  
  # Create batch assignments using study names (not numeric IDs)
  batch <- rep(train_names, times = sapply(dat_lst[train_names], ncol))
  batches_ind <- lapply(train_names, function(name) which(batch == name))
  names(batches_ind) <- train_names
  
  group <- do.call(c, label_lst[train_names])
  
  # Test data
  dat_test <- dat_lst[[test_name]]
  group_test <- label_lst[[test_name]]
  
  cat(sprintf("Training: %d samples from %d batches\n", ncol(dat), length(train_names)))
  cat(sprintf("  Batches: %s\n", paste(train_names, collapse = ", ")))
  cat(sprintf("Test: %d samples from %s\n", ncol(dat_test), test_name))
  
  list(
    dat = dat,
    batch = batch,
    batches_ind = batches_ind,
    batch_names = train_names,
    group = group,
    dat_test = dat_test,
    group_test = group_test,
    test_name = test_name
  )
}

#' Reduce features to top N most variable genes
#' @param dat Training data matrix
#' @param dat_test Test data matrix
#' @param n_genes Number of genes to select
#' @return List with reduced datasets
reduce_features <- function(dat, dat_test, n_genes = 1000) {
  genes_sel_idx <- order(rowVars(dat), decreasing = TRUE)[1:n_genes]
  
  cat(sprintf("Selected top %d most variable genes\n", n_genes))
  
  list(
    dat = dat[genes_sel_idx, ],
    dat_test = dat_test[genes_sel_idx, ]
  )
}


# ====================================================================
# BATCH CORRECTION (Single Responsibility: Data Transformation)
# ====================================================================

# ====================================================================
# RANK ADJUSTMENT HELPER FUNCTIONS
# ====================================================================

rank_normalized <- function(matrix_, dim) {
  if (dim < 1 || dim > 2) {
    stop("Invalid dimension. Must be 1 for rows or 2 for columns.")
  }
  ranked = apply(matrix_, dim, rank, ties.method = "average")
  
  # apply() transposes the result when dim=1, so we need to transpose it back
  # When dim=1: apply ranks across columns (samples) for each row (feature)
  # When dim=2: apply ranks across rows (features) for each column (sample)
  if (dim == 1 && is.matrix(ranked)) {
    ranked = t(ranked)
  }
  
  return(ranked / max(ranked, na.rm = TRUE))
}

adjust_ranked_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  #' Normalize sample-wise by ranking the genes within the sample, and then by batch.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @return The adjusted matrix (features x samples).
  
  cat("Adjusting with ranked with batch info.\n")
  ranked = rank_normalized(matrix_, 1)
  
  if (debug) {
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
    cat("DEBUG: ranked dimensions: ", nrow(ranked), " x ", ncol(ranked), "\n")
  }
  
  batch_levels <- unique(batch)
  ranked2 <- matrix(NA, nrow = nrow(ranked), ncol = ncol(ranked))
  
  for (b in batch_levels) {
    # For each batch, we rank by sample.
    batch_indices <- which(batch == b)
    batch_data <- ranked[, batch_indices, drop = FALSE]
    
    if (debug) {
      cat("DEBUG: Processing batch '", b, "' with ", length(batch_indices), " samples\n")
      cat("DEBUG: batch_data dimensions: ", nrow(batch_data), " x ", ncol(batch_data), "\n")
    }
    
    # Only apply ranking if there's more than one sample in the batch
    if (ncol(batch_data) > 1) {
      batch_ranked <- rank_normalized(batch_data, 2)
      if (debug) {
        cat("DEBUG: batch_ranked dimensions: ", nrow(batch_ranked), " x ", ncol(batch_ranked), "\n")
      }
      ranked2[, batch_indices] <- batch_ranked
    } else {
      # For single-sample batches, just use the original ranked values
      ranked2[, batch_indices] <- batch_data
    }
  }
  
  # Handle any remaining NA values
  if (any(is.na(ranked2))) {
    cat("WARNING: Found NA values in ranked2 matrix. Replacing with original ranked values.\n")
    ranked2[is.na(ranked2)] <- ranked[is.na(ranked2)]
  }
  
  max_val <- max(ranked2, na.rm = TRUE)
  if (max_val == 0) {
    cat("WARNING: Maximum value in ranked2 is 0. Using 1 as denominator.\n")
    max_val <- 1
  }
  
  return(ranked2 / max_val)
}

adjust_npn <- function(matrix_, batch, debug = FALSE) {
  #' Adjust matrix using Nonparanormal (NPN) transformation.
  #' If batch is NULL, the entire matrix is adjusted at once.
  #' Assumes the batch vector contains no NA values.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @return The adjusted matrix.
  
  if (is.null(batch)) {
    # If no batch is provided, adjust the whole matrix.
    cat("Batch is NULL. Adjusting entire matrix with NPN transformation.\n")
    
    # Transpose to (samples x features) for huge.npn.
    matrix_t <- t(matrix_)
    
    npn_transformed_t <- huge::huge.npn(matrix_t, verbose = FALSE)
    
    return(t(npn_transformed_t))
    
  } else {
    cat("Adjusting using Nonparanormal (NPN) transformation by batch.\n")
    
    # Split the matrix by batch.
    batch_levels <- unique(batch)
    matrix_by_batch <- list()
    
    for (b in batch_levels) {
      batch_indices <- which(batch == b)
      if (length(batch_indices) > 0) {
        matrix_by_batch[[as.character(b)]] <- matrix_[, batch_indices, drop = FALSE]
      }
    }
    
    # Apply NPN transformation to each batch.
    for (b in names(matrix_by_batch)) {
      matrix_t <- t(matrix_by_batch[[b]])
      npn_transformed_t <- huge::huge.npn(matrix_t, verbose = FALSE)
      matrix_by_batch[[b]] <- t(npn_transformed_t)
    }
    
    # Reassemble the matrix from the adjusted batches.
    result_matrix <- matrix_
    for (b in names(matrix_by_batch)) {
      batch_indices <- which(batch == as.character(b))
      result_matrix[, batch_indices] <- matrix_by_batch[[b]]
    }
    
    return(result_matrix)
  }
}

adjust_ranked_samples_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  #' Rank samples within each gene (across samples), then merge batches with normalized ranks.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @return The adjusted matrix (features x samples).
  
  cat("Adjusting with rank_samples (rank samples within genes, batch-aware).\n")
  
  batch_levels <- unique(batch)
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  
  for (b in batch_levels) {
    # For each batch, rank samples within each gene
    batch_indices <- which(batch == b)
    batch_data <- matrix_[, batch_indices, drop = FALSE]
    
    if (debug) {
      cat("DEBUG: Processing batch '", b, "' with ", length(batch_indices), " samples\n")
    }
    
    # Only apply ranking if there's more than one sample in the batch
    if (ncol(batch_data) > 1) {
      # Rank samples within each gene (dim=2 means rank across columns/samples for each row/gene)
      batch_ranked <- rank_normalized(batch_data, 2)
      result_matrix[, batch_indices] <- batch_ranked
    } else {
      # For single-sample batches, just normalize to [0,1]
      result_matrix[, batch_indices] <- batch_data / max(batch_data, na.rm = TRUE)
    }
  }
  
  # Handle any remaining NA values
  if (any(is.na(result_matrix))) {
    cat("WARNING: Found NA values in result_matrix. Replacing with original values.\n")
    result_matrix[is.na(result_matrix)] <- matrix_[is.na(result_matrix)]
  }
  
  max_val <- max(result_matrix, na.rm = TRUE)
  if (max_val == 0) {
    max_val <- 1
  }
  
  return(result_matrix / max_val)
}

adjust_ranked_twice_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  #' Double ranking: first rank genes within samples, then rank samples within genes, batch-aware.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @return The adjusted matrix (features x samples).
  
  cat("Adjusting with rank_twice (genes within samples, then samples within genes, batch-aware).\n")
  
  batch_levels <- unique(batch)
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  
  for (b in batch_levels) {
    # For each batch, apply double ranking
    batch_indices <- which(batch == b)
    batch_data <- matrix_[, batch_indices, drop = FALSE]
    
    if (debug) {
      cat("DEBUG: Processing batch '", b, "' with ", length(batch_indices), " samples\n")
    }
    
    # Only apply ranking if there's more than one sample in the batch
    if (ncol(batch_data) > 1) {
      # First: rank genes within each sample (dim=1)
      # Then: rank samples within each gene (dim=2)
      batch_ranked <- rank_normalized(rank_normalized(batch_data, 1), 2)
      result_matrix[, batch_indices] <- batch_ranked
    } else {
      # For single-sample batches, just apply single ranking (genes within sample)
      batch_ranked <- rank_normalized(batch_data, 1)
      result_matrix[, batch_indices] <- batch_ranked
    }
  }
  
  # Handle any remaining NA values
  if (any(is.na(result_matrix))) {
    cat("WARNING: Found NA values in result_matrix. Replacing with original values.\n")
    result_matrix[is.na(result_matrix)] <- matrix_[is.na(result_matrix)]
  }
  
  max_val <- max(result_matrix, na.rm = TRUE)
  if (max_val == 0) {
    max_val <- 1
  }
  
  return(result_matrix / max_val)
}

#' Apply batch correction method
#' @param dat Training data matrix
#' @param batch Batch assignments (character vector with study names)
#' @param group Sample labels
#' @param dat_test Test data matrix
#' @param method Correction method: "unadjusted", "combat", or "mnn"
#' @return List with corrected training and test data
apply_batch_correction <- function(dat, batch, group, dat_test, method) {
  cat(sprintf("Applying batch correction: %s\n", method))
  
  if (method == "unadjusted") {
    return(list(
      dat_corrected = dat,
      dat_test_corrected = dat_test
    ))
  }
  
  if (method == "naive") {
    # Naive correction: match means and variances across batches
    cat(sprintf("Applying naive batch correction (mean/variance matching)\n"))
    
    # Step 1: Correct training data by standardizing each batch
    dat_corrected <- dat
    unique_batches <- unique(batch)
    
    # Calculate overall statistics across all training data
    overall_mean <- rowMeans(dat)
    overall_var <- apply(dat, 1, var)
    
    # For each batch, standardize to match overall statistics
    for (b in unique_batches) {
      batch_idx <- which(batch == b)
      batch_data <- dat[, batch_idx, drop = FALSE]
      
      # Calculate batch-specific statistics
      batch_mean <- rowMeans(batch_data)
      batch_var <- apply(batch_data, 1, var)
      
      # Avoid division by zero for genes with no variance
      batch_sd <- sqrt(pmax(batch_var, 1e-10))
      overall_sd <- sqrt(pmax(overall_var, 1e-10))
      
      # Standardize: (x - batch_mean) / batch_sd * overall_sd + overall_mean
      for (i in 1:nrow(dat)) {
        if (batch_sd[i] > 1e-10) {
          dat_corrected[i, batch_idx] <- (batch_data[i, ] - batch_mean[i]) / batch_sd[i] * overall_sd[i] + overall_mean[i]
        } else {
          # If no variance in batch, just shift to overall mean
          dat_corrected[i, batch_idx] <- overall_mean[i]
        }
      }
    }
    
    # Step 2: Correct test data to match training data distribution
    test_mean <- rowMeans(dat_test)
    test_var <- apply(dat_test, 1, var)
    test_sd <- sqrt(pmax(test_var, 1e-10))
    
    dat_test_corrected <- dat_test
    for (i in 1:nrow(dat_test)) {
      if (test_sd[i] > 1e-10) {
        dat_test_corrected[i, ] <- (dat_test[i, ] - test_mean[i]) / test_sd[i] * overall_sd[i] + overall_mean[i]
      } else {
        # If no variance in test data, just shift to overall training mean
        dat_test_corrected[i, ] <- overall_mean[i]
      }
    }
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "rank_samples") {
    # Rank samples within genes, batch-aware
    cat(sprintf("Applying rank_samples batch correction\n"))
    
    # Step 1: Apply rank_samples adjustment to training data
    dat_corrected <- adjust_ranked_samples_with_batch_info(dat, batch, debug = FALSE)
    
    # Step 2: Apply rank_samples adjustment to test data
    # Treat test data as a single batch for ranking
    test_batch <- rep(1, ncol(dat_test))
    dat_test_corrected <- adjust_ranked_samples_with_batch_info(dat_test, test_batch, debug = FALSE)
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "rank_twice") {
    # Double ranking: genes within samples, then samples within genes, batch-aware
    cat(sprintf("Applying rank_twice batch correction\n"))
    
    # Step 1: Apply rank_twice adjustment to training data
    dat_corrected <- adjust_ranked_twice_with_batch_info(dat, batch, debug = FALSE)
    
    # Step 2: Apply rank_twice adjustment to test data
    # Treat test data as a single batch for ranking
    test_batch <- rep(1, ncol(dat_test))
    dat_test_corrected <- adjust_ranked_twice_with_batch_info(dat_test, test_batch, debug = FALSE)
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "npn") {
    # NPN (Nonparanormal) quantile normalization
    cat(sprintf("Applying Nonparanormal quantile normalization\n"))
    
    # Load required library
    if (!requireNamespace("huge", quietly = TRUE)) {
      stop("Package 'huge' is required for NPN adjustment but is not installed.")
    }
    
    # Step 1: Apply NPN adjustment to training data
    dat_corrected <- adjust_npn(dat, batch, debug = FALSE)
    
    # Step 2: Apply NPN adjustment to test data
    # Treat test data as a single batch for NPN transformation
    test_batch <- rep(1, ncol(dat_test))
    dat_test_corrected <- adjust_npn(dat_test, test_batch, debug = FALSE)
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "combat") {
    library(sva, quietly = TRUE)
    
    # ComBat correction without labels (unsupervised)
    # Step 1: Correct batch effects within training data without using labels
    dat_corrected <- ComBat(dat, batch=batch, mod=NULL)
    
    # Step 2: Adjust test data to match corrected training distribution
    combined_dat <- cbind(dat_corrected, dat_test)
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                       rep(test_batch_id, ncol(dat_test)))
    
    combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                             mod=NULL, ref.batch=ref_batch_id)
    
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "combat_mean") {
    library(sva, quietly = TRUE)
    
    # ComBat correction with mean adjustment only (no variance adjustment)
    dat_corrected <- ComBat(dat, batch=batch, mod=NULL, mean.only=TRUE)
    
    # Step 2: Adjust test data to match corrected training distribution
    combined_dat <- cbind(dat_corrected, dat_test)
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                       rep(test_batch_id, ncol(dat_test)))
    
    combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                             mod=NULL, ref.batch=ref_batch_id, mean.only=TRUE)
    
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "combat_sup") {
    library(sva, quietly = TRUE)
    
    # ComBat correction with labels (supervised)
    # Step 1: Correct batch effects within training data while preserving biological signal
    dat_corrected <- ComBat(dat, batch=batch, mod=model.matrix(~group))
    
    # Step 2: Adjust test data to match corrected training distribution
    combined_dat <- cbind(dat_corrected, dat_test)
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                       rep(test_batch_id, ncol(dat_test)))
    
    combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                             mod=NULL, ref.batch=ref_batch_id)
    
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "mnn") {
    library(batchelor, quietly = TRUE)
    library(SummarizedExperiment, quietly = TRUE)
    
    # MNN without pre-centering
    combined_dat <- cbind(dat, dat_test)
    # Test set gets a unique batch ID
    test_id <- "TEST_SET"
    combined_batch <- c(batch, rep(test_id, ncol(dat_test)))
    
    # Determine merge order: training batches by size, test last
    u_batches <- unique(batch)
    train_sizes <- table(batch)[u_batches]
    train_ord <- order(train_sizes, decreasing = TRUE)
    merge_ord <- c(u_batches[train_ord], test_id)
    
    # Apply MNN correction
    mnn_object <- batchelor::mnnCorrect(
      combined_dat, 
      batch = combined_batch, 
      merge.order = merge_ord
    )
    mnn_matrix <- SummarizedExperiment::assay(mnn_object)
    
    # Split back
    dat_corrected <- mnn_matrix[, 1:ncol(dat)]
    dat_test_corrected <- mnn_matrix[, (ncol(dat) + 1):ncol(mnn_matrix)]
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "fast_mnn") {
    library(batchelor, quietly = TRUE)
    library(SummarizedExperiment, quietly = TRUE)
    
    # FastMNN correction using batchelor::fastMNN
    cat("Applying FastMNN correction...\n")
    
    # Prepare data for fastMNN - it expects matrices with same number of rows
    # Split training data by batch
    batch_names <- unique(batch)
    batch_matrices <- list()
    
    for (b in batch_names) {
      batch_idx <- which(batch == b)
      batch_matrices[[paste0("batch_", b)]] <- dat[, batch_idx, drop = FALSE]
    }
    
    # Add test data as a separate batch
    batch_matrices[["test_batch"]] <- dat_test
    
    cat(sprintf("FastMNN: Processing %d training batches + 1 test batch\n", length(batch_names)))
    
    # Apply fastMNN - pass matrices directly, not as a list
    fastmnn_result <- do.call(batchelor::fastMNN, c(batch_matrices, list(verbose = FALSE)))
    corrected_matrix <- SummarizedExperiment::assay(fastmnn_result, "corrected")
    
    # Split corrected data back into training and test
    n_train_samples <- ncol(dat)
    dat_corrected <- corrected_matrix[, 1:n_train_samples]
    dat_test_corrected <- corrected_matrix[, (n_train_samples + 1):ncol(corrected_matrix)]
    
    cat("FastMNN correction completed successfully\n")
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "ruvr") {
    # RUVr: Remove Unwanted Variation using Residuals
    # Custom implementation without ruv package dependency
    cat("Applying RUVr correction...\n")
    
    # Step 1: Fit initial GLM on training data to get residuals
    design <- model.matrix(~ group + batch)
    
    # Fit gene-wise linear models
    cat("Fitting initial GLM to estimate residuals...\n")
    residuals <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))
    for (i in 1:nrow(dat)) {
      fit <- lm(dat[i, ] ~ group + batch)
      residuals[i, ] <- residuals(fit)
    }
    
    # Step 2: Estimate unwanted variation factors from residuals using SVD
    k <- 3  # Number of unwanted variation factors
    cat(sprintf("Estimating %d unwanted variation factors...\n", k))
    
    svd_res <- svd(residuals)
    W <- svd_res$u[, 1:k, drop = FALSE]  # Factor loadings (genes x k)
    alpha <- svd_res$v[, 1:k, drop = FALSE] %*% diag(svd_res$d[1:k])  # Factor scores (samples x k)
    
    # Step 3: Correct training data by regressing out the factors
    dat_corrected <- dat
    for (i in 1:nrow(dat)) {
      fit <- lm(dat[i, ] ~ alpha)
      dat_corrected[i, ] <- residuals(fit) + mean(dat[i, ])
    }
    
    # Step 4: Project test data onto the learned factors and correct
    cat("Projecting test data onto learned factors...\n")
    alpha_test <- t(dat_test) %*% W
    
    dat_test_corrected <- dat_test
    for (i in 1:nrow(dat_test)) {
      fit <- lm(dat_test[i, ] ~ alpha_test)
      dat_test_corrected[i, ] <- residuals(fit) + mean(dat_test[i, ])
    }
    
    cat("RUVr correction complete\n")
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "gmm") {
    # GMM adjustment: fits 2-component GMM to each gene within each batch
    cat("Applying GMM adjustment...\n")
    
    # Source the GMM adjustment function (use absolute path from workspace root)
    gmm_script <- file.path(getwd(), "..", "..", "adjust", "gmm_adjust.R")
    if (!file.exists(gmm_script)) {
      gmm_script <- "scripts/adjust/gmm_adjust.R"  # Fallback to relative from workspace root
    }
    source(gmm_script)
    
    # Apply GMM to training data
    dat_corrected <- gmm_adjust(
      data = dat,
      batch = batch,
      genes_are_columns = FALSE,
      mean_mean_zero = TRUE,
      unit_var = TRUE,
      log_transform = FALSE,
      debug = FALSE,
      num_workers = 1
    )
    
    # Apply GMM to test data (single batch)
    dat_test_corrected <- gmm_adjust(
      data = dat_test,
      batch = rep(1, ncol(dat_test)),
      genes_are_columns = FALSE,
      mean_mean_zero = TRUE,
      unit_var = TRUE,
      log_transform = FALSE,
      debug = FALSE,
      num_workers = 1
    )
    
    cat("GMM adjustment complete\n")
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (startsWith(method, "pace_")) {
    # PACE: Pathway-Aware Consensus Estimator with different parameter settings
    pace_variant <- sub("pace_", "", method)
    cat(sprintf("Applying PACE adjustment (%s variant)...\n", pace_variant))
    
    # Load reticulate for Python integration
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("reticulate package is required for PACE method. Please install it.")
    }
    library(reticulate, quietly = TRUE)
    
    # Source the PACE wrapper
    pace_wrapper_path <- "scripts/pace_wrapper.py"
    if (!file.exists(pace_wrapper_path)) {
      stop(sprintf("PACE wrapper not found at: %s", pace_wrapper_path))
    }
    
    tryCatch({
      # Import the PACE wrapper module
      source_python(pace_wrapper_path)
      
      # Apply PACE correction
      pace_result <- apply_pace_correction(
        train_data = dat,
        test_data = dat_test,
        train_batch = batch,
        pace_variant = pace_variant,
        gene_names = rownames(dat)
      )
      
      dat_corrected <- pace_result$train_corrected
      dat_test_corrected <- pace_result$test_corrected
      
      # Log detailed PACE metrics if available
      if (!is.null(pace_result$metrics)) {
        metrics <- pace_result$metrics
        cat(sprintf("PACE %s detailed metrics:\n", pace_variant))
        cat(sprintf("  S_factor=%.4f\n", ifelse(is.null(metrics$S_factor), 0, metrics$S_factor)))
        cat(sprintf("  Alpha: mean=%.4f, std=%.4f, min=%.4f, max=%.4f\n",
                    ifelse(is.null(metrics$alpha_mean), 0, metrics$alpha_mean),
                    ifelse(is.null(metrics$alpha_std), 0, metrics$alpha_std),
                    ifelse(is.null(metrics$alpha_min), 0, metrics$alpha_min),
                    ifelse(is.null(metrics$alpha_max), 0, metrics$alpha_max)))
        cat(sprintf("  Beta: mean=%.4f, std=%.4f\n",
                    ifelse(is.null(metrics$beta_mean), 0, metrics$beta_mean),
                    ifelse(is.null(metrics$beta_std), 0, metrics$beta_std)))
        cat(sprintf("  Genes: common=%d, unique=%d\n",
                    ifelse(is.null(metrics$n_common_genes), 0, metrics$n_common_genes),
                    ifelse(is.null(metrics$n_unique_genes), 0, metrics$n_unique_genes)))
        cat(sprintf("  Pathways used=%d, prior_strength=%.1f\n",
                    ifelse(is.null(metrics$n_pathways_used), 0, metrics$n_pathways_used),
                    ifelse(is.null(metrics$prior_strength), 1, metrics$prior_strength)))
        
        # Diagnostic interpretation
        alpha_mean <- ifelse(is.null(metrics$alpha_mean), 0, metrics$alpha_mean)
        alpha_std <- ifelse(is.null(metrics$alpha_std), 0, metrics$alpha_std)
        if (alpha_mean < 0.2 && alpha_std < 0.1) {
          cat("  [DIAGNOSTIC] Algorithm appears TOO CONSERVATIVE - alphas near minimum clamp\n")
        } else if (alpha_mean > 2.0) {
          cat("  [DIAGNOSTIC] Algorithm appears AGGRESSIVE - high variance corrections\n")
        } else {
          cat("  [DIAGNOSTIC] Algorithm appears BALANCED\n")
        }
      }
      
      cat(sprintf("PACE %s adjustment complete\n", pace_variant))
      
    }, error = function(e) {
      cat(sprintf("[WARNING] PACE %s correction failed: %s\n", pace_variant, e$message))
      cat("[WARNING] Falling back to unadjusted data\n")
      dat_corrected <- dat
      dat_test_corrected <- dat_test
    })
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  stop(sprintf("Unknown batch correction method: %s", method))
}

#' Global scaling: scale entire dataset to have overall variance = 1
#' Preserves relative gene importance while putting data on consistent scale
#' @param dat_train Training data matrix
#' @param dat_test Test data matrix
#' @return List with scaled training and test data
global_scale <- function(dat_train, dat_test) {
  # Compute global mean and SD from training data
  train_mean <- mean(dat_train, na.rm = TRUE)
  train_sd <- sd(as.vector(dat_train), na.rm = TRUE)
  
  # Check for invalid scaling parameters
  if (is.na(train_mean) || is.na(train_sd) || train_sd == 0 || !is.finite(train_sd) || train_sd < 1e-10) {
    cat(sprintf("Global scaling: mean=%s, sd=%s\n", 
                ifelse(is.na(train_mean), "NA", sprintf("%.4f", train_mean)),
                ifelse(is.na(train_sd), "NA", sprintf("%.4f", train_sd))))
    
    # If standard deviation is too small, use unit scaling (subtract mean only)
    if (!is.na(train_mean) && is.finite(train_mean)) {
      cat("[WARNING] Using unit scaling (mean centering only) due to low variance\n")
      dat_train_scaled <- dat_train - train_mean
      dat_test_scaled <- dat_test - train_mean
      
      return(list(
        dat_train = dat_train_scaled,
        dat_test = dat_test_scaled
      ))
    } else {
      stop("Scaling produced invalid values (NA, 0, or infinite standard deviation)")
    }
  }
  
  # Apply same transformation to both train and test
  dat_train_scaled <- (dat_train - train_mean) / train_sd
  dat_test_scaled <- (dat_test - train_mean) / train_sd
  
  # Check for NaN/Inf in results
  if (any(!is.finite(dat_train_scaled)) || any(!is.finite(dat_test_scaled))) {
    stop("Scaling produced non-finite values in output data")
  }
  
  cat(sprintf("Global scaling: mean=%.4f, sd=%.4f\n", train_mean, train_sd))
  
  list(
    dat_train = dat_train_scaled,
    dat_test = dat_test_scaled
  )
}

# ====================================================================
# DIMENSIONALITY REDUCTION (Single Responsibility: Transformation)
# ====================================================================

#' Perform PCA on training data and project test data
#' @param dat_train Training data matrix (genes x samples)
#' @param dat_test Test data matrix (genes x samples)
#' @return List with PCA results
compute_pca <- function(dat_train, dat_test) {
  cat("Computing PCA (fit on training, project test)...\n")
  
  # Fit PCA on training data only
  # Transpose: PCA expects samples x features
  pca_fit <- prcomp(t(dat_train), center = TRUE, scale. = TRUE)
  
  # Project training data
  train_coords <- pca_fit$x[, 1:2]
  
  # Project test data using the same transformation
  test_coords <- predict(pca_fit, newdata = t(dat_test))[, 1:2]
  
  # Combine coordinates
  combined_coords <- rbind(train_coords, test_coords)
  
  # Calculate variance explained
  var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  
  list(
    coords = combined_coords,
    var_explained = var_explained[1:2],
    method = "PCA"
  )
}

#' Perform LDA on training data and project test data
#' @param dat_train Training data matrix (genes x samples)
#' @param dat_test Test data matrix (genes x samples)
#' @param labels_train Training sample labels
#' @param batch_train Training batch assignments
#' @return List with LDA results
compute_lda <- function(dat_train, dat_test, labels_train, batch_train) {
  cat("Computing LDA (fit on training, project test)...\n")
  
  library(MASS, quietly = TRUE)
  
  # Transpose: LDA expects samples x features
  dat_train_t <- t(dat_train)
  dat_test_t <- t(dat_test)
  
  # Create combined grouping for training: label + batch
  # This creates groups like "0_1", "0_2", "1_1", "1_2" etc.
  combined_group_train <- paste(labels_train, batch_train, sep = "_")
  
  # LDA requires at least 2 classes
  if (length(unique(combined_group_train)) < 2) {
    warning("LDA requires at least 2 classes. Skipping.")
    return(NULL)
  }
  
  # Fit LDA on training data only
  lda_fit <- lda(dat_train_t, grouping = as.factor(combined_group_train))
  
  # Project training data
  train_coords <- predict(lda_fit, dat_train_t)$x
  
  # Project test data using the same transformation
  test_coords <- predict(lda_fit, dat_test_t)$x
  
  # Combine coordinates
  combined_coords <- rbind(train_coords, test_coords)
  
  # Handle 1D case (only 1 discriminant)
  if (is.null(dim(combined_coords))) {
    combined_coords <- cbind(combined_coords, rep(0, length(combined_coords)))
    colnames(combined_coords) <- c("LD1", "LD2")
  } else if (ncol(combined_coords) == 1) {
    combined_coords <- cbind(combined_coords, rep(0, nrow(combined_coords)))
    colnames(combined_coords) <- c("LD1", "LD2")
  }
  
  list(
    coords = combined_coords[, 1:2],
    method = "LDA"
  )
}

#' Perform UMAP on data
#' @param dat Data matrix (genes x samples)
#' @return List with UMAP results
compute_umap <- function(dat) {
  cat("Computing UMAP...\n")
  
  library(umap, quietly = TRUE)
  
  # Use fixed seed for reproducibility
  set.seed(42)
  
  # Transpose: UMAP expects samples x features
  # Use custom config for stability
  custom_config <- umap.defaults
  custom_config$random_state <- 42
  custom_config$n_neighbors <- min(15, ncol(dat) - 1)
  
  umap_result <- umap(t(dat), config = custom_config)
  
  list(
    coords = umap_result$layout,
    method = "UMAP"
  )
}


# ====================================================================
# VISUALIZATION (Single Responsibility: Plotting)
# ====================================================================

#' Create a dimensionality reduction plot
#' @param coords 2D coordinates matrix
#' @param batch Batch assignments
#' @param labels Sample labels
#' @param method Method name (PCA, LDA, UMAP)
#' @param var_explained Variance explained (for PCA)
#' @param title Plot title
#' @return ggplot object
create_reduction_plot <- function(coords, batch, labels, dataset_type, method, 
                                 var_explained = NULL, title = NULL) {
  
  # Create data frame for plotting
  plot_df <- data.frame(
    Dim1 = coords[, 1],
    Dim2 = coords[, 2],
    Batch = as.factor(batch),
    Label = as.factor(labels),
    DatasetType = as.factor(dataset_type)
  )
  
  # Calculate sample counts per batch
  batch_counts <- table(batch)
  
  # Create new labels with sample counts: "StudyName (n=XX)"
  batch_labels <- sapply(levels(plot_df$Batch), function(b) {
    sprintf("%s (n=%d)", b, batch_counts[b])
  })
  names(batch_labels) <- levels(plot_df$Batch)
  
  # Determine axis labels
  if (method == "PCA" && !is.null(var_explained)) {
    xlab <- sprintf("PC1 (%.1f%%)", var_explained[1] * 100)
    ylab <- sprintf("PC2 (%.1f%%)", var_explained[2] * 100)
  } else if (method == "LDA") {
    xlab <- "LD1"
    ylab <- "LD2"
  } else if (method == "UMAP") {
    xlab <- "UMAP1"
    ylab <- "UMAP2"
  } else {
    xlab <- "Dimension 1"
    ylab <- "Dimension 2"
  }
  
  # Create base plot
  p <- ggplot(plot_df, aes(x = Dim1, y = Dim2))
  
  # Add points: color by batch, shape by label, border by dataset type
  # Use larger size for test set to make it stand out
  p <- p + geom_point(aes(color = Batch, shape = Label, size = DatasetType, alpha = DatasetType))
  p <- p + scale_size_manual(values = c("Training" = 2.5, "Test" = 3.2), guide = "none")
  p <- p + scale_alpha_manual(values = c("Training" = 0.6, "Test" = 0.75))
  
  # Styling
  p <- p + theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold")
    ) +
    labs(
      x = xlab,
      y = ylab,
      title = title,
      color = "Study",
      shape = "Status"
    ) +
    scale_shape_manual(values = c(16, 17), labels = c("Non-Progressing", "Progressing")) +
    # Use a colorblind-friendly palette with enough colors
    scale_color_brewer(palette = "Set2", labels = batch_labels)
  
  p
}

#' Save plot to file
#' @param plot ggplot object
#' @param filepath Output file path
#' @param width Plot width in inches
#' @param height Plot height in inches
save_plot <- function(plot, filepath, width = 8, height = 6) {
  ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  cat(sprintf("Saved plot: %s\n", filepath))
}

#' Create and save all three reduction plots
#' @param dat_train Training data matrix
#' @param dat_test Test data matrix
#' @param batch_train Training batch assignments
#' @param labels_train Training sample labels
#' @param labels_test Test sample labels
#' @param test_study Test study name
#' @param adjuster Adjuster name
#' @param num_datasets Number of datasets
#' @param output_dir Output directory
create_all_visualizations <- function(dat_train, dat_test, batch_train, 
                                     labels_train, labels_test, test_study,
                                     adjuster, num_datasets, output_dir) {
  
  # Combine training and test data
  dat_combined <- cbind(dat_train, dat_test)
  
  # Create combined batch vector (test set gets its own label)
  batch_combined <- c(batch_train, rep(test_study, ncol(dat_test)))
  
  # Create combined labels
  labels_combined <- c(labels_train, labels_test)
  
  # Create dataset type indicator (train vs test)
  dataset_type <- c(rep("Training", ncol(dat_train)), rep("Test", ncol(dat_test)))
  
  base_name <- sprintf("%s_n%d_test%s", adjuster, num_datasets, test_study)
  
  # PCA (fit on training, project test)
  pca_result <- compute_pca(dat_train, dat_test)
  pca_plot <- create_reduction_plot(
    coords = pca_result$coords,
    batch = batch_combined,
    labels = labels_combined,
    dataset_type = dataset_type,
    method = "PCA",
    var_explained = pca_result$var_explained,
    title = sprintf("PCA: %s (n=%d, test=%s)", 
                   tools::toTitleCase(adjuster), num_datasets, test_study)
  )
  save_plot(
    pca_plot,
    file.path(output_dir, "pca", paste0(base_name, ".png"))
  )
  
  # LDA (fit on training, project test)
  lda_result <- compute_lda(dat_train, dat_test, labels_train, batch_train)
  if (!is.null(lda_result)) {
    lda_plot <- create_reduction_plot(
      coords = lda_result$coords,
      batch = batch_combined,
      labels = labels_combined,
      dataset_type = dataset_type,
      method = "LDA",
      title = sprintf("LDA: %s (n=%d, test=%s)", 
                     tools::toTitleCase(adjuster), num_datasets, test_study)
    )
    save_plot(
      lda_plot,
      file.path(output_dir, "lda", paste0(base_name, ".png"))
    )
  }
  
  # UMAP (still uses combined data - cannot project separately)
  dat_combined <- cbind(dat_train, dat_test)
  umap_result <- compute_umap(dat_combined)
  umap_plot <- create_reduction_plot(
    coords = umap_result$coords,
    batch = batch_combined,
    labels = labels_combined,
    dataset_type = dataset_type,
    method = "UMAP",
    title = sprintf("UMAP: %s (n=%d, test=%s)", 
                   tools::toTitleCase(adjuster), num_datasets, test_study)
  )
  save_plot(
    umap_plot,
    file.path(output_dir, "umap", paste0(base_name, ".png"))
  )
  
  cat("All visualizations created successfully\n")
}


# ====================================================================
# MAIN EXECUTION
# ====================================================================

main <- function() {
  tryCatch({
    cat("=== BATCH ADJUSTMENT VISUALIZATION ===\n")
    cat(sprintf("Adjuster: %s\n", adjuster))
    cat(sprintf("Num datasets: %d\n", num_datasets))
    cat(sprintf("Test study: %s\n", test_study))
    cat(sprintf("Output directory: %s\n", output_dir))
    cat("======================================\n\n")
    
    # Load and prepare data
    cat("Step 1: Loading data...\n")
    data <- load_and_filter_data(num_datasets)
    
    cat("\nStep 2: Validating test study...\n")
    test_name <- validate_test_study(data$study_names, test_study)
    
    cat("\nStep 3: Preparing train/test split...\n")
    datasets <- prepare_train_test_split(
      data$dat_lst, 
      data$label_lst, 
      test_name, 
      data$study_names
    )
    
    cat("\nStep 4: Reducing features...\n")
    if (reduce == 0) {
      reduced = list(dat = datasets$dat, dat_test = datasets$dat_test)
    }
    else {
      reduced <- reduce_features(datasets$dat, datasets$dat_test, n_genes = reduce)
    }
    
    cat("\nStep 5: Applying batch correction...\n")
    corrected <- apply_batch_correction(
      reduced$dat,
      datasets$batch,
      datasets$group,
      reduced$dat_test,
      adjuster
    )
    
    cat("\nStep 6: Applying global scaling...\n")
    
    # Check if batch correction produced valid data
    if(all(corrected$dat_corrected == corrected$dat_corrected[1,1]) || all(corrected$dat_test_corrected == corrected$dat_test_corrected[1,1])) {
      cat("[WARNING] Batch correction produced constant data - using original data instead\n")
      corrected$dat_corrected <- reduced$dat
      corrected$dat_test_corrected <- reduced$dat_test
    }
    
    scaled <- global_scale(corrected$dat_corrected, corrected$dat_test_corrected)
    
    cat("\nStep 7: Creating visualizations...\n")
    create_all_visualizations(
      dat_train = scaled$dat_train,
      dat_test = scaled$dat_test,
      batch_train = datasets$batch,
      labels_train = datasets$group,
      labels_test = datasets$group_test,
      test_study = test_study,
      adjuster = adjuster,
      num_datasets = num_datasets,
      output_dir = output_dir
    )
    
    cat("\n=== VISUALIZATION COMPLETE ===\n")
    cat(sprintf("Output directory: %s\n", output_dir))
    cat("Files organized by method: pca/, lda/, umap/\n")
    
    return(0)
    
  }, error = function(e) {
    cat(sprintf("\n[ERROR] Visualization failed: %s\n", e$message), file = stderr())
    cat(sprintf("[ERROR] Traceback:\n"), file = stderr())
    traceback_lines <- capture.output(traceback())
    for (line in traceback_lines) {
      cat(sprintf("[ERROR] %s\n", line), file = stderr())
    }
    return(1)
  })
}

# Execute main function
exit_code <- main()
quit(status = exit_code)
