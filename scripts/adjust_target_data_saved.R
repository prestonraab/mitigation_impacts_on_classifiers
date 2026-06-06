#!/usr/bin/env Rscript

# Suppress warnings and messages
options(warn = -1)
options(repos = c(CRAN = "https://cloud.r-project.org"))
# Prevent BiocManager from trying to phone home on compute nodes (no internet)
options(BiocManager.check_repositories = FALSE)

suppressMessages(suppressWarnings({
  required_packages <- c("argparse", "dplyr", "purrr", "sva", "batchelor", 
                        "SummarizedExperiment", "MASS")
  sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
}))

# Configure reticulate to use the pixi environment Python
if (requireNamespace("reticulate", quietly = TRUE)) {
  Sys.setenv(RETICULATE_AUTOCREATE_VENV = "FALSE")
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
  pixi_python <- file.path(getwd(), ".pixi/envs/default/bin/python")
  if (file.exists(pixi_python)) {
    reticulate::use_python(pixi_python, required = TRUE)
  }
}

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Perform batch correction and save adjusted target dataset")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, naive, rank_samples, rank_twice, npn, combat, combat_mean, combat_sup, mnn, fast_mnn, ruvg, yugene, cublock, angel, tdm, rnabc, shambhala2, coconut, rankin, recombat, or recombat_sup")
parser$add_argument("--output-data", type = "character", required = TRUE,
                   help = "Output CSV file path to save the adjusted target matrix")
parser$add_argument("--num-datasets", type = "integer", default = NULL,
                   help = "Number of datasets to use as reference (optional, for SA/UK use 1, otherwise null)")
parser$add_argument("--test-study", type = "character", default = "India",
                   help = "Test/target study name (default: India for backward compatibility)")

args <- parser$parse_args()

valid_adjusters <- c("unadjusted", "naive", "rank_samples", "rank_twice", "npn",
                     "combat", "combat_mean", "combat_sup", "combat_sup_ca", "combat_sup_mg",
                     "combat_sup_nat", "combat_sup_mo",
                     "combat_sup_nat_v2", "combat_sup_mo_v2", "mnn", "fast_mnn",
                     "ruvg", "ruvr", "yugene", "cublock", "angel", "tdm", "rnabc",
                     "shambhala2", "coconut", "rankin", "recombat", "recombat_sup")

if (!args$adjuster %in% valid_adjusters) {
  stop(sprintf("Invalid adjuster '%s'. Must be one of: %s", 
               args$adjuster, paste(valid_adjusters, collapse = ", ")))
}

adjuster <- args$adjuster
output_data <- args$output_data
num_datasets <- args$num_datasets
test_study <- args$test_study

# Create output directories
dir.create(dirname(output_data), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Adjusting with %s (test_study=%s, num_datasets=%s)\n", adjuster, test_study,
           ifelse(is.null(num_datasets), "NULL", num_datasets)))

# Load helper functions
source("scripts/helper.R")

# ====================================================================
# ADJUSTMENT HELPER FUNCTIONS
# ====================================================================

rank_normalized <- function(matrix_, dim) {
  if (dim < 1 || dim > 2) {
    stop("Invalid dimension. Must be 1 for rows or 2 for columns.")
  }
  ranked = apply(matrix_, dim, rank, ties.method = "average")
  if (dim == 1 && is.matrix(ranked)) {
    ranked = t(ranked)
  }
  return(ranked / max(ranked, na.rm = TRUE))
}

adjust_ranked_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  ranked = rank_normalized(matrix_, 1)
  batch_levels <- unique(batch)
  ranked2 <- matrix(NA, nrow = nrow(ranked), ncol = ncol(ranked))
  for (b in batch_levels) {
    batch_indices <- which(batch == b)
    batch_data <- ranked[, batch_indices, drop = FALSE]
    if (ncol(batch_data) > 1) {
      batch_ranked <- rank_normalized(batch_data, 2)
      ranked2[, batch_indices] <- batch_ranked
    } else {
      ranked2[, batch_indices] <- batch_data
    }
  }
  if (any(is.na(ranked2))) {
    ranked2[is.na(ranked2)] <- ranked[is.na(ranked2)]
  }
  max_val <- max(ranked2, na.rm = TRUE)
  if (max_val == 0) max_val <- 1
  return(ranked2 / max_val)
}

adjust_npn <- function(matrix_, batch, debug = FALSE) {
  if (!requireNamespace("huge", quietly = TRUE)) {
    stop("Package 'huge' is required for NPN adjustment but is not installed.")
  }
  if (is.null(batch)) {
    matrix_t <- t(matrix_)
    npn_transformed_t <- huge::huge.npn(matrix_t, verbose = FALSE)
    return(t(npn_transformed_t))
  } else {
    batch_levels <- unique(batch)
    matrix_by_batch <- list()
    for (b in batch_levels) {
      batch_indices <- which(batch == b)
      if (length(batch_indices) > 0) {
        matrix_by_batch[[as.character(b)]] <- matrix_[, batch_indices, drop = FALSE]
      }
    }
    for (b in names(matrix_by_batch)) {
      matrix_t <- t(matrix_by_batch[[b]])
      npn_transformed_t <- huge::huge.npn(matrix_t, verbose = FALSE)
      matrix_by_batch[[b]] <- t(npn_transformed_t)
    }
    result_matrix <- matrix_
    for (b in names(matrix_by_batch)) {
      batch_indices <- which(batch == as.character(b))
      result_matrix[, batch_indices] <- matrix_by_batch[[b]]
    }
    return(result_matrix)
  }
}

adjust_ranked_samples_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  batch_levels <- unique(batch)
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  for (b in batch_levels) {
    batch_indices <- which(batch == b)
    batch_data <- matrix_[, batch_indices, drop = FALSE]
    if (ncol(batch_data) > 1) {
      batch_ranked <- rank_normalized(batch_data, 2)
      result_matrix[, batch_indices] <- batch_ranked
    } else {
      result_matrix[, batch_indices] <- batch_data / max(batch_data, na.rm = TRUE)
    }
  }
  if (any(is.na(result_matrix))) {
    result_matrix[is.na(result_matrix)] <- matrix_[is.na(result_matrix)]
  }
  max_val <- max(result_matrix, na.rm = TRUE)
  if (max_val == 0) max_val <- 1
  return(result_matrix / max_val)
}

adjust_ranked_twice_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  batch_levels <- unique(batch)
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  for (b in batch_levels) {
    batch_indices <- which(batch == b)
    batch_data <- matrix_[, batch_indices, drop = FALSE]
    if (ncol(batch_data) > 1) {
      batch_ranked <- rank_normalized(rank_normalized(batch_data, 1), 2)
      result_matrix[, batch_indices] <- batch_ranked
    } else {
      batch_ranked <- rank_normalized(batch_data, 1)
      result_matrix[, batch_indices] <- batch_ranked
    }
  }
  if (any(is.na(result_matrix))) {
    result_matrix[is.na(result_matrix)] <- matrix_[is.na(result_matrix)]
  }
  max_val <- max(result_matrix, na.rm = TRUE)
  if (max_val == 0) max_val <- 1
  return(result_matrix / max_val)
}

adjust_yugene <- function(matrix_, debug = FALSE) {
  if (!requireNamespace("YuGene", quietly = TRUE)) {
    stop("Package 'YuGene' is required but not installed.")
  }
  result_matrix <- YuGene::YuGene(matrix_)
  result_matrix <- as.matrix(unclass(result_matrix))
  return(result_matrix)
}

adjust_cublock <- function(matrix_, batch, debug = FALSE) {
  input_file <- tempfile(fileext = ".csv")
  output_file <- tempfile(fileext = ".csv")
  write.table(matrix_, input_file, sep=",", row.names=FALSE, col.names=FALSE)
  pixi_cmd <- "pixi"
  user_pixi <- "/home/phr23/.pixi/bin/pixi"
  if (file.exists(user_pixi)) {
    pixi_cmd <- user_pixi
  }
  octave_eval <- sprintf("addpath(\".\"); pkg load statistics; data = csvread(\"%s\"); [norm_data] = CuBlock(data); csvwrite(\"%s\", norm_data);", input_file, output_file)
  sys_res <- system2(pixi_cmd, args = c("run", "octave", "--eval", shQuote(octave_eval)), wait = TRUE)
  if (sys_res != 0) {
    stop("Octave system call failed. Verify Octave is installed and CuBlock.m is in the working directory.")
  }
  result_matrix <- as.matrix(read.csv(output_file, header=FALSE))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  unlink(input_file)
  unlink(output_file)
  return(result_matrix)
}

adjust_angel <- function(matrix_, debug = FALSE) {
  result_matrix <- apply(matrix_, 2, function(x) {
    r <- rank(x, ties.method = "average")
    (r - min(r)) / (max(r) - min(r))
  })
  if (!is.matrix(result_matrix)) {
      result_matrix <- as.matrix(result_matrix)
  }
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  return(result_matrix)
}

adjust_tdm <- function(dat_train, dat_test, debug = FALSE) {
  if (!requireNamespace("TDM", quietly = TRUE)) {
    stop("Package 'TDM' is required but not installed.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for TDM.")
  }
  train_dt <- data.table::as.data.table(dat_train, keep.rownames = "gene")
  test_dt <- data.table::as.data.table(dat_test, keep.rownames = "gene")
  res_dt <- TDM::tdm_transform(target_data = test_dt, ref_data = train_dt)
  gene_names <- res_dt[[1]]
  result_matrix <- as.matrix(res_dt[, -1, with = FALSE])
  rownames(result_matrix) <- gene_names
  return(result_matrix)
}

adjust_rnabc <- function(dat_train, dat_test, debug = FALSE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required but not installed.")
  }
  target_dist <- rowMeans(dat_train)
  combined_for_qnorm <- cbind(dat_test, target_dist)
  qnorm_res <- limma::normalizeQuantiles(as.matrix(combined_for_qnorm))
  dat_test_qnorm <- qnorm_res[, 1:ncol(dat_test), drop = FALSE]
  rownames(dat_test_qnorm) <- rownames(dat_test)
  colnames(dat_test_qnorm) <- colnames(dat_test)
  combined_dat <- cbind(dat_train, dat_test_qnorm)
  batch_vec <- c(rep(1, ncol(dat_train)), rep(2, ncol(dat_test_qnorm)))
  combined_bc <- sva::ComBat(combined_dat, batch = batch_vec, ref.batch = 1)
  dat_test_corrected <- combined_bc[, (ncol(dat_train) + 1):ncol(combined_bc), drop = FALSE]
  return(dat_test_corrected)
}

adjust_shambhala2 <- function(matrix_, calib_P, ref_Q, debug = FALSE) {
  if (!exists("shambhala2_harmonize")) {
    if (file.exists("Shambhala2/Shambhala2.R")) {
      source("Shambhala2/Shambhala2.R")
    } else if (file.exists("Shambhala2.R")) {
      source("Shambhala2.R")
    } else {
      stop("Shambhala-2 functions not found. Source the Shambhala2.R script.")
    }
  }
  result_matrix <- shambhala2_harmonize(matrix_, P_dataset = calib_P, Q_dataset = ref_Q)
  return(result_matrix)
}

adjust_coconut <- function(matrix_, batch, group, debug = FALSE) {
  if (!requireNamespace("COCONUT", quietly = TRUE)) {
    stop("Package 'COCONUT' is required but not installed.")
  }
  gse_list <- list()
  for (b in unique(batch)) {
    idx <- which(batch == b)
    disease_vec <- as.numeric(group[idx])
    pheno_df <- data.frame(
      disease_state = disease_vec,
      dummy = 1,
      row.names = colnames(matrix_[, idx])
    )
    gse_list[[as.character(b)]] <- list(
      pheno = pheno_df,
      genes = matrix_[, idx, drop = FALSE]
    )
  }
  res <- COCONUT::COCONUT(GSEs = gse_list, control.0.col = "disease_state")
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  for (b in names(res$COCONUTList)) {
    disease_cols <- colnames(res$COCONUTList[[b]]$genes)
    result_matrix[, disease_cols] <- as.matrix(res$COCONUTList[[b]]$genes)
  }
  for (b in names(res$controlList$GSEs)) {
    control_cols <- colnames(res$controlList$GSEs[[b]]$genes)
    result_matrix[, control_cols] <- as.matrix(res$controlList$GSEs[[b]]$genes)
  }
  return(result_matrix)
}

adjust_rankin <- function(matrix_, n_svd = 1L, train_indices = NULL, debug = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for Rank-In.")
  }
  rankin_script <- local({
    candidates <- c(
      "scripts/rankin.py",
      "rankin.py",
      "/home/phr23/confounded_analysis/scripts/evaluations/book_chapter/scripts/rankin.py"
    )
    found <- Filter(file.exists, candidates)
    if (length(found) == 0) {
      stop(sprintf("Rank-In python script not found. Searched: %s", paste(candidates, collapse = ", ")))
    }
    found[[1]]
  })
  reticulate::source_python(rankin_script)
  result_list <- rank_in_from_r(
    unname(as.list(as.data.frame(t(matrix_)))), 
    train_indices = train_indices,
    n_svd = as.integer(n_svd)
  )
  result_matrix <- matrix(unlist(result_list), nrow = nrow(matrix_), ncol = ncol(matrix_), byrow = TRUE)
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  return(result_matrix)
}

adjust_recombat <- function(matrix_, batch, debug = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }
  data_t <- t(matrix_)
  pd <- reticulate::import("pandas", convert = FALSE)
  recombat_pkg <- reticulate::import("reComBat", convert = FALSE)
  data_pd <- pd$DataFrame(data_t)
  batch_pd <- pd$Series(batch)
  combat_model <- recombat_pkg$reComBat(parametric = TRUE, model = "elastic_net")
  res_pd <- combat_model$fit_transform(data_pd, batch_pd)
  res_matrix_t <- reticulate::py_to_r(res_pd)
  result_matrix <- t(as.matrix(res_matrix_t))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  return(result_matrix)
}

adjust_recombat_sup <- function(matrix_, batch, group, debug = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }
  data_t <- t(matrix_)
  pd <- reticulate::import("pandas", convert = FALSE)
  recombat_pkg <- reticulate::import("reComBat", convert = FALSE)
  data_pd <- pd$DataFrame(data_t)
  batch_pd <- pd$Series(batch)
  X_pd <- pd$DataFrame(list(disease = as.integer(group)))
  combat_model <- recombat_pkg$reComBat(parametric = TRUE, model = "elastic_net")
  res_pd <- combat_model$fit_transform(data_pd, batch_pd, X = X_pd)
  res_matrix_t <- reticulate::py_to_r(res_pd)
  result_matrix <- t(as.matrix(res_matrix_t))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  return(result_matrix)
}

# ====================================================================
# MAIN BATCH CORRECTION WRAPPER
# ====================================================================

apply_batch_correction <- function(dat, batch, group, dat_test, method, group_test) {
  if (method == "unadjusted") {
    return(list(dat_corrected = dat, dat_test_corrected = dat_test))
    
  } else if (method == "naive") {
    dat_corrected <- dat
    unique_batches <- unique(batch)
    overall_mean <- rowMeans(dat)
    overall_var <- apply(dat, 1, var)
    for (b in unique_batches) {
      batch_idx <- which(batch == b)
      batch_data <- dat[, batch_idx, drop = FALSE]
      batch_mean <- rowMeans(batch_data)
      batch_var <- apply(batch_data, 1, var)
      batch_sd <- sqrt(pmax(batch_var, 1e-10))
      overall_sd <- sqrt(pmax(overall_var, 1e-10))
      for (i in 1:nrow(dat)) {
        if (batch_sd[i] > 1e-10) {
          dat_corrected[i, batch_idx] <- (batch_data[i, ] - batch_mean[i]) / batch_sd[i] * overall_sd[i] + overall_mean[i]
        } else {
          dat_corrected[i, batch_idx] <- overall_mean[i]
        }
      }
    }
    test_mean <- rowMeans(dat_test)
    test_var <- apply(dat_test, 1, var)
    test_sd <- sqrt(pmax(test_var, 1e-10))
    dat_test_corrected <- dat_test
    overall_sd <- sqrt(pmax(overall_var, 1e-10))
    for (i in 1:nrow(dat_test)) {
      if (test_sd[i] > 1e-10) {
        dat_test_corrected[i, ] <- (dat_test[i, ] - test_mean[i]) / test_sd[i] * overall_sd[i] + overall_mean[i]
      } else {
        dat_test_corrected[i, ] <- overall_mean[i]
      }
    }
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "rank_samples") {
    dat_corrected <- adjust_ranked_samples_with_batch_info(dat, batch, debug = FALSE)
    test_batch <- rep(1, ncol(dat_test))
    dat_test_corrected <- adjust_ranked_samples_with_batch_info(dat_test, test_batch, debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "rank_twice") {
    dat_corrected <- adjust_ranked_twice_with_batch_info(dat, batch, debug = FALSE)
    test_batch <- rep(1, ncol(dat_test))
    dat_test_corrected <- adjust_ranked_twice_with_batch_info(dat_test, test_batch, debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "npn") {
    dat_corrected <- adjust_npn(dat, batch, debug = FALSE)
    test_batch <- rep(1, ncol(dat_test))
    dat_test_corrected <- adjust_npn(dat_test, test_batch, debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat") {
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- ComBat(dat, batch=batch, mod=NULL)
    }
    combined_dat <- cbind(dat_corrected, dat_test)
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), rep(test_batch_id, ncol(dat_test)))
    combat_combined <- ComBat(combined_dat, batch=combined_batch, mod=NULL, ref.batch=ref_batch_id)
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat_mean") {
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- ComBat(dat, batch=batch, mod=NULL, mean.only=TRUE)
    }
    combined_dat <- cbind(dat_corrected, dat_test)
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), rep(test_batch_id, ncol(dat_test)))
    combat_combined <- ComBat(combined_dat, batch=combined_batch, mod=NULL, ref.batch=ref_batch_id, mean.only=TRUE)
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat_sup") {
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- ComBat(dat, batch=batch, mod=model.matrix(~group))
    }
    combined_dat <- cbind(dat_corrected, dat_test)
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), rep(test_batch_id, ncol(dat_test)))
    # Step 2 is unsupervised: test labels are unavailable at correction time.
    combat_combined <- ComBat(combined_dat, batch=combined_batch, mod=NULL, ref.batch=ref_batch_id)
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat_sup_ca") {
    # Supervised ComBat on TRAINING (protects class signal during batch estimation),
    # then project the TEST batch against CLASS-AGNOSTIC reference statistics.
    #
    # The original combat_sup aligns the test batch to the frozen, class-bimodal
    # training distribution — without test labels it cannot know which mode each
    # test sample belongs to, so samples scatter into the wrong class cluster and
    # their variance is shrunk to the tiny within-class residual (KNN purity ~0.19).
    #
    # Here we instead align the test batch's marginal mean and scale to the
    # training's grand mean (class-marginalized) and pooled SD (which retains the
    # full between-class spread). A single per-gene affine map preserves the test
    # set's own internal class structure rather than forcing it onto training modes.
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- ComBat(dat, batch=batch, mod=model.matrix(~group))
    }
    # Class-agnostic reference statistics from the corrected training data
    ref_mean <- rowMeans(dat_corrected)
    ref_sd   <- sqrt(pmax(apply(dat_corrected, 1, var), 1e-10))
    # Test batch's own marginal statistics (no labels used)
    test_mean <- rowMeans(dat_test)
    test_sd   <- sqrt(pmax(apply(dat_test, 1, var), 1e-10))
    # Per-gene affine alignment of test marginal to training class-agnostic marginal
    dat_test_corrected <- (dat_test - test_mean) / test_sd * ref_sd + ref_mean
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat_sup_mg") {
    # Supervised ComBat on TRAINING + mean-only per-gene step-2 + global scale.
    #
    # Problem with combat_sup step-2: per-gene delta_g = var_test_g / var_train_g
    # where var_train_g is the tiny within-class residual after supervised correction.
    # Large delta_g values for class-associated genes reweight the class axis enough
    # to invert the test class ordering (verified empirically: separation = -5.65
    # vs +24.5 for plain combat).
    #
    # Mean-only step-2 removes per-gene scale correction entirely, preserving the
    # class ordering for all batch types. But alone it can't handle platform-level
    # dynamic range differences (e.g. microarray counts ~0-65000, RNA-seq log-counts
    # live in a different range entirely).
    #
    # Global scale factor: compute a single positive scalar
    #   k = sd(vec(dat_corrected)) / sd(vec(dat_test_step2))
    # and apply it uniformly. A positive scalar cannot invert any direction, so
    # class ordering is guaranteed to survive, while cross-platform scale differences
    # are corrected.
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- ComBat(dat, batch=batch, mod=model.matrix(~group))
    }
    # Step 2a: mean-only per-gene alignment (removes location shift, no delta)
    combined_dat   <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test)))
    step2 <- ComBat(combined_dat, batch=combined_batch, mod=NULL,
                    ref.batch=1L, mean.only=TRUE)
    dat_test_step2 <- step2[, (ncol(dat_corrected) + 1):ncol(step2), drop=FALSE]
    # Step 2b: global scale correction
    ref_global_sd  <- sd(as.vector(dat_corrected))
    test_global_sd <- sd(as.vector(dat_test_step2))
    if (ref_global_sd > 1e-10 && test_global_sd > 1e-10) {
      dat_test_corrected <- dat_test_step2 * (ref_global_sd / test_global_sd)
    } else {
      dat_test_corrected <- dat_test_step2
    }
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat_sup_nat") {
    # Class-protected batch MEAN correction, no variance modification.
    #
    # The ComBat_nat (var.pooled decoupling) approach inflates training variance by
    # sqrt(var.pooled.nat / var.pooled) ~ 5-7x, creating a scale mismatch with the
    # natural-scale test data that breaks step-2 projection even more severely.
    #
    # This version instead uses the correct formulation of the hypothesis: use
    # class-protected OLS to estimate batch means (so class signal is not absorbed
    # as batch when batches are confounded), then apply ONLY the mean shift —
    # identical to limma::removeBatchEffect with class covariate.
    # Training and test then share the same natural variance scale, and step-2
    # mean.only alignment should preserve the test class ordering.
    if (!requireNamespace("limma", quietly = TRUE)) stop("limma required for combat_sup_nat")
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- limma::removeBatchEffect(
        dat,
        batch  = batch,
        design = model.matrix(~group)
      )
    }
    combined_dat   <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test)))
    step2 <- ComBat(combined_dat, batch = combined_batch, mod = NULL,
                    ref.batch = 1L, mean.only = TRUE)
    dat_test_corrected <- step2[, (ncol(dat_corrected) + 1):ncol(step2), drop = FALSE]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat_sup_mo") {
    # Supervised ComBat mean.only=TRUE for step-1 (EB-shrunk class-protected batch
    # mean correction, no delta/variance correction), then mean.only step-2 for test.
    # Nearly identical to combat_sup_nat (limma) but uses ComBat's EB shrinkage on
    # gamma rather than raw OLS; difference is small with 2-3 batches.
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- ComBat(dat, batch = batch,
                              mod = model.matrix(~group), mean.only = TRUE)
    }
    combined_dat   <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test)))
    step2 <- ComBat(combined_dat, batch = combined_batch, mod = NULL,
                    ref.batch = 1L, mean.only = TRUE)
    dat_test_corrected <- step2[, (ncol(dat_corrected) + 1):ncol(step2), drop = FALSE]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat_sup_nat_v2") {
    # Like combat_sup_nat but step-2 uses full ComBat (mean + variance correction).
    # Step-1 (limma mean-only) preserves natural variance in training, so step-2
    # sees comparable scales in both batches and delta.star is well-conditioned.
    if (!requireNamespace("limma", quietly = TRUE)) stop("limma required for combat_sup_nat_v2")
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- limma::removeBatchEffect(
        dat,
        batch  = batch,
        design = model.matrix(~group)
      )
    }
    combined_dat   <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test)))
    step2 <- ComBat(combined_dat, batch = combined_batch, mod = NULL,
                    ref.batch = 1L, mean.only = FALSE)
    dat_test_corrected <- step2[, (ncol(dat_corrected) + 1):ncol(step2), drop = FALSE]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "combat_sup_mo_v2") {
    # Like combat_sup_mo but step-2 uses full ComBat (mean + variance correction).
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- ComBat(dat, batch = batch,
                              mod = model.matrix(~group), mean.only = TRUE)
    }
    combined_dat   <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test)))
    step2 <- ComBat(combined_dat, batch = combined_batch, mod = NULL,
                    ref.batch = 1L, mean.only = FALSE)
    dat_test_corrected <- step2[, (ncol(dat_corrected) + 1):ncol(step2), drop = FALSE]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "mnn") {
    combined_dat <- cbind(dat, dat_test)
    test_id <- max(as.numeric(factor(batch))) + 1
    combined_batch <- c(as.numeric(factor(batch)), rep(test_id, ncol(dat_test)))
    u_batches <- sort(unique(as.numeric(factor(batch))))
    train_sizes <- table(as.numeric(factor(batch)))[as.character(u_batches)]
    train_ord <- order(train_sizes, decreasing = TRUE)
    merge_ord <- c(u_batches[train_ord], test_id)
    mnn_object <- batchelor::mnnCorrect(combined_dat, batch = combined_batch, merge.order = merge_ord)
    mnn_matrix <- SummarizedExperiment::assay(mnn_object)
    dat_corrected <- mnn_matrix[, 1:ncol(dat)]
    dat_test_corrected <- mnn_matrix[, (ncol(dat) + 1):ncol(mnn_matrix)]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "fast_mnn") {
    batch_names <- unique(batch)
    batch_matrices <- list()
    for (b in batch_names) {
      batch_idx <- which(batch == b)
      batch_matrices[[length(batch_matrices) + 1]] <- dat[, batch_idx, drop = FALSE]
    }
    batch_matrices[[length(batch_matrices) + 1]] <- dat_test
    fastmnn_result <- do.call(batchelor::fastMNN, batch_matrices)
    corrected_matrix <- SummarizedExperiment::assay(fastmnn_result, "reconstructed")
    corrected_matrix <- as.matrix(corrected_matrix)
    n_train_samples <- ncol(dat)
    dat_corrected <- corrected_matrix[, 1:n_train_samples]
    dat_test_corrected <- corrected_matrix[, (n_train_samples + 1):ncol(corrected_matrix)]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "ruvg") {
    if (!requireNamespace("RUVSeq", quietly = TRUE)) {
      stop("Package 'RUVSeq' is required but not installed.")
    }
    housekeeping_genes <- c("PPIA", "YWHAZ", "HPRT1", "RPLP0")
    available_hk <- intersect(housekeeping_genes, rownames(dat))
    if (length(available_hk) < 2) {
      cat("WARNING: Fewer than 2 housekeeping genes found. Falling back to unadjusted data.\n")
      return(list(dat_corrected = dat, dat_test_corrected = dat_test))
    }
    k <- min(3, length(available_hk) - 1)
    ruv_res <- RUVSeq::RUVg(as.matrix(dat), available_hk, k=k, isLog=TRUE)
    dat_corrected <- ruv_res$normalizedCounts
    W <- ruv_res$W
    gene_means <- rowMeans(dat)
    dat_centered <- dat - gene_means
    alpha <- solve(t(W) %*% W) %*% t(W) %*% t(dat_centered)
    alpha_hk <- alpha[, which(rownames(dat) %in% available_hk), drop=FALSE]
    dat_test_hk_centered <- dat_test[available_hk, , drop=FALSE] - gene_means[available_hk]
    W_test <- t(dat_test_hk_centered) %*% t(alpha_hk) %*% solve(alpha_hk %*% t(alpha_hk))
    dat_test_centered <- dat_test - gene_means
    dat_test_corrected <- dat_test_centered - t(W_test %*% alpha)
    dat_test_corrected <- dat_test_corrected + gene_means
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "ruvr") {
    if (!requireNamespace("RUVSeq", quietly = TRUE)) {
      stop("Package 'RUVSeq' is required but not installed.")
    }
    # Fit initial model to get residuals using group/labels
    design <- model.matrix(~ group)
    residuals <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))
    for (i in 1:nrow(dat)) {
      fit <- lm(dat[i, ] ~ group)
      residuals[i, ] <- residuals(fit)
    }
    k <- 3
    ruv_res <- RUVSeq::RUVr(as.matrix(dat), c(1:nrow(dat)), k=k, residuals=residuals, isLog=TRUE)
    dat_corrected <- ruv_res$normalizedCounts
    W <- ruv_res$W
    gene_means <- rowMeans(dat)
    dat_centered <- dat - gene_means
    alpha <- solve(t(W) %*% W) %*% t(W) %*% t(dat_centered)
    dat_test_centered <- dat_test - gene_means
    W_test <- t(dat_test_centered) %*% t(alpha) %*% solve(alpha %*% t(alpha))
    dat_test_corrected <- dat_test_centered - t(W_test %*% alpha)
    dat_test_corrected <- dat_test_corrected + gene_means
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "yugene") {
    dat_corrected <- adjust_yugene(dat, debug = FALSE)
    dat_test_corrected <- adjust_yugene(dat_test, debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "cublock") {
    dat_corrected <- adjust_cublock(dat, batch, debug = FALSE)
    dat_test_corrected <- adjust_cublock(dat_test, rep(1, ncol(dat_test)), debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "angel") {
    dat_corrected <- adjust_angel(dat, debug = FALSE)
    dat_test_corrected <- adjust_angel(dat_test, debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "tdm") {
    dat_corrected <- dat
    dat_test_corrected <- adjust_tdm(dat_train = dat, dat_test = dat_test, debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "rnabc") {
    dat_corrected <- dat
    dat_test_corrected <- adjust_rnabc(dat_train = dat, dat_test = dat_test, debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "shambhala2") {
    P_ref_target <- dat
    dat_corrected <- adjust_shambhala2(dat, P_ref_target, P_ref_target, debug = FALSE)
    dat_test_corrected <- adjust_shambhala2(dat_test, P_ref_target, P_ref_target, debug = FALSE)
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "recombat") {
    combined_dat <- cbind(dat, dat_test)
    combined_batch <- c(rep(1, ncol(dat)), rep(2, ncol(dat_test)))
    combined_corrected <- adjust_recombat(combined_dat, combined_batch, debug = FALSE)
    dat_corrected <- combined_corrected[, 1:ncol(dat), drop = FALSE]
    dat_test_corrected <- combined_corrected[, (ncol(dat) + 1):ncol(combined_corrected), drop = FALSE]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "recombat_sup") {
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- adjust_recombat_sup(dat, batch, group, debug = FALSE)
    }
    combined_dat <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test)))
    combat_combined <- ComBat(combined_dat, batch = combined_batch, mod = NULL, ref.batch = 1L)
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined), drop = FALSE]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "coconut") {
    if (length(unique(batch)) < 2) {
      dat_corrected <- dat
    } else {
      dat_corrected <- adjust_coconut(dat, batch, group, debug = FALSE)
    }
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_dat <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), rep(test_batch_id, ncol(dat_test)))
    combined_group <- c(group, group_test)
    combat_combined <- ComBat(combined_dat, batch=combined_batch, mod=model.matrix(~combined_group), ref.batch=ref_batch_id)
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else if (method == "rankin") {
    combined_dat <- cbind(dat, dat_test)
    train_idx <- 1:ncol(dat)
    combined_corrected <- adjust_rankin(combined_dat, train_indices = train_idx, debug = FALSE)
    dat_corrected <- combined_corrected[, 1:ncol(dat), drop = FALSE]
    dat_test_corrected <- combined_corrected[, (ncol(dat) + 1):ncol(combined_corrected), drop = FALSE]
    return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

  } else {
    stop(sprintf("Unknown batch correction method: %s", method))
  }
}

# ====================================================================
# GLOBAL SCALING
# ====================================================================

global_scale <- function(dat_train, dat_test) {
  train_mean <- mean(dat_train, na.rm = TRUE)
  train_sd <- sd(as.vector(dat_train), na.rm = TRUE)
  if (is.na(train_mean) || is.na(train_sd) || train_sd == 0 || !is.finite(train_sd) || train_sd < 1e-10) {
    if (!is.na(train_mean) && is.finite(train_mean)) {
      dat_train_scaled <- dat_train - train_mean
      dat_test_scaled <- dat_test - train_mean
      return(list(dat_train = dat_train_scaled, dat_test = dat_test_scaled))
    } else {
      stop("Scaling produced invalid values (NA, 0, or infinite standard deviation)")
    }
  }
  dat_train_scaled <- (dat_train - train_mean) / train_sd
  dat_test_scaled <- (dat_test - train_mean) / train_sd
  if (any(is.nan(dat_train_scaled) | is.infinite(dat_train_scaled), na.rm = TRUE) || 
      any(is.nan(dat_test_scaled) | is.infinite(dat_test_scaled), na.rm = TRUE)) {
    stop("Scaling produced non-finite values (NaN or Inf) in output data")
  }
  list(dat_train = dat_train_scaled, dat_test = dat_test_scaled)
}

# ====================================================================
# MAIN SYSTEM CALL
# ====================================================================

main <- function() {
  tryCatch({
    cat("Step 1: Loading TB data...\n")
    data_path <- "data/TB_real_data.RData"
    if (!file.exists(data_path)) {
      stop(sprintf("Data file not found: %s", data_path))
    }
    load(data_path)
    
    # Define study order
    all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")

    # Determine reference and target studies
    if (is.null(num_datasets)) {
      # Default behavior: use GSE37250_SA as reference, India as target (for backward compatibility)
      ref_studies <- c("GSE37250_SA")
      cat("Step 2: Using default configuration (GSE37250_SA -> India)...\n")
    } else {
      # General case: use first num_datasets studies as reference, excluding the test study
      non_test_studies <- all_studies[all_studies != test_study]
      if (num_datasets > length(non_test_studies)) {
        stop(sprintf("num_datasets=%d exceeds available non-test studies (%d)", num_datasets, length(non_test_studies)))
      }
      ref_studies <- non_test_studies[1:num_datasets]
      cat(sprintf("Step 2: Using first %d non-test studies as reference: %s\n", num_datasets, paste(ref_studies, collapse=", ")))
      cat(sprintf("        Using '%s' as target study (excluded from step-1 training)\n", test_study))
    }

    # Find common genes across all reference and target studies
    target_idx <- which(all_studies == test_study)
    if (length(target_idx) == 0) {
      stop(sprintf("Test study '%s' not found in all_studies", test_study))
    }

    common_genes <- rownames(dat_lst[[ref_studies[1]]])
    for (study in c(ref_studies, test_study)) {
      common_genes <- intersect(common_genes, rownames(dat_lst[[study]]))
    }
    cat(sprintf("  Common genes count: %d genes\n", length(common_genes)))

    # Combine reference data
    dat <- NULL
    batch <- NULL
    group <- NULL
    for (study in ref_studies) {
      study_data <- dat_lst[[study]][common_genes, , drop=FALSE]
      study_labels <- label_lst[[study]]
      if (is.null(dat)) {
        dat <- study_data
        batch <- rep(study, ncol(study_data))
        group <- study_labels
      } else {
        dat <- cbind(dat, study_data)
        batch <- c(batch, rep(study, ncol(study_data)))
        group <- c(group, study_labels)
      }
    }

    dat_test <- dat_lst[[test_study]][common_genes, , drop=FALSE]
    group_test <- label_lst[[test_study]]
    
    # Format supervised groups
    group_binary <- ifelse(group == "Control" | group == "0", 0, 1)
    group_test_binary <- ifelse(group_test == "Control" | group_test == "0", 0, 1)
    
    # DATA PREPROCESSING: LOG TRANSFORMATION FOR RAW INTENSITY DATA
    # Detect if data appears to be raw intensity values that need log transformation
    needs_log_transform_train <- max(dat, na.rm = TRUE) > 100 || mean(dat, na.rm = TRUE) > 20 || (max(dat, na.rm = TRUE) / median(dat, na.rm = TRUE)) > 50
    needs_log_transform_test <- max(dat_test, na.rm = TRUE) > 100 || mean(dat_test, na.rm = TRUE) > 20 || (max(dat_test, na.rm = TRUE) / median(dat_test, na.rm = TRUE)) > 50
    
    if (needs_log_transform_train) {
      cat(sprintf("[PREPROCESSING] Applying log transformation to TRAINING data (GSE37250_SA)\n"))
      train_min <- min(dat, na.rm = TRUE)
      if (train_min < 0) {
        dat <- dat - train_min
      }
      dat <- dat + 1
      dat <- log2(dat)
    }
    
    if (needs_log_transform_test) {
      cat(sprintf("[PREPROCESSING] Applying log transformation to TEST data (India)\n"))
      test_min <- min(dat_test, na.rm = TRUE)
      if (test_min < 0) {
        dat_test <- dat_test - test_min
      }
      dat_test <- dat_test + 1
      dat_test <- log2(dat_test)
    }
    
    cat("Step 3: Performing batch correction...\n")
    corrected <- apply_batch_correction(
      dat = dat,
      batch = batch,
      group = group_binary,
      dat_test = dat_test,
      method = adjuster,
      group_test = group_test_binary
    )
    
    cat("Step 4: Performing global scaling...\n")
    is_constant_train <- FALSE
    is_constant_test <- FALSE
    tryCatch({
      if (nrow(corrected$dat_corrected) > 0 && ncol(corrected$dat_corrected) > 0) {
        val_train <- corrected$dat_corrected[1,1]
        if (is.finite(val_train)) {
          is_constant_train <- all(corrected$dat_corrected == val_train, na.rm = TRUE)
        }
      }
      if (nrow(corrected$dat_test_corrected) > 0 && ncol(corrected$dat_test_corrected) > 0) {
        val_test <- corrected$dat_test_corrected[1,1]
        if (is.finite(val_test)) {
          is_constant_test <- all(corrected$dat_test_corrected == val_test, na.rm = TRUE)
        }
      }
    }, error = function(e) {})
    
    if (is_constant_train || is_constant_test) {
      cat("[WARNING] Batch correction produced constant data - using original data instead\n")
      corrected$dat_corrected <- dat
      corrected$dat_test_corrected <- dat_test
    }
    
    scaled <- global_scale(corrected$dat_corrected, corrected$dat_test_corrected)
    
    # Save the unadjusted reference dataset (scaled and unscaled) to the same output directory
    ref_out_unscaled <- file.path(dirname(output_data), "reference_unadjusted.csv")
    ref_out_scaled <- file.path(dirname(output_data), "reference_unadjusted_scaled.csv")
    if (!file.exists(ref_out_unscaled)) {
      write.csv(dat, file = ref_out_unscaled, row.names = TRUE)
      cat(sprintf("✓ Successfully saved unadjusted reference data to: %s\n", ref_out_unscaled))
    }
    if (!file.exists(ref_out_scaled)) {
      scaled_baseline <- global_scale(dat, dat_test)
      write.csv(scaled_baseline$dat_train, file = ref_out_scaled, row.names = TRUE)
      cat(sprintf("✓ Successfully saved scaled unadjusted reference data to: %s\n", ref_out_scaled))
    }
    
    # Save the sample metadata (labels/groups and study/batch names) to the same output directory
    target_meta_out <- file.path(dirname(output_data), "target_metadata.csv")
    ref_meta_out <- file.path(dirname(output_data), "reference_metadata.csv")
    if (!file.exists(target_meta_out)) {
      target_meta <- data.frame(
        Sample_ID = colnames(dat_test),
        Batch = "India",
        Disease = as.character(group_test),
        stringsAsFactors = FALSE
      )
      write.csv(target_meta, file = target_meta_out, row.names = FALSE)
      cat(sprintf("✓ Successfully saved target metadata to: %s\n", target_meta_out))
    }
    if (!file.exists(ref_meta_out)) {
      ref_meta <- data.frame(
        Sample_ID = colnames(dat),
        Batch = as.character(batch),
        Disease = as.character(group),
        stringsAsFactors = FALSE
      )
      write.csv(ref_meta, file = ref_meta_out, row.names = FALSE)
      cat(sprintf("✓ Successfully saved reference metadata to: %s\n", ref_meta_out))
    }
    
    cat("Step 5: Writing adjusted target CSV to output file...\n")
    write.csv(scaled$dat_test, file = output_data, row.names = TRUE)
    cat(sprintf("✓ Successfully saved adjusted target data to: %s\n", output_data))
    
    # Save adjusted reference dataset
    output_ref <- gsub("_target\\.csv$", "_reference.csv", output_data)
    cat("Step 6: Writing adjusted reference CSV to output file...\n")
    write.csv(scaled$dat_train, file = output_ref, row.names = TRUE)
    cat(sprintf("✓ Successfully saved adjusted reference data to: %s\n", output_ref))
    
    # Save combined reference and target dataset
    output_combined <- gsub("_target\\.csv$", "_combined.csv", output_data)
    cat("Step 7: Writing adjusted combined CSV to output file...\n")
    combined_scaled <- cbind(scaled$dat_train, scaled$dat_test)
    write.csv(combined_scaled, file = output_combined, row.names = TRUE)
    cat(sprintf("✓ Successfully saved adjusted combined data to: %s\n", output_combined))
    
    return(0)
    
  }, error = function(e) {
    cat(sprintf("\n[ERROR] Pipeline failed: %s\n", e$message), file = stderr())
    return(1)
  })
}

exit_code <- main()
quit(status = exit_code)
