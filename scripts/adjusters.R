# adjusters.R
# Shared batch-adjustment functions used by classify_adjusters.R and
# classify_class_imbalanced.R. Source this file; do not run it directly.

# ── Low-level helpers ─────────────────────────────────────────────────────────

rank_normalized <- function(matrix_, dim) {
  if (dim < 1 || dim > 2) stop("Invalid dimension. Must be 1 (rows) or 2 (columns).")
  ranked <- apply(matrix_, dim, rank, ties.method = "average")
  # apply() transposes when dim=1; undo that
  if (dim == 1 && is.matrix(ranked)) ranked <- t(ranked)
  ranked / max(ranked, na.rm = TRUE)
}

# Row-wise standardize mat to match target_mean / target_sd.
# Genes with effectively zero variance are set to target_mean.
standardize_rows <- function(mat, target_mean, target_sd) {
  src_var <- apply(mat, 1, var)
  src_sd  <- sqrt(pmax(src_var, 1e-10))
  no_var  <- src_var <= 1e-10
  result  <- sweep(
    sweep(mat, 1, rowMeans(mat), `-`) / src_sd * target_sd,
    1, target_mean, `+`
  )
  result[no_var, ] <- target_mean[no_var]
  result
}

# Log2-transform matrix only when the data looks like raw intensities.
# Returns mat unchanged if it already appears log-scaled.
maybe_log2_transform <- function(mat) {
  if (max(mat) > 100 || mean(mat) > 20 || (max(mat) / median(mat)) > 50) {
    mat_min <- min(mat)
    if (mat_min < 0) mat <- mat - mat_min
    log2(mat + 1)
  } else {
    mat
  }
}

# ── Rank-based methods ────────────────────────────────────────────────────────

adjust_ranked_with_batch_info <- function(matrix_, batch) {
  ranked  <- rank_normalized(matrix_, 1)
  ranked2 <- matrix(NA, nrow = nrow(ranked), ncol = ncol(ranked))
  for (b in unique(batch)) {
    idx  <- which(batch == b)
    data <- ranked[, idx, drop = FALSE]
    ranked2[, idx] <- if (ncol(data) > 1) rank_normalized(data, 2) else data
  }
  if (any(is.na(ranked2))) ranked2[is.na(ranked2)] <- ranked[is.na(ranked2)]
  max_val <- max(ranked2, na.rm = TRUE)
  ranked2 / if (max_val == 0) 1 else max_val
}

adjust_ranked_samples_with_batch_info <- function(matrix_, batch) {
  result <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  for (b in unique(batch)) {
    idx  <- which(batch == b)
    data <- matrix_[, idx, drop = FALSE]
    result[, idx] <- if (ncol(data) > 1) {
      rank_normalized(data, 2)
    } else {
      data / max(data, na.rm = TRUE)
    }
  }
  if (any(is.na(result))) result[is.na(result)] <- matrix_[is.na(result)]
  max_val <- max(result, na.rm = TRUE)
  result / if (max_val == 0) 1 else max_val
}

adjust_ranked_twice_with_batch_info <- function(matrix_, batch) {
  result <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  for (b in unique(batch)) {
    idx  <- which(batch == b)
    data <- matrix_[, idx, drop = FALSE]
    result[, idx] <- if (ncol(data) > 1) {
      rank_normalized(rank_normalized(data, 1), 2)
    } else {
      rank_normalized(data, 1)
    }
  }
  if (any(is.na(result))) result[is.na(result)] <- matrix_[is.na(result)]
  max_val <- max(result, na.rm = TRUE)
  result / if (max_val == 0) 1 else max_val
}

# ── Other normalisation methods ───────────────────────────────────────────────

adjust_npn <- function(matrix_, batch) {
  if (is.null(batch)) {
    return(t(huge::huge.npn(t(matrix_), verbose = FALSE)))
  }
  result <- matrix_
  for (b in unique(batch)) {
    idx <- which(batch == as.character(b))
    if (length(idx) == 0) idx <- which(batch == b)
    result[, idx] <- t(huge::huge.npn(t(matrix_[, idx, drop = FALSE]), verbose = FALSE))
  }
  result
}

adjust_yugene <- function(matrix_) {
  result <- YuGene::YuGene(matrix_)
  as.matrix(unclass(result))
}

adjust_cublock <- function(matrix_, batch) {
  input_f  <- tempfile(fileext = ".csv")
  output_f <- tempfile(fileext = ".csv")
  write.table(matrix_, input_f, sep = ",", row.names = FALSE, col.names = FALSE)
  cmd <- sprintf(
    "octave --eval 'addpath(\".\"); pkg load statistics; data = csvread(\"%s\"); [norm_data] = CuBlock(data); csvwrite(\"%s\", norm_data);'",
    input_f, output_f
  )
  if (system(cmd) != 0) stop("Octave/CuBlock call failed. Verify Octave is installed and CuBlock.m is in the working directory.")
  result <- as.matrix(read.csv(output_f, header = FALSE))
  rownames(result) <- rownames(matrix_)
  colnames(result) <- colnames(matrix_)
  unlink(c(input_f, output_f))
  result
}

adjust_angel <- function(matrix_) {
  result <- apply(matrix_, 2, function(x) {
    r <- rank(x, ties.method = "average")
    (r - min(r)) / (max(r) - min(r))
  })
  if (!is.matrix(result)) result <- as.matrix(result)
  rownames(result) <- rownames(matrix_)
  colnames(result) <- colnames(matrix_)
  result
}

adjust_tdm <- function(dat_train, dat_test) {
  if (!requireNamespace("TDM",        quietly = TRUE)) stop("Package 'TDM' is required.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required for TDM.")
  train_dt <- data.table::as.data.table(dat_train, keep.rownames = "gene")
  test_dt  <- data.table::as.data.table(dat_test,  keep.rownames = "gene")
  res_dt   <- TDM::tdm_transform(target_data = test_dt, ref_data = train_dt)
  result   <- as.matrix(res_dt[, -1, with = FALSE])
  rownames(result) <- res_dt[[1]]
  result
}

adjust_rnabc <- function(dat_train, dat_test) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  target_dist     <- rowMeans(dat_train)
  combined        <- cbind(dat_test, target_dist)
  qnorm_res       <- limma::normalizeQuantiles(as.matrix(combined))
  dat_test_qnorm  <- qnorm_res[, 1:ncol(dat_test), drop = FALSE]
  rownames(dat_test_qnorm) <- rownames(dat_test)
  colnames(dat_test_qnorm) <- colnames(dat_test)
  combined_dat  <- cbind(dat_train, dat_test_qnorm)
  batch_vec     <- c(rep(1, ncol(dat_train)), rep(2, ncol(dat_test_qnorm)))
  combined_bc   <- sva::ComBat(combined_dat, batch = batch_vec, ref.batch = 1)
  combined_bc[, (ncol(dat_train) + 1):ncol(combined_bc), drop = FALSE]
}

adjust_shambhala2 <- function(matrix_, calib_P, ref_Q) {
  if (!exists("shambhala2_harmonize")) {
    if      (file.exists("Shambhala2/Shambhala2.R")) source("Shambhala2/Shambhala2.R")
    else if (file.exists("Shambhala2.R"))             source("Shambhala2.R")
    else stop("Shambhala-2 functions not found. Source Shambhala2.R from the cloned repository.")
  }
  shambhala2_harmonize(matrix_, P_dataset = calib_P, Q_dataset = ref_Q)
}

# ── Supervised methods ────────────────────────────────────────────────────────

adjust_coconut <- function(matrix_, batch, group) {
  if (!requireNamespace("COCONUT", quietly = TRUE)) stop("Package 'COCONUT' is required.")
  gse_list <- list()
  for (b in unique(batch)) {
    idx <- which(batch == b)
    pheno_df <- data.frame(
      disease_state = as.numeric(group[idx]),
      dummy = 1,
      row.names = colnames(matrix_[, idx])
    )
    gse_list[[as.character(b)]] <- list(pheno = pheno_df, genes = matrix_[, idx, drop = FALSE])
  }
  res    <- COCONUT::COCONUT(GSEs = gse_list, control.0.col = "disease_state")
  result <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_),
                   dimnames = dimnames(matrix_))
  for (b in names(res$COCONUTList)) {
    cols <- colnames(res$COCONUTList[[b]]$genes)
    result[, cols] <- as.matrix(res$COCONUTList[[b]]$genes)
  }
  for (b in names(res$controlList$GSEs)) {
    cols <- colnames(res$controlList$GSEs[[b]]$genes)
    result[, cols] <- as.matrix(res$controlList$GSEs[[b]]$genes)
  }
  result
}

adjust_rankin <- function(matrix_, n_svd = 1L, train_indices = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Package 'reticulate' is required for Rank-In.")
  rankin_script <- local({
    candidates <- c("scripts/rankin.py", "rankin.py")
    found <- Filter(file.exists, candidates)
    if (length(found) == 0) stop(sprintf("Rank-In script not found. Searched: %s", paste(candidates, collapse = ", ")))
    found[[1]]
  })
  reticulate::source_python(rankin_script)
  result_list <- rank_in_from_r(
    unname(as.list(as.data.frame(t(matrix_)))),
    train_indices = train_indices,
    n_svd = as.integer(n_svd)
  )
  result <- matrix(unlist(result_list), nrow = nrow(matrix_), ncol = ncol(matrix_), byrow = TRUE)
  rownames(result) <- rownames(matrix_)
  colnames(result) <- colnames(matrix_)
  result
}

adjust_recombat <- function(matrix_, batch) {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Package 'reticulate' is required.")
  pd           <- reticulate::import("pandas",   convert = FALSE)
  recombat_pkg <- reticulate::import("reComBat", convert = FALSE)
  data_pd      <- pd$DataFrame(t(matrix_))
  batch_pd     <- pd$Series(batch)
  model        <- recombat_pkg$reComBat(parametric = TRUE, model = "elastic_net")
  res_pd       <- model$fit_transform(data_pd, batch_pd)
  result       <- t(as.matrix(reticulate::py_to_r(res_pd)))
  rownames(result) <- rownames(matrix_)
  colnames(result) <- colnames(matrix_)
  result
}

adjust_recombat_sup <- function(matrix_, batch, group) {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Package 'reticulate' is required.")
  pd           <- reticulate::import("pandas",   convert = FALSE)
  recombat_pkg <- reticulate::import("reComBat", convert = FALSE)
  data_pd      <- pd$DataFrame(t(matrix_))
  batch_pd     <- pd$Series(batch)
  X_pd         <- pd$DataFrame(list(disease = as.integer(group)))
  model        <- recombat_pkg$reComBat(parametric = TRUE, model = "elastic_net")
  res_pd       <- model$fit_transform(data_pd, batch_pd, X = X_pd)
  result       <- t(as.matrix(reticulate::py_to_r(res_pd)))
  rownames(result) <- rownames(matrix_)
  colnames(result) <- colnames(matrix_)
  result
}

# ── Dispatch table ────────────────────────────────────────────────────────────
# Each entry: function(dat, dat_test, batch, group, ...)
# Returns:    list(dat_corrected, dat_test_corrected)

BATCH_CORRECTION_METHODS <- list(

  unadjusted = function(dat, dat_test, ...) {
    list(dat_corrected = dat, dat_test_corrected = dat_test)
  },

  naive = function(dat, dat_test, batch, ...) {
    overall_mean <- rowMeans(dat)
    overall_sd   <- sqrt(apply(dat, 1, var))
    dat_corrected <- dat
    for (b in unique(batch)) {
      idx <- which(batch == b)
      dat_corrected[, idx] <- standardize_rows(dat[, idx, drop = FALSE], overall_mean, overall_sd)
    }
    dat_test_corrected <- standardize_rows(dat_test, overall_mean, overall_sd)
    list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected)
  },

  rank_samples = function(dat, dat_test, batch, ...) {
    list(
      dat_corrected      = adjust_ranked_samples_with_batch_info(dat, batch),
      dat_test_corrected = adjust_ranked_samples_with_batch_info(dat_test, rep(1L, ncol(dat_test)))
    )
  },

  rank_twice = function(dat, dat_test, batch, ...) {
    list(
      dat_corrected      = adjust_ranked_twice_with_batch_info(dat, batch),
      dat_test_corrected = adjust_ranked_twice_with_batch_info(dat_test, rep(1L, ncol(dat_test)))
    )
  },

  npn = function(dat, dat_test, batch, ...) {
    list(
      dat_corrected      = adjust_npn(dat, batch),
      dat_test_corrected = adjust_npn(dat_test, NULL)
    )
  },

  combat = function(dat, dat_test, batch, ...) {
    dat_corrected <- ComBat(dat, batch = batch, mod = NULL)
    combined      <- cbind(dat_corrected, dat_test)
    combined_bc   <- ComBat(combined,
                            batch     = c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test))),
                            mod       = NULL,
                            ref.batch = 1L)
    list(dat_corrected      = dat_corrected,
         dat_test_corrected = combined_bc[, (ncol(dat_corrected) + 1):ncol(combined_bc)])
  },

  combat_mean = function(dat, dat_test, batch, ...) {
    dat_corrected <- ComBat(dat, batch = batch, mod = NULL, mean.only = TRUE)
    combined      <- cbind(dat_corrected, dat_test)
    combined_bc   <- ComBat(combined,
                            batch     = c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test))),
                            mod       = NULL, ref.batch = 1L, mean.only = TRUE)
    list(dat_corrected      = dat_corrected,
         dat_test_corrected = combined_bc[, (ncol(dat_corrected) + 1):ncol(combined_bc)])
  },

  combat_sup = function(dat, dat_test, batch, group, ...) {
    dat_corrected <- ComBat(dat, batch = batch, mod = model.matrix(~group))
    combined      <- cbind(dat_corrected, dat_test)
    combined_bc   <- ComBat(combined,
                            batch     = c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test))),
                            mod       = NULL,
                            ref.batch = 1L)
    list(dat_corrected      = dat_corrected,
         dat_test_corrected = combined_bc[, (ncol(dat_corrected) + 1):ncol(combined_bc)])
  },

  mnn = function(dat, dat_test, batch, ...) {
    suppressMessages({ library(batchelor, quietly = TRUE); library(SummarizedExperiment, quietly = TRUE) })
    test_id      <- max(batch) + 1L
    combined_bat <- c(batch, rep(test_id, ncol(dat_test)))
    u_batches    <- sort(unique(batch))
    merge_ord    <- c(u_batches[order(table(batch)[as.character(u_batches)], decreasing = TRUE)], test_id)
    mnn_out      <- tryCatch(
      batchelor::mnnCorrect(cbind(dat, dat_test), batch = combined_bat, merge.order = merge_ord),
      error = function(e) stop(sprintf("MNN correction failed: %s", e$message))
    )
    corrected <- SummarizedExperiment::assay(mnn_out)
    list(dat_corrected      = corrected[, 1:ncol(dat)],
         dat_test_corrected = corrected[, (ncol(dat) + 1):ncol(corrected)])
  },

  fast_mnn = function(dat, dat_test, batch, ...) {
    suppressMessages({ library(batchelor, quietly = TRUE); library(SummarizedExperiment, quietly = TRUE) })
    batch_mats <- lapply(sort(unique(batch)), function(b) dat[, which(batch == b), drop = FALSE])
    batch_mats[[length(batch_mats) + 1]] <- dat_test
    result <- tryCatch(
      do.call(batchelor::fastMNN, batch_mats),
      error = function(e) stop(sprintf("FastMNN failed: %s", e$message))
    )
    corrected <- as.matrix(SummarizedExperiment::assay(result, "reconstructed"))
    list(dat_corrected      = corrected[, 1:ncol(dat)],
         dat_test_corrected = corrected[, (ncol(dat) + 1):ncol(corrected)])
  },

  ruvr = function(dat, dat_test, batch, group, ...) {
    if (!requireNamespace("RUVSeq", quietly = TRUE)) stop("Package 'RUVSeq' is required.")
    residuals <- do.call(rbind, lapply(1:nrow(dat), function(i) {
      fit <- lm(dat[i, ] ~ group)
      residuals(fit)
    }))
    ruv_res      <- RUVSeq::RUVr(as.matrix(dat), 1:nrow(dat), k = 3, residuals = residuals, isLog = TRUE)
    dat_corrected <- ruv_res$normalizedCounts
    W            <- ruv_res$W
    gene_means   <- rowMeans(dat)
    alpha        <- solve(t(W) %*% W) %*% t(W) %*% t(dat - gene_means)
    W_test       <- t(dat_test - gene_means) %*% t(alpha) %*% solve(alpha %*% t(alpha))
    dat_test_corrected <- (dat_test - gene_means) - t(W_test %*% alpha) + gene_means
    list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected)
  },

  ruvg = function(dat, dat_test, batch, ...) {
    if (!requireNamespace("RUVSeq", quietly = TRUE)) stop("Package 'RUVSeq' is required.")
    hk_genes <- intersect(c("PPIA", "YWHAZ", "HPRT1", "RPLP0"), rownames(dat))
    if (length(hk_genes) < 2) {
      warning("Fewer than 2 housekeeping genes found; falling back to unadjusted.")
      return(list(dat_corrected = dat, dat_test_corrected = dat_test))
    }
    k        <- min(3, length(hk_genes) - 1)
    ruv_res  <- RUVSeq::RUVg(as.matrix(dat), hk_genes, k = k, isLog = TRUE)
    dat_corrected <- ruv_res$normalizedCounts
    W          <- ruv_res$W
    gene_means <- rowMeans(dat)
    alpha      <- solve(t(W) %*% W) %*% t(W) %*% t(dat - gene_means)
    alpha_hk   <- alpha[, which(rownames(dat) %in% hk_genes), drop = FALSE]
    W_test     <- t(dat_test[hk_genes, ] - gene_means[hk_genes]) %*% t(alpha_hk) %*% solve(alpha_hk %*% t(alpha_hk))
    dat_test_corrected <- (dat_test - gene_means) - t(W_test %*% alpha) + gene_means
    list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected)
  },

  yugene = function(dat, dat_test, ...) {
    list(
      dat_corrected      = adjust_yugene(dat),
      dat_test_corrected = adjust_yugene(dat_test)
    )
  },

  cublock = function(dat, dat_test, batch, ...) {
    list(
      dat_corrected      = adjust_cublock(dat, batch),
      dat_test_corrected = adjust_cublock(dat_test, rep(1L, ncol(dat_test)))
    )
  },

  angel = function(dat, dat_test, ...) {
    list(
      dat_corrected      = adjust_angel(dat),
      dat_test_corrected = adjust_angel(dat_test)
    )
  },

  tdm = function(dat, dat_test, ...) {
    list(
      dat_corrected      = dat,
      dat_test_corrected = adjust_tdm(dat_train = dat, dat_test = dat_test)
    )
  },

  rnabc = function(dat, dat_test, ...) {
    list(
      dat_corrected      = dat,
      dat_test_corrected = adjust_rnabc(dat_train = dat, dat_test = dat_test)
    )
  },

  shambhala2 = function(dat, dat_test, batch, ...) {
    batch_counts  <- table(batch)
    ref_idx       <- as.numeric(names(batch_counts)[which.max(batch_counts)])
    P_ref         <- dat[, which(batch == ref_idx), drop = FALSE]
    dat_corrected <- dat
    for (b in unique(batch)) {
      idx <- which(batch == b)
      dat_corrected[, idx] <- adjust_shambhala2(dat[, idx, drop = FALSE], P_ref, P_ref)
    }
    list(dat_corrected      = dat_corrected,
         dat_test_corrected = adjust_shambhala2(dat_test, P_ref, P_ref))
  },

  recombat = function(dat, dat_test, batch, ...) {
    combined   <- cbind(dat, dat_test)
    comb_batch <- c(batch, rep(max(batch) + 1L, ncol(dat_test)))
    corrected  <- adjust_recombat(combined, comb_batch)
    list(dat_corrected      = corrected[, 1:ncol(dat), drop = FALSE],
         dat_test_corrected = corrected[, (ncol(dat) + 1):ncol(corrected), drop = FALSE])
  },

  recombat_sup = function(dat, dat_test, batch, group, ...) {
    dat_corrected  <- adjust_recombat_sup(dat, batch, group)
    combined       <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test)))
    combined_bc    <- ComBat(combined, batch = combined_batch, mod = NULL, ref.batch = 1L)
    list(dat_corrected      = dat_corrected,
         dat_test_corrected = combined_bc[, (ncol(dat_corrected) + 1):ncol(combined_bc), drop = FALSE])
  },

  coconut = function(dat, dat_test, batch, group, ...) {
    dat_corrected <- adjust_coconut(dat, batch, group)
    combined      <- cbind(dat_corrected, dat_test)
    combined_bc   <- ComBat(combined,
                            batch     = c(rep(1L, ncol(dat_corrected)), rep(2L, ncol(dat_test))),
                            mod       = NULL,
                            ref.batch = 1L)
    list(dat_corrected      = dat_corrected,
         dat_test_corrected = combined_bc[, (ncol(dat_corrected) + 1):ncol(combined_bc)])
  },

  rankin = function(dat, dat_test, ...) {
    combined  <- cbind(dat, dat_test)
    corrected <- adjust_rankin(combined, train_indices = 1:ncol(dat))
    list(dat_corrected      = corrected[, 1:ncol(dat), drop = FALSE],
         dat_test_corrected = corrected[, (ncol(dat) + 1):ncol(corrected), drop = FALSE])
  }
)
