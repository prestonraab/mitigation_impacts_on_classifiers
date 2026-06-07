#!/usr/bin/env Rscript
# make_figure1_mechanism.R
# Five-panel mechanism figure for the supervised ComBat x KNN chapter section.
#
# Panel A: Ridgeplot of per-gene delta^2 across training studies in the canonical failing
#          scenario (n=5, test=South Africa). Dashed line at delta^2=1 marks the flip
#          threshold; ~60% of genes fall below it.
#
# Panel B: Ridgeplot of the biology fraction -- the share of each gene's delta^2
#          attributable to class-effect heterogeneity (Delta_beta^2 * Var(X) / sigma_p^2)
#          / delta^2 -- per training study.
#
# Panels C1-C3: Per-study class signal projected onto the consensus beta direction.
#   C1: Observed (pre-ComBat) effect size per study.
#   C2: Per-study deviation from the pooled beta.
#   C3: Post-ComBat effect size after delta-reweighting.

# pixi run Rscript scripts/make_figure1_mechanism.R -i data/TB_real_data.RData -o outputs/diagnostics/hypothesis_tests/figure1_mechanism.png --utils /home/phr23/mitigation_impacts_on_classifiers/scripts/adjuster_plot_utils.R 2>&1

suppressMessages({
  library(sva); library(matrixStats); library(ggplot2); library(cowplot)
  library(dplyr); library(ggridges); library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("-i", "--input",  required = TRUE,  help = "TB_real_data.RData")
parser$add_argument("-o", "--output", required = TRUE,  help = "Output PNG path")
parser$add_argument("--utils", default = NULL, help = "Path to adjuster_plot_utils.R")
args <- parser$parse_args()

if (!is.null(args$utils) && file.exists(args$utils)) source(args$utils) else
  format_study_label <- function(x) x

load(args$input)
out_png <- args$output
out_pdf <- sub("\\.png$", ".pdf", out_png)

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin      <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
log_safe <- function(m) {
  if (max(m, na.rm=TRUE) > 100) {
    mn <- min(m, na.rm=TRUE); if (mn < 0) m <- m - mn; m <- log2(m + 1)
  }; m
}

# ── Canonical failing scenario ────────────────────────────────────────────────
n <- 5; test <- "GSE37250_SA"
refs <- ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]
cg   <- Reduce(intersect, lapply(c(refs, test), function(s) rownames(dat_lst[[s]])))
dat  <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg,])))
bat  <- as.factor(do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]])))))
lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
ubat <- levels(bat); nb <- length(ubat); N <- ncol(dat); G <- nrow(dat)
study_labels <- format_study_label(ubat)

# ── ComBat internals ──────────────────────────────────────────────────────────
design      <- cbind(model.matrix(~ -1 + bat), class = lab)
B           <- solve(crossprod(design), crossprod(design, t(dat)))
class_effect <- B["class",]                             # β_g
nbi         <- as.numeric(table(bat))
grand       <- as.numeric(crossprod(nbi / N, B[1:nb,]))
fitted      <- t(design %*% B)
vp          <- rowMeans((dat - fitted)^2)               # pooled variance
stand_mean  <- outer(grand, rep(1, N)) + outer(class_effect, lab)
s_data      <- (dat - stand_mean) / sqrt(pmax(vp, 1e-8))

# Per-batch quantities
gap_i  <- sapply(ubat, function(b) {
  i <- bat == b
  rowMeans(dat[, i & lab == 1, drop=FALSE]) - rowMeans(dat[, i & lab == 0, drop=FALSE])
})  # G × nb
n1 <- sapply(ubat, function(b) sum(bat == b & lab == 1))
n0 <- sapply(ubat, function(b) sum(bat == b & lab == 0))

# δ_raw per gene × batch
delta_raw <- sapply(ubat, function(b) rowVars(s_data[, bat == b, drop=FALSE]))  # G × nb

# Per-batch: within-class variance (vl) and mismatch term
vl_mat       <- matrix(NA_real_, G, nb)
mismatch_mat <- matrix(NA_real_, G, nb)
for (bi in seq_len(nb)) {
  i  <- bat == ubat[bi]
  v1 <- rowVars(dat[, i & lab == 1, drop=FALSE])
  v0 <- rowVars(dat[, i & lab == 0, drop=FALSE])
  vl_mat[, bi]       <- pmax((pmax(n1[bi]-1,1)*v1 + pmax(n0[bi]-1,1)*v0) /
                               pmax(n1[bi]+n0[bi]-2, 1), 1e-8)
  p_b                 <- n1[bi] / (n1[bi] + n0[bi])
  mismatch_mat[, bi] <- (gap_i[, bi] - class_effect)^2 * p_b * (1 - p_b)
}

# Normalise by pooled variance so both axes are dimensionless ratios
observed  <- as.vector(delta_raw)                                           # δ²
predicted <- as.vector(sweep(vl_mat + mismatch_mat, 1, pmax(vp, 1e-8), "/")) # (vl + mm) / vp
flip_vec  <- observed < 1
study_vec <- factor(rep(study_labels, each = G), levels = study_labels)

# ── Panel A: δ² ridgeplot, one ridge per training study ───────────────────────
clr_flip  <- "#e07b39"
clr_nflip <- "#4e79a7"

dr_long <- data.frame(
  d2    = as.vector(delta_raw),
  study = factor(rep(study_labels, each = G), levels = rev(study_labels))
)

pA <- ggplot(dr_long, aes(x = d2, y = study, fill = study, colour = study)) +
  geom_hline(yintercept = seq_len(nb), linewidth = 0.3, colour = "grey85") +
  geom_density_ridges(alpha = 0.70, scale = 0.90, rel_min_height = 0.01) +
  geom_vline(xintercept = 1, linewidth = 0.8, colour = "grey20", linetype = "dashed") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  scale_x_continuous(limits = c(0, 4.5), expand = c(0, 0)) +
  labs(x = expression(hat(delta)^2 ~ "(batch / pooled gene variance)"),
       y = NULL,
       title = "A  Genewise Distribution of\nEstimated Noise Inflation") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 10))

# ── Panel B: ridgeplot of biology fraction of δ² ─────────────────────────────
# Fraction = (Δβ²·Var(X)/σp²) / δ²: how much of each gene's variance scaling
# is class-gap heterogeneity rather than within-class noise.
bio_frac_mat <- sweep(mismatch_mat, 1, pmax(vp, 1e-8), "/") / pmax(delta_raw, 1e-8)
bio_frac_mat <- pmin(pmax(bio_frac_mat, 0), 1)

bf_long <- data.frame(
  frac  = as.vector(bio_frac_mat),
  study = factor(rep(study_labels, each = G), levels = rev(study_labels))
)

pB <- ggplot(bf_long, aes(x = frac, y = study, fill = study, colour = study)) +
  geom_hline(yintercept = seq_len(nb), linewidth = 0.3, colour = "grey85") +
  geom_density_ridges(alpha = 0.70, scale = 0.90, rel_min_height = 0.01) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = expression((Delta*beta^2 * Var(X) / sigma[p]^2) / hat(delta)^2),
       y = NULL,
       title = "B  Fraction of δ² Explained By\nEffect Size Differences") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 10))

# ── Panel C: three-stage view of the class signal ────────────────────────────
# All values are per-study scalars: projection of each study's class gap onto
# the consensus direction β / ‖β‖², so β itself maps to 1.0.
beta_norm2      <- sum(class_effect^2)
interaction_mat <- sweep(gap_i, 1, class_effect, "-")            # G × nb: gap_i − β

proj_raw_i        <- colSums(class_effect * gap_i)                                / beta_norm2
corrected_gap_mat <- class_effect + interaction_mat / sqrt(pmax(delta_raw, 1e-8)) # G × nb
proj_corrected_i  <- colSums(class_effect * corrected_gap_mat)                   / beta_norm2
proj_dev_i        <- proj_raw_i - 1                               # deviation from β

# Does the original Beta align with a new estimate of the class effect after ComBat correction?
post_combat_B <- solve(crossprod(design), crossprod(design, t(dat - fitted)))
post_combat_beta <- post_combat_B["class",]
# Cosine similarity between original and post-ComBat beta
cos_sim <- sum(class_effect * post_combat_beta) / (sqrt(sum(class_effect^2)) * sqrt(sum(post_combat_beta^2)))
cat("Cosine similarity between original and post-ComBat beta:", round(cos_sim, 4), "\n")

# Named color vector: same assignment as A/B (levels = rev(study_labels) → Set2 1:nb)
study_cols <- setNames(rev(RColorBrewer::brewer.pal(nb, "Set2")), study_labels)

make_c_panel <- function(vals, title, ref_y, ylab = "") {
  df <- data.frame(study = factor(study_labels, levels = study_labels), value = vals)
  ggplot(df, aes(x = study, y = value, fill = study)) +
    geom_col(colour = "grey30", linewidth = 0.25, width = 0.7) +
    scale_fill_manual(values = study_cols, guide = "none") +
    geom_hline(yintercept = ref_y, linewidth = 0.7, colour = "#1a9850", linetype = "dashed") +
    geom_hline(yintercept = 0,     linewidth = 0.4, colour = "grey40") +
    annotate("text", x = -Inf, y = ref_y, hjust = -0.1, vjust = -0.4,
             label = if (abs(ref_y - 1) < 0.01) "β" else "0",
             size = 3.2, colour = "#1a9850") +
    labs(x = NULL, y = ylab, title = title) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
          plot.title  = element_text(face = "bold", size = 10))
}

# Does the original Beta align with a new estimate of the class effect after ComBat correction?
# We must regress the fully adjusted ComBat data, not the residuals
mod <- model.matrix(~ lab)
combat_dat <- sva::ComBat(dat = dat, batch = bat, mod = mod)
post_combat_B <- solve(crossprod(design), crossprod(design, t(combat_dat)))
post_combat_beta <- post_combat_B["class",]

# Cosine similarity between original and post-ComBat beta
cos_sim <- sum(class_effect * post_combat_beta) / (sqrt(sum(class_effect^2)) * sqrt(sum(post_combat_beta^2)))
cat("\nCosine similarity between original and post-ComBat beta:", round(cos_sim, 4), "\n\n")



pC1 <- make_c_panel(proj_raw_i,      "C1  Observed Effect Size",   ref_y = 1,
                    ylab = "Projected onto β direction")
pC2 <- make_c_panel(proj_dev_i,      "C2  Deviation From β",       ref_y = 0)
pC3 <- make_c_panel(proj_corrected_i,"C3  Post-ComBat Effect Size", ref_y = 1)

bottom_row <- plot_grid(pC1, pC2, pC3, nrow = 1)

# Does the original Beta align with a new estimate of the class effect after ComBat correction?
# We must regress the fully adjusted ComBat data, not the residuals
mod <- model.matrix(~ lab)
combat_dat <- sva::ComBat(dat = dat, batch = bat, mod = mod)
post_combat_B <- solve(crossprod(design), crossprod(design, t(combat_dat)))
post_combat_beta <- post_combat_B["class",]

# Cosine similarity between original and post-ComBat beta
cos_sim <- sum(class_effect * post_combat_beta) / (sqrt(sum(class_effect^2)) * sqrt(sum(post_combat_beta^2)))
cat("\nCosine similarity between original and post-ComBat beta:", round(cos_sim, 4), "\n\n")

# ── Assemble ─────────────────────────────────────────────────────────────────
top_row <- plot_grid(pA, pB, ncol = 2, labels = NULL)
fig1    <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 0.9))

ggsave(out_png, fig1, width = 9, height = 7, dpi = 150)
ggsave(out_pdf, fig1, width = 9, height = 7)
cat("->", out_png, "\n")
