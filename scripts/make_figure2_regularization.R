#!/usr/bin/env Rscript
# make_figure2_regularization.R
# Figure showing how regularization ignores flipped/distorted genes.

suppressMessages({
  library(sva); library(matrixStats); library(ggplot2)
  library(dplyr); library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("-i", "--input",  default = "data/TB_real_data.RData", help = "TB_real_data.RData")
parser$add_argument("-o", "--output", default = "results/figures/figure2_regularization.png", help = "Output PNG path")
args <- parser$parse_args()

load(args$input)
out_png <- args$output
out_pdf <- sub("\\.png$", ".pdf", out_png)
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin      <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
log_safe <- function(m) {
  if (max(m, na.rm=TRUE) > 100) {
    mn <- min(m, na.rm=TRUE); if (mn < 0) m <- m - mn; m <- log2(m + 1)
  }; m
}

# Setup Canonical scenario
n <- 5; test <- "GSE37250_SA"
refs <- ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]
cg   <- Reduce(intersect, lapply(c(refs, test), function(s) rownames(dat_lst[[s]])))
dat  <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg,])))
bat  <- as.factor(do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]])))))
lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))

N <- ncol(dat); G <- nrow(dat)

# Calculate initial beta
design <- cbind(model.matrix(~ -1 + bat), class = lab)
B <- solve(crossprod(design), crossprod(design, t(dat)))
initial_beta <- B["class",]

# Calculate initial pooled variance and Cohen's d
initial_fitted <- t(design %*% B)
vp_initial <- rowMeans((dat - initial_fitted)^2)
cohens_d <- initial_beta / sqrt(pmax(vp_initial, 1e-8))

# Run supervised ComBat
mod <- model.matrix(~ lab)
combat_dat <- ComBat(dat = dat, batch = bat, mod = mod)

# Calculate post-ComBat beta
post_combat_B <- solve(crossprod(design), crossprod(design, t(combat_dat)))
post_combat_beta <- post_combat_B["class",]

# Calculate post-ComBat Variance (Pooled Variance)
fitted <- t(design %*% post_combat_B)
vp_adjusted <- rowMeans((combat_dat - fitted)^2)

# RDA uses diagonal inverse-variance. 
# Weight = abs(beta) / var
rda_weight <- abs(post_combat_beta) / vp_adjusted

df <- data.frame(
  gene = rownames(dat),
  initial_beta = initial_beta,
  post_beta = post_combat_beta,
  var_adj = vp_adjusted,
  rda_weight = rda_weight,
  cohens_d = cohens_d
)

df$new_cohens_d <- df$post_beta / sqrt(pmax(df$var_adj, 1e-8))

# Sort so high-weight points are plotted on top
df <- df %>% arrange(rda_weight)

pB <- ggplot(df, aes(x = cohens_d, y = new_cohens_d)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey70") +
  geom_point(aes(color = rda_weight), alpha = 0.7, size = 1.2) +
  scale_color_viridis_c(option = "magma", name = "RDA Feature Weight\n(|β*| / Var*)") +
  labs(x = "Original Effect Size (Cohen's d)", 
       y = "Post-ComBat Effect Size (Cohen's d)",
       title = "Classifier Ignores Flipped Genes") +
  coord_cartesian(xlim = quantile(df$cohens_d, c(0.005, 0.995), na.rm=TRUE),
                  ylim = quantile(df$new_cohens_d, c(0.005, 0.995), na.rm=TRUE)) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

ggsave(out_png, pB, width = 7, height = 5, dpi = 150)
ggsave(out_pdf, pB, width = 7, height = 5)
cat("->", out_png, "\n")
