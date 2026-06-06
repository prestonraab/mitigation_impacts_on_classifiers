#!/usr/bin/env Rscript
# make_figure2_regularization.R
# Three-panel figure showing how regularization ignores flipped genes.

suppressMessages({
  library(sva); library(matrixStats); library(ggplot2); library(cowplot)
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
# Weight = beta / var
rda_weight <- abs(post_combat_beta) / vp_adjusted

df <- data.frame(
  gene = rownames(dat),
  initial_beta = initial_beta,
  post_beta = post_combat_beta,
  var_adj = vp_adjusted,
  rda_weight = rda_weight,
  is_flipped = sign(initial_beta) != sign(post_combat_beta) & abs(initial_beta) > 0.5
)

# Identify a severely flipped gene for Panel C
# We want one with a large positive initial beta, but very negative post beta
flip_candidates <- df %>% 
  filter(initial_beta > 1, post_beta < -5) %>%
  arrange(post_beta)
prime_gene <- flip_candidates$gene[1]

# Panel A: Post-ComBat Beta vs Variance
pA <- ggplot(df, aes(x = post_beta, y = log2(var_adj))) +
  geom_point(aes(color = is_flipped), alpha = 0.5, size = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#d73027"), guide="none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  labs(x = "Post-ComBat Effect Size (β*)", 
       y = expression(Log[2]~"Adjusted Variance"),
       title = "A  Variance Explosion for Flipped Genes") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 10))

# Panel B: Post-ComBat Beta vs RDA Weight
pB <- ggplot(df, aes(x = post_beta, y = rda_weight)) +
  geom_point(aes(color = is_flipped), alpha = 0.5, size = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#d73027"), guide="none") +
  geom_smooth(method = "loess", color = "black", linewidth = 0.8, se = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  labs(x = "Post-ComBat Effect Size (β*)", 
       y = "RDA Feature Weight (|β*| / Var*)",
       title = "B  Classifier Ignores High-Variance Flips") +
  coord_cartesian(ylim = c(0, quantile(df$rda_weight, 0.99, na.rm=TRUE))) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 10))

# Panel C: Boxplot of the prime gene
expr_df <- data.frame(
  expr = combat_dat[prime_gene, ],
  batch = bat,
  class = factor(lab, levels=c(0,1), labels=c("Control", "Active TB"))
)

pC <- ggplot(expr_df, aes(x = batch, y = expr, fill = class)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  geom_point(position = position_jitterdodge(jitter.width=0.1), size=0.5, alpha=0.5, aes(color=class)) +
  scale_fill_manual(values=c("Control"="#4e79a7", "Active TB"="#e07b39")) +
  scale_color_manual(values=c("Control"="#4e79a7", "Active TB"="#e07b39"), guide="none") +
  labs(x = "Training Cohort", y = "Adjusted Expression", 
       title = sprintf("C  Why it's ignored: Gene '%s'", prime_gene),
       fill = "Status") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 10),
        legend.position = "bottom")

top_row <- plot_grid(pA, pB, ncol = 2, align = "h")
fig2 <- plot_grid(top_row, pC, ncol = 1, rel_heights = c(1, 1.1))

ggsave(out_png, fig2, width = 9, height = 7, dpi = 150)
ggsave(out_pdf, fig2, width = 9, height = 7)
cat("->", out_png, "\n")
