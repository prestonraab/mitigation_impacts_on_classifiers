#!/usr/bin/env Rscript
# Two-panel plot of adjuster performance averaged over classifiers.
# Left panel:  gap vs within-study CV  (adjusters can't match no-batch-effect scenario)
# Right panel: gap vs unadjusted       (adjusters do help compared to doing nothing)

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

parser <- ArgumentParser(description = "Plot relative performance averaged over classifiers")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV (adjusters_on_classifiers.csv)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")
parser$add_argument("--show-log-only", action = "store_true", default = FALSE,
                    help = "Include the lower panel (vs. Log Only)")
opt <- parser$parse_args()

source("scripts/adjuster_plot_utils.R")

data <- read.csv(opt$input, stringsAsFactors = FALSE)
data$test_study <- format_study_label(data$test_study)

if (!is.null(opt$adjusters)) {
  selected <- unique(c(trimws(strsplit(opt$adjusters, ",")[[1]]), "within_study_cv"))
  data <- data[data$adjuster %in% selected, ]
}

mxe <- data[data$metric == "mcc" & !is.na(data$n_datasets) & data$classifier != "logistic", ]
mxe$n_train <- mxe$n_datasets - 1

# --- Helper: compute gaps relative to a reference adjuster, averaged over classifiers ---
compute_panel_data <- function(mxe, ref_adjuster, exclude_adjusters = c()) {
  baselines <- mxe %>%
    filter(adjuster == ref_adjuster) %>%
    select(classifier, n_datasets, test_study, value) %>%
    rename(baseline = value)

  rel <- mxe %>%
    filter(!adjuster %in% c(ref_adjuster, exclude_adjusters)) %>%
    inner_join(baselines, by = c("classifier", "n_datasets", "test_study")) %>%
    mutate(gap = value - baseline)

  # Average over classifiers
  avg <- rel %>%
    group_by(adjuster, n_datasets, test_study) %>%
    summarise(mean_gap = mean(gap, na.rm = TRUE), .groups = "drop") %>%
    mutate(n_train = n_datasets - 1)

  # Order by overall mean gap (descending)
  adj_order <- avg %>%
    group_by(adjuster) %>%
    summarise(overall = mean(mean_gap, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(overall))

  adj_labels <- sapply(adj_order$adjuster, format_adjuster_label)
  avg$adjuster_label <- factor(avg$adjuster, levels = adj_order$adjuster, labels = adj_labels)

  # Order test studies by average performance (best first)
  study_order <- avg %>%
    group_by(test_study) %>%
    summarise(avg_gap = mean(mean_gap, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg_gap))
  avg$test_study <- factor(avg$test_study, levels = study_order$test_study)

  grand <- avg %>%
    group_by(adjuster_label) %>%
    summarise(grand_mean = mean(mean_gap, na.rm = TRUE), .groups = "drop")

  # Significance: one-sample t-test per adjuster
  sig <- data.frame()
  alt <- if (all(avg$mean_gap <= 0, na.rm = TRUE)) "less" else "two.sided"
  for (lbl in adj_labels) {
    vals <- avg$mean_gap[avg$adjuster_label == lbl]
    if (length(vals) > 1) {
      tt <- t.test(vals, mu = 0, alternative = alt)
      sig <- rbind(sig, data.frame(adjuster_label = lbl, p_value = tt$p.value,
                                    mean_diff = mean(vals), n_obs = length(vals),
                                    stringsAsFactors = FALSE))
    }
  }
  if (nrow(sig) > 0) {
    sig$p_adj <- p.adjust(sig$p_value, method = "bonferroni")
    sig$sig_label <- sapply(sig$p_adj, function(p) {
      if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    })
  }

  list(avg = avg, grand = grand, sig = sig)
}

# --- Helper: build one panel ---
make_panel <- function(panel, title, y_label, show_sig = TRUE) {
  avg   <- panel$avg
  grand <- panel$grand
  sig   <- panel$sig

  y_max <- max(avg$mean_gap, na.rm = TRUE)
  y_min <- min(avg$mean_gap, na.rm = TRUE)
  y_range <- max(y_max - y_min, 0.01)

  p <- ggplot(avg, aes(x = adjuster_label, y = mean_gap)) +
    geom_point(aes(color = test_study, shape = as.factor(n_train)), size = 3.5, alpha = 0.5) +
    geom_crossbar(data = grand,
                  aes(x = adjuster_label, y = grand_mean, ymin = grand_mean, ymax = grand_mean),
                  color = "gray30", linewidth = 0.5, width = 0.5, fatten = 2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6) +
    scale_x_discrete(drop = FALSE, limits = rev) +
    scale_y_reverse() +
    coord_flip(clip = "off") +
    scale_color_brewer(palette = "Set2", name = "Test Study") +
    scale_shape_manual(values = c(15, 16, 17, 18), name = "# Train\nStudies") +
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 2, title.theme = element_text(size = 11, face = "bold", margin = margin(t = 12)))
    ) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      axis.ticks = element_line(color = "grey70", linewidth = 0.3),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "grey40", linewidth = 0.4),
      axis.line.y = element_blank(),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.55, "cm"),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.spacing.y = unit(0.15, "cm"),
      legend.box.spacing = unit(0.4, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.position = "right",
      plot.margin = margin(5, 5, 5, 5),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 13, hjust = 0.5, face = "bold")
    ) +
    labs(y = y_label, title = title)

  if (show_sig && nrow(sig) > 0) {
    sig$y_pos <- y_max + 0.06 * y_range
    p <- p + geom_text(data = sig,
                       aes(x = adjuster_label, y = y_pos, label = sig_label),
                       inherit.aes = FALSE, size = 4.5, fontface = "bold", hjust = 0)
  }
  p
}

# --- Build both panels ---
panel_cv  <- compute_panel_data(mxe, ref_adjuster = "within_study_cv")
panel_unadj <- compute_panel_data(mxe, ref_adjuster = "unadjusted",
                                   exclude_adjusters = "within_study_cv")

p_left  <- make_panel(panel_cv,
                       title = if (opt$show_log_only) "Adjuster Performance vs. Within-Study CV" else "",
                       y_label = "\u2190 Better                    Difference in Mean MCC (Adjuster - Within-Study CV)          ",
                       show_sig = opt$show_log_only)

if (opt$show_log_only) {
  p_right <- make_panel(panel_unadj,
                         title = "Adjuster Performance vs. Log Only Transformation",
                         y_label = "\u2190 Better | Difference in Mean MCC")

  arrow_label <- ggdraw()

  combined <- plot_grid(
    p_left  + theme(legend.position = "none"),
    arrow_label,
    p_right + theme(legend.position = "none"),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 0.06, 1)
  )
  legend <- get_legend(p_left + theme(legend.position = "right"))
  final  <- plot_grid(combined, legend, rel_widths = c(1, 0.15))

  ggsave(filename = opt$output, plot = final,
         width = 10, height = 8, dpi = 300, units = "in", bg = "white")
} else {
  legend <- get_legend(p_left + theme(legend.position = "right"))
  final  <- plot_grid(
    p_left + theme(legend.position = "none"),
    legend, rel_widths = c(1, 0.15)
  )

  ggsave(filename = opt$output, plot = final,
         width = 10, height = 5, dpi = 300, units = "in", bg = "white")
}
cat("Saved:", opt$output, "\n")

# Save significance results from both panels
all_sig <- rbind(
  if (nrow(panel_cv$sig) > 0) cbind(panel_cv$sig, panel = "vs_within_study_cv") else NULL,
  if (nrow(panel_unadj$sig) > 0) cbind(panel_unadj$sig, panel = "vs_unadjusted") else NULL
)
if (!is.null(all_sig) && nrow(all_sig) > 0) {
  sig_out <- sub("\\.png$", "_significance.csv", opt$output)
  write.csv(all_sig, sig_out, row.names = FALSE)
  cat("Significance results saved to:", sig_out, "\n")
}
