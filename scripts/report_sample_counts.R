#!/usr/bin/env Rscript
# Report active/latent sample counts per study from TB_real_data.RData

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Export sample count report")
parser$add_argument("-i", "--input", type = "character", default = "data/TB_real_data.RData",
                    help = "Input RData file")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output CSV file")
opt <- parser$parse_args()

load(opt$input)

report <- do.call(rbind, lapply(names(label_lst), function(name) {
  labs <- label_lst[[name]]
  data.frame(
    study = name,
    n_samples = length(labs),
    n_active = sum(labs == 1),
    n_latent = sum(labs == 0),
    stringsAsFactors = FALSE
  )
}))

write.csv(report, opt$output, row.names = FALSE)
cat("Saved:", opt$output, "\n")
