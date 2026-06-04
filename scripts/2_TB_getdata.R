# ==============================================================================
# TB Data Retrieval & Harmonization Script
# ==============================================================================
# Purpose: Integrates heterogeneous TB gene expression datasets for batch effect analysis.
#
# COHORT DEMOGRAPHICS & SOURCES:
# 1. Batch 1 (USA): GSE73408 - Walter et al. (2016). 
#    - Population: US Adults (Denver, CO). Mixed demographics.
#    - Platform: Affymetrix HuGene 1.1 ST.
#
# 2. Batch 2 (Africa): GSE79362 - Zak et al. (2016).
#    - Population: South African Adolescents (Western Cape). Predominantly Black African/Mixed.
#    - Platform: Illumina RNA-seq.
#
# 3. Batch 3 (India Proxy): GSE107994 - Leong et al. (2018).
#    - Population: UK Residents of South Asian descent (Leicester). 
#    - Note: Used as a proxy for India due to genetic ancestry.
#    - Platform: Illumina RNA-seq (HT-12 v4 equivalent).
#
# 4. Validation (SA/Malawi): GSE37250 & GSE39941 - Kaforou/Berry.
#    - Population: South African & Malawian adults.
#    - Platform: Illumina HumanHT-12 v4.
# ==============================================================================

rm(list=ls())

# --- Dependencies ---
required_packages <- c("GEOquery", "annotate", "hugene11sttranscriptcluster.db",
                       "SummarizedExperiment", "limma", "BatchQC", "ggplot2", 
                       "readxl", "matrixStats", "illuminaHumanv4.db")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Package missing:", pkg, "- Please install via BiocManager::install()"))
  }
}

set.seed(123)
options(warn = 1) # Print warnings as they occur

# ==============================================================================
# 1. Process Batch 1: USA (GSE73408 - Walter et al.)
# ==============================================================================
cat("\n--- Processing Batch 1: USA (GSE73408) ---\n")

# Download Series Matrix
gse <- getGEO("GSE73408", destdir="data", GSEMatrix=TRUE)[[1]]

# Annotate Probes to Gene Symbols
x <- hugene11sttranscriptclusterSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])

# Map and filter NA symbols
gene_symbols <- sapply(featureNames(gse), function(s) ifelse(s %in% names(xx), xx[[s]], NA))
gse <- gse[!is.na(gene_symbols), ]
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

dat_usa <- exprs(gse)
rownames(dat_usa) <- gene_symbols

# Metadata Cleaning
pheno <- pData(gse)
disease_status <- pheno$characteristics_ch1.2 # "clinical group: TB/LTBI/PNA"

# Filter: Keep Active TB (1) and LTBI (0), Remove Pneumonia
is_tb <- grepl("clinical group: TB$", disease_status, ignore.case=TRUE)
is_control <- grepl("clinical group: LTBI", disease_status, ignore.case=TRUE)
is_pna <- grepl("clinical group: PNA", disease_status, ignore.case=TRUE)

keep_usa <- !is_pna & (is_tb | is_control)
dat_usa <- dat_usa[, keep_usa]
group_usa <- as.numeric(is_tb[keep_usa]) # 1=TB, 0=Control

cat(sprintf("  > Retained %d samples (%d TB, %d Control)\n", 
            ncol(dat_usa), sum(group_usa==1), sum(group_usa==0)))

# ==============================================================================
# 2. Process Batch 2: Africa (GSE79362 - Zak et al.)
# ==============================================================================
cat("\n--- Processing Batch 2: Africa (GSE79362) ---\n")

if (!file.exists("data/combined_sub.RData")) {
  stop("Critical Error: 'data/combined_sub.RData' is missing. This file contains the processed RNA-seq data.")
}

load("data/combined_sub.RData")

# Extract Training set (Zak et al.)
dat_africa <- train_expr
group_africa <- as.numeric(y_train) - 1 # Convert factor levels 1/2 -> 0/1

# Ensure matrix format
if (is(dat_africa, "SummarizedExperiment")) dat_africa <- assay(dat_africa, 1)

cat(sprintf("  > Loaded %d samples (%d TB, %d Control)\n", 
            ncol(dat_africa), sum(group_africa==1), sum(group_africa==0)))

# ==============================================================================
# 3. Process Batch 3: India Proxy (GSE107994 - Leicester UK)
# ==============================================================================
cat("\n--- Processing Batch 3: India Proxy (Leicester - GSE107994) ---\n")

excel_path <- "data/GSE107994/GSE107994_edgeR_normalized_Leicester_with_progressor_longitudinal.xlsx"
meta_path <- "data/GSE107994_sample_info.csv"

if (!file.exists(excel_path) || !file.exists(meta_path)) {
  stop("Missing GSE107994 files. Ensure the Excel and CSV are in data/GSE107994/")
}

# Load Data
expr_data <- read_excel(excel_path, sheet = 1)
meta_data <- read.csv(meta_path, stringsAsFactors = FALSE)

# Formatting Matrix
# Assumption: Columns 1-3 are metadata (Genes, Gene_name, Gene_biotype)
gene_names <- expr_data$Gene_name
dat_india_raw <- as.matrix(expr_data[, -(1:3)])
rownames(dat_india_raw) <- gene_names

# Identify Samples (Active TB vs Control)
active_tb_ids <- meta_data$title[grepl("Active_TB", meta_data$group)]
control_ids <- meta_data$title[grepl("^group: Control$", meta_data$group)]
target_samples <- c(active_tb_ids, control_ids)

# Filter Matrix columns
valid_cols <- colnames(dat_india_raw) %in% target_samples
dat_india <- dat_india_raw[, valid_cols]

# Generate Labels
# Re-match column names to metadata to ensure order is correct
matched_meta <- meta_data[match(colnames(dat_india), meta_data$title), ]
group_india <- as.numeric(grepl("Active_TB", matched_meta$group))

cat(sprintf("  > Retained %d samples (%d TB, %d Control)\n", 
            ncol(dat_india), sum(group_india==1), sum(group_india==0)))

# ==============================================================================
# 4. Harmonization & Intersection
# ==============================================================================
cat("\n--- Harmonizing Gene Sets ---\n")

common_genes <- Reduce(intersect, list(rownames(dat_usa), rownames(dat_africa), rownames(dat_india)))
cat(sprintf("  > Overlapping genes across all 3 training batches: %d\n", length(common_genes)))

# Subset
dat_usa <- dat_usa[common_genes, ]
dat_africa <- dat_africa[common_genes, ]
dat_india <- dat_india[common_genes, ]

# Variance Filtering (Quality Control)
# Keep genes with variance > 0 and detected in > 2 samples per batch
pass_usa <- rowVars(dat_usa) > 0 & rowSums(dat_usa != 0) > 2
pass_afr <- rowVars(dat_africa) > 0 & rowSums(dat_africa != 0) > 2
pass_ind <- rowVars(dat_india) > 0 & rowSums(dat_india != 0) > 2

keep_final <- pass_usa & pass_afr & pass_ind
cat(sprintf("  > Genes passing variance filter in ALL batches: %d\n", sum(keep_final)))

# Fallback if intersection is too strict
if (sum(keep_final) < 5000) {
  cat("  ! Warning: Strict filter dropped too many genes. Relaxing to 'present in 2/3 batches'.\n")
  keep_final <- (pass_usa & pass_afr) | (pass_afr & pass_ind) | (pass_usa & pass_ind)
  cat(sprintf("  > Genes passing relaxed filter: %d\n", sum(keep_final)))
}

# Build Final List
dat_lst <- list(
  USA = dat_usa[keep_final, ], 
  Africa = dat_africa[keep_final, ], 
  India = dat_india[keep_final, ]
)
label_lst <- list(
  USA = group_usa, 
  Africa = group_africa, 
  India = group_india
)

# ==============================================================================
# 5. Process Validation Sets (GSE37250 / GSE39941)
# ==============================================================================
cat("\n--- Processing Validation Batches ---\n")

# Helper to annotate Illumina V4
annotate_illumina <- function(dat) {
  mapped <- mapIds(illuminaHumanv4.db, keys=rownames(dat), column="SYMBOL", keytype="PROBEID", multiVals="first")
  dat <- dat[!is.na(mapped), ]
  rownames(dat) <- mapped[!is.na(mapped)]
  return(dat)
}

# Download GSE37250
gse37 <- getGEO("GSE37250", destdir="data", GSEMatrix=TRUE)[[1]]
p37 <- pData(gse37)

# Filter SA & Malawi
idx_sa <- p37$`hiv status:ch1`=="HIV negative" & 
  p37$`geographical region:ch1`=="South Africa" & 
  p37$`disease state:ch1` %in% c("active tuberculosis", "latent TB infection")

idx_mw <- p37$`hiv status:ch1`=="HIV negative" & 
  p37$`geographical region:ch1`=="Malawi" & 
  p37$`disease state:ch1` %in% c("active tuberculosis", "latent TB infection")

dat_sa_val <- annotate_illumina(exprs(gse37[, idx_sa]))
dat_mw_val <- annotate_illumina(exprs(gse37[, idx_mw]))

grp_sa_val <- as.numeric(p37$`disease state:ch1`[idx_sa] == "active tuberculosis")
grp_mw_val <- as.numeric(p37$`disease state:ch1`[idx_mw] == "active tuberculosis")

# Download GSE39941
gse39 <- getGEO("GSE39941", destdir="data", GSEMatrix=TRUE)[[1]]
p39 <- pData(gse39)

idx_mw2 <- p39$`hiv status:ch1`=="HIV negative" & 
  p39$`geographical region:ch1`=="Malawi" & 
  p39$`disease status:ch1` %in% c("active tuberculosis", "latent TB infection")

dat_mw_val2 <- annotate_illumina(exprs(gse39[, idx_mw2]))
grp_mw_val2 <- as.numeric(p39$`disease status:ch1`[idx_mw2] == "active tuberculosis")

# ==============================================================================
# 6. Final Merge and Save
# ==============================================================================
cat("\n--- Finalizing Data Object ---\n")

# Align Validation sets to the training genes
train_genes <- rownames(dat_lst$USA)
dat_sa_val <- dat_sa_val[intersect(rownames(dat_sa_val), train_genes), ]
dat_mw_val <- dat_mw_val[intersect(rownames(dat_mw_val), train_genes), ]
dat_mw_val2 <- dat_mw_val2[intersect(rownames(dat_mw_val2), train_genes), ]

# Note: This effectively drops genes not found in validation sets. 
# A robust pipeline might impute or handle NAs, but for now we intersect.
final_genes <- Reduce(intersect, list(train_genes, rownames(dat_sa_val), rownames(dat_mw_val), rownames(dat_mw_val2)))

# Subset everything one last time
dat_lst <- lapply(dat_lst, function(x) x[final_genes, ])
dat_lst$GSE37250_SA <- dat_sa_val[final_genes, ]
dat_lst$GSE37250_M  <- dat_mw_val[final_genes, ]
dat_lst$GSE39941_M  <- dat_mw_val2[final_genes, ]

label_lst$GSE37250_SA <- grp_sa_val
label_lst$GSE37250_M  <- grp_mw_val
label_lst$GSE39941_M  <- grp_mw_val2

cat(sprintf("Final gene count: %d\n", length(final_genes)))
cat("Sample distribution:\n")
print(sapply(label_lst, table))

save(dat_lst, label_lst, file="data/TB_real_data.RData")
cat("\n✅ Success! Saved to 'data/TB_real_data.RData'\n")
