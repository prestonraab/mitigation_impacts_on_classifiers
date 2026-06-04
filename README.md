# Mitigation Impacts on Classifiers

Does correcting for batch effects help or hurt cross-cohort classification? This pipeline benchmarks 22 batch correction methods against 7 classifiers on tuberculosis gene expression data, where each dataset was collected in a different country and on different platforms.

The core experiment: train on 2–5 cohorts combined, correct for batch effects, then test on a held-out cohort. Repeat across all valid train/test splits and both corrected and uncorrected data. The output is a performance table (AUC, MCC, accuracy, etc.) for every adjuster × classifier × split combination, plus visualizations that reveal which correction methods are broadly beneficial, which are harmful, and which interact badly with specific classifiers.

A secondary analysis (`classify_class_imbalanced.R`) tests the same methods under deliberate class imbalance, examining whether unequal case/control ratios in training data change the relative ranking of correction methods.

## Datasets

Six TB cohorts, all downloaded automatically from NCBI GEO:

| Study | Cohort | Reference |
|---|---|---|
| GSE37250 (SA) | South Africa | Kaforou / Berry |
| GSE73408 | USA | Walter et al. 2016 |
| GSE107994 | India | Leong et al. 2018 |
| GSE37250 (M) | Malawi | Kaforou / Berry |
| GSE79362 | Africa | Zak et al. 2016 |
| GSE39941 | Malawi | Berry |

## Batch correction methods evaluated

unadjusted, naive, rank_samples, rank_twice, NPN, ComBat, ComBat (mean-only), ComBat (supervised), MNN, FastMNN, RUVg, RUVr, YuGene, CuBlock, Angel, TDM, RNABC, Shambhala2, COCONUT, Rank-In, reComBat, reComBat (supervised)

## Quickstart

```bash
# 1. Install environment (first time only)
#    The recombat package requires this flag due to a deprecated transitive dependency.
SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True pixi install
pixi run install-all

# 2. Download and process data from GEO (pre-processed inputs are in data/)
pixi run Rscript scripts/2_TB_getdata.R

# 3. Run the full pipeline
pixi run snakemake --cores 8
```

The output folder defaults to `results/`. To write elsewhere, copy `.env.example` to `.env` and set `OUTPUT_FOLDER`.

## Repository structure

```
scripts/
  adjusters.R                   # All batch correction functions + dispatch table
  classify_adjusters.R          # Main classification job (one adjuster × classifier × split)
  classify_class_imbalanced.R   # Class-imbalance variant
  helper.R                      # Classifier training/evaluation utilities
  aggregate_results.R           # Combine per-job CSVs into summary table
  generate_main_plot.R          # Primary results figure
  plot_*.R                      # Additional figures
  2_TB_getdata.R                # Data acquisition from GEO
Snakefile                       # Workflow definition
config.yaml                     # Adjusters, classifiers, dataset sizes, study names
data/README.md                  # How to acquire data
```

## Leakage-free design

Test data is never seen during batch correction fitting. For methods that require a combined matrix (MNN, reComBat, Rank-In), training indices are passed explicitly so SVD/model fitting uses only training samples. For methods that need a reference (ComBat projection, COCONUT), the corrected training set is used as the reference batch and test data is projected in without using test labels.

The `--num-datasets` argument to `classify_adjusters.R` specifies the number of **training** cohorts; the test cohort is always held out and selected separately.
