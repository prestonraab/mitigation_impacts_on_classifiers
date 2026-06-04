# adjuster_plot_utils.R
# Common utility functions for adjuster plotting

format_adjuster_label <- function(adjuster_name) {
  if (adjuster_name == "within_study_cv") return("Within-Study CV")
  if (adjuster_name == "unadjusted") return("Log Only")
  if (adjuster_name == "naive") return("Naive")
  if (adjuster_name == "rank_samples") return("Rank Features")
  if (adjuster_name == "rank_twice") return("Rank Twice")
  if (adjuster_name == "npn") return("NPN")
  if (adjuster_name == "combat") return("ComBat")
  if (adjuster_name == "combat_sup") return("ComBat Sup.")
  if (adjuster_name == "combat_mean") return("ComBat Mean")
  if (adjuster_name == "mnn") return("MNN")
  if (adjuster_name == "fast_mnn") return("FastMNN")
  if (adjuster_name == "ruvg") return("RUVg")
  if (adjuster_name == "cublock") return("CuBlock")
  if (adjuster_name == "angel") return("Angel")
  if (adjuster_name == "tdm") return("TDM")
  if (adjuster_name == "rnabc") return("RNAbc")
  if (adjuster_name == "shambhala2") return("Shambhala2")
  if (adjuster_name == "coconut") return("COCONUT")
  if (adjuster_name == "rankin") return("Rank-In")
  if (adjuster_name == "recombat") return("ReComBat")
  if (adjuster_name == "yugene") return("YuGene")
  return(tools::toTitleCase(adjuster_name))
}

get_classifier_ordering <- function(mxe_data, classifier_name) {
  ordering <- mxe_data %>%
    filter(classifier_label == classifier_name & !is.na(classifier_label)) %>%
    group_by(adjuster) %>%
    summarise(mean_mcc = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_mcc))
  
  # Ensure within_study_cv is always first (it's the reference baseline)
  if ("within_study_cv" %in% ordering$adjuster) {
    within_study_row <- ordering[ordering$adjuster == "within_study_cv", ]
    other_rows <- ordering[ordering$adjuster != "within_study_cv", ]
    ordering <- rbind(within_study_row, other_rows)
  }
  
  return(ordering)
}

create_adjuster_labels <- function(adjuster_order) {
  sapply(adjuster_order$adjuster, format_adjuster_label)
}

format_study_label <- function(study_name) {
  labels <- c(
    "GSE37250_SA" = "South Africa",
    "GSE37250_M"  = "Malawi",
    "GSE39941_M"  = "Children",
    "Africa"      = "Adolescents",
    "India"       = "UK"
  )
  ifelse(study_name %in% names(labels), labels[study_name], study_name)
}
