# generate_main_plot.R
# Generate the main adjuster effectiveness plot with classifier-specific ordering

generate_main_plot <- function(mxe_data, output_file, width = 20, height = 16, dpi = 300, share_y_axis = TRUE) {
  
  cat("Creating main figure\n")
  
  # Ensure dataset_label exists
  if (!"dataset_label" %in% colnames(mxe_data)) {
    mxe_data$dataset_label <- factor(paste(mxe_data$n_datasets, "studies"),
                                    levels = c("3 studies", "4 studies", "5 studies", "6 studies"))
  }
  
  # Add adjuster type
  mxe_data$adjuster_type <- "Batch Correction"
  mxe_data$adjuster_type[mxe_data$adjuster == "unadjusted"] <- "Original Data"
  mxe_data$adjuster_type[mxe_data$adjuster == "naive"] <- "Naive Correction"
  mxe_data$adjuster_type[mxe_data$adjuster == "rank_samples"] <- "Rank Adjustment"
  mxe_data$adjuster_type[mxe_data$adjuster == "rank_twice"] <- "Rank Adjustment"
  mxe_data$adjuster_type[mxe_data$adjuster == "npn"] <- "Quantile Normalization"
  mxe_data$adjuster_type[mxe_data$adjuster == "fast_mnn"] <- "FastMNN"
  
  # Calculate summary statistics
  sumstats <- mxe_data %>%
    group_by(adjuster_label, classifier_label, dataset_label, adjuster_type) %>%
    summarise(
      Avg = mean(value),
      Up = quantile(value, 0.975),
      Down = quantile(value, 0.025),
      .groups = "drop"
    )
  
  # Calculate frequency of best method for annotations
  freq_data <- mxe_data %>%
    group_by(classifier_label, dataset_label, test_study) %>%
    summarise(
      best_adjuster = adjuster_label[which.max(value)],
      .groups = "drop"
    ) %>%
    group_by(classifier_label, dataset_label, best_adjuster) %>%
    summarise(
      freq = n(),
      .groups = "drop"
    ) %>%
    group_by(classifier_label, dataset_label) %>%
    mutate(
      total = sum(freq),
      pct = freq / total
    ) %>%
    ungroup()
  
  # Add frequency annotations to summary stats
  sumstats <- sumstats %>%
    left_join(
      freq_data %>% 
        select(classifier_label, dataset_label, best_adjuster, pct) %>%
        rename(adjuster_label = best_adjuster),
      by = c("classifier_label", "dataset_label", "adjuster_label")
    ) %>%
    mutate(
      annot = ifelse(is.na(pct), "", scales::percent(pct, accuracy = 1))
    )
  
  # Determine which classifiers have data
  classifiers_with_data <- sumstats %>%
    filter(!is.na(classifier_label)) %>%
    group_by(classifier_label) %>%
    summarise(has_data = n() > 0, .groups = "drop") %>%
    filter(has_data) %>%
    pull(classifier_label)
  
  cat("Classifiers with data:", paste(classifiers_with_data, collapse = ", "), "\n")
  
  # Calculate global y-axis limits if sharing
  if (share_y_axis) {
    global_y_min <- 0
    global_y_max <- max(mxe_data$value, na.rm = TRUE)
    global_y_range <- global_y_max - global_y_min
    global_y_limits <- c(global_y_min - 0.05 * global_y_range, 
                         global_y_max + 0.15 * global_y_range)
    cat("Using shared y-axis limits:", round(global_y_limits[1], 3), "to", round(global_y_limits[2], 3), "\n")
  }
  
  # Create individual plots for each classifier
  plot_list <- list()
  
  for (classifier in classifiers_with_data) {
    raw_data <- mxe_data[mxe_data$classifier_label == classifier & !is.na(mxe_data$classifier_label), ]
    
    if (nrow(raw_data) == 0) {
      cat("Skipping", classifier, "- no data\n")
      next
    }
    
    # Get classifier-specific ordering
    classifier_adjuster_order <- get_classifier_ordering(mxe_data, classifier)
    classifier_specific_labels <- create_adjuster_labels(classifier_adjuster_order)
    
    cat("Classifier:", classifier, "\n")
    for (i in 1:nrow(classifier_adjuster_order)) {
      cat(sprintf("  %d. %s (mean MCC: %.4f)\n", 
                  i, classifier_specific_labels[i], classifier_adjuster_order$mean_mcc[i]))
    }
    
    # Apply ordering
    raw_data$adjuster_label <- factor(raw_data$adjuster,
                                      levels = classifier_adjuster_order$adjuster,
                                      labels = classifier_specific_labels)
    raw_data$dataset_label <- factor(raw_data$dataset_label, 
                                     levels = c("3 studies", "4 studies", "5 studies", "6 studies"))
    
    # Get plot_data with same ordering
    plot_data <- sumstats[sumstats$classifier_label == classifier & !is.na(sumstats$classifier_label), ]
    plot_data$adjuster_label <- factor(plot_data$adjuster_label, levels = classifier_specific_labels)
    plot_data$dataset_label <- factor(plot_data$dataset_label, 
                                     levels = c("3 studies", "4 studies", "5 studies", "6 studies"))
    
    # Calculate annotation position
    if (share_y_axis) {
      annotation_y <- global_y_limits[2] * 0.98
    } else {
      local_y_max <- max(raw_data$value, na.rm = TRUE)
      local_y_min <- min(raw_data$value, na.rm = TRUE)
      local_y_range <- local_y_max - local_y_min
      annotation_y <- local_y_max + 0.12 * local_y_range
    }
    
    # Calculate means
    mean_data <- raw_data %>%
      group_by(adjuster_label, dataset_label, adjuster_type) %>%
      summarise(mean_value = mean(value), .groups = "drop")
    
    # Create plot
    p <- ggplot(raw_data, aes(x = adjuster_label, y = value)) +
      geom_segment(data = mean_data, 
                  aes(x = as.numeric(adjuster_label) - 0.3, 
                      xend = as.numeric(adjuster_label) + 0.3,
                      y = mean_value, yend = mean_value),
                  color = "gray40", linewidth = 0.8) +
      geom_point(aes(color = test_study, shape = test_study), size = 3, alpha = 0.5) +
      geom_text(data = plot_data, aes(x = adjuster_label, y = annotation_y, label = annot), 
                color = "black", size = 2.0, vjust = 1) +
      facet_wrap(~ dataset_label, scales = "fixed", ncol = 4) +
      scale_x_discrete() +
      {if (share_y_axis) {
        scale_y_continuous(limits = global_y_limits, expand = expansion(mult = c(0, 0)))
      } else {
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
      }} +
      scale_color_brewer(palette = "Set2", name = "Test Study") +
      scale_shape_manual(values = c(15, 16, 17, 25, 18, 19), name = "Test Study") +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 9, face = "bold"),
        legend.position = "right",
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9),
        plot.title = element_text(size = 12, hjust = 0.5)
      ) +
      labs(y = "Matthews Correlation Coefficient", title = classifier)
    
    plot_list[[as.character(classifier)]] <- p
  }
  
  # Arrange and save
  cat("Arranging plots in grid\n")
  final_plot <- gridExtra::grid.arrange(grobs = unname(plot_list), ncol = 2)
  
  cat("Saving plot to:", output_file, "\n")
  ggplot2::ggsave(filename = output_file, plot = final_plot, width = width, height = height, dpi = dpi, units = "in", bg = "white")
  
  cat("Plot saved successfully!\n")
}
