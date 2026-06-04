# generate_relative_plot_aggregated.R
# Generate relative performance plot with data aggregated across n_datasets
# One panel per classifier, p-values calculated across all n_datasets

#' Generate aggregated relative performance plot
#' @param mxe_data Data frame with MCC values
#' @param top_adjuster Name of top performing adjuster (not used, kept for compatibility)
#' @param top_adjuster_label Formatted label (not used, kept for compatibility)
#' @param output_file Path to save plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Plot resolution
generate_relative_plot_aggregated <- function(mxe_data, top_adjuster, top_adjuster_label, 
                                             output_file, width = 20, height = 16, dpi = 300) {
  
  cat("\nCreating aggregated relative performance plot (one panel per classifier)...\n")
  
  # Get classifiers with data
  classifiers_with_data <- unique(mxe_data$classifier_label[!is.na(mxe_data$classifier_label)])
  
  # Create individual plots for each classifier
  relative_plot_list <- list()
  all_sig_results <- data.frame()
  
  for (classifier in classifiers_with_data) {
    # Get data for this classifier
    classifier_data <- mxe_data[mxe_data$classifier_label == classifier, ]
    
    if (nrow(classifier_data) == 0) {
      cat(sprintf("Skipping %s - no data\n", classifier))
      next
    }
    
    # Get classifier-specific ordering (within_study_cv will be first)
    classifier_adjuster_order <- get_classifier_ordering(mxe_data, classifier)
    classifier_specific_labels <- create_adjuster_labels(classifier_adjuster_order)
    reference_adjuster <- "within_study_cv"
    reference_label <- "Within-Study CV"
    
    cat(sprintf("Classifier %s: using %s as reference baseline\n", classifier, reference_label))
    
    # Calculate performance gap (negative = worse than within-study CV baseline)
    raw_data <- classifier_data %>%
      select(adjuster, classifier, classifier_label, n_datasets, test_study, value) %>%
      mutate(condition_id = paste(classifier, n_datasets, test_study, sep = "_")) %>%
      group_by(condition_id) %>%
      mutate(baseline_value = value[adjuster == reference_adjuster]) %>%
      ungroup() %>%
      mutate(relative_value = value - baseline_value) %>%
      select(-baseline_value, -condition_id)
    
    # Apply ordering
    raw_data$adjuster_label <- factor(raw_data$adjuster,
                                      levels = classifier_adjuster_order$adjuster,
                                      labels = classifier_specific_labels)
    
    # Statistical testing: one p-value per adjuster (pooling across all n_datasets)
    pval_data <- data.frame()
    
    for (adj_label in classifier_specific_labels[-1]) {  # Skip first (within_study_cv)
      adj_data <- raw_data[raw_data$adjuster_label == adj_label, ]
      
      if (nrow(adj_data) > 1) {
        # One-sample t-test: is relative_value significantly < 0?
        # (negative means this adjuster is worse than within-study CV baseline)
        t_result <- t.test(adj_data$relative_value, mu = 0, alternative = "less")
        
        pval_data <- rbind(pval_data, data.frame(
          adjuster_label = adj_label,
          p_value = t_result$p.value,
          n_obs = nrow(adj_data),
          mean_diff = mean(adj_data$relative_value, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Apply Bonferroni correction within this classifier
    if (nrow(pval_data) > 0) {
      cat(sprintf("  Applying Bonferroni correction across %d comparisons for %s\n", 
                  nrow(pval_data), classifier))
      pval_data$p_adj <- p.adjust(pval_data$p_value, method = "bonferroni")
      
      # Add significance labels
      pval_data$sig_label <- sapply(pval_data$p_adj, function(p) {
        if (p < 0.001) return("***")
        if (p < 0.01) return("**")
        if (p < 0.05) return("*")
        return("")
      })
      
      # Add classifier info
      pval_data$classifier <- classifier
      pval_data$reference_adjuster <- reference_label
      
      all_sig_results <- rbind(all_sig_results, pval_data)
    }
    
    # Calculate means
    mean_data <- raw_data %>%
      group_by(adjuster_label) %>%
      summarise(mean_value = mean(relative_value, na.rm = TRUE), .groups = "drop")
    
    # Calculate y-axis limits
    local_y_max <- max(raw_data$relative_value, na.rm = TRUE)
    local_y_min <- min(raw_data$relative_value, na.rm = TRUE)
    local_y_range <- local_y_max - local_y_min
    y_limits <- c(local_y_min - 0.05 * local_y_range, 
                  local_y_max + 0.15 * local_y_range)
    
    # Create plot (no faceting by dataset_label)
    p <- ggplot(raw_data, aes(x = adjuster_label, y = relative_value)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6, alpha = 0.7) +
      geom_segment(data = mean_data, 
                  aes(x = as.numeric(adjuster_label) - 0.3, 
                      xend = as.numeric(adjuster_label) + 0.3,
                      y = mean_value, yend = mean_value),
                  color = "gray40", linewidth = 0.8) +
      geom_point(aes(color = test_study, shape = as.factor(n_datasets)), size = 3, alpha = 0.5) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0, 0))) +
      scale_color_brewer(palette = "Set2", name = "Test Study") +
      scale_shape_manual(values = c(15, 16, 17, 18), name = "# Studies") +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
      ) +
      labs(
        y = sprintf("Performance Gap\n(vs. %s)", reference_label),
        title = classifier
      )
    
    # Add p-value annotations
    if (nrow(pval_data) > 0) {
      annotation_y <- local_y_max + 0.08 * local_y_range
      pval_data$y_pos <- annotation_y
      
      p <- p + geom_text(data = pval_data,
                        aes(x = adjuster_label, y = y_pos, label = sig_label),
                        inherit.aes = FALSE,
                        size = 5, fontface = "bold", vjust = 0)
    }
    
    relative_plot_list[[classifier]] <- p
  }
  
  # Arrange and save
  cat("Arranging aggregated relative plots in grid\n")
  relative_final_plot <- grid.arrange(grobs = unname(relative_plot_list), ncol = 2)
  
  cat("Saving aggregated relative performance plot to:", output_file, "\n")
  ggsave(filename = output_file, plot = relative_final_plot, 
         width = width, height = height, dpi = dpi, units = "in", bg = "white")
  
  cat("Aggregated relative performance plot saved successfully!\n")
  
  # Save significance results
  if (nrow(all_sig_results) > 0) {
    sig_output <- sub("\\.png$", "_significance.csv", output_file)
    write.csv(all_sig_results, sig_output, row.names = FALSE)
    cat("Significance results saved to:", sig_output, "\n")
  }
}
