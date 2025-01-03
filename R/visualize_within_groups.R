#' Visualize Within-Group Normalization Effects
#'
#' @param original_data List of original methylation data
#' @param coverage_normalized List of coverage-normalized data
#' @param fully_normalized List of fully normalized data
#' @param group_names Character vector of group names
#' @param sample_size Number of sites to sample for visualization (default = 100000)
#' @param save_plot Logical indicating whether to save plots (default = FALSE)
#' @param output_dir Directory for saving plots (default = ".")
#' @param plot_width Width of saved plot (default: 15)
#' @param plot_height Height of saved plot (default: 12)
#' @param show_correlation_heatmaps Logical indicating whether to show correlation heatmaps (default: FALSE)
#'
#' @return List containing plot data and plots
#'
#' @importFrom data.table setDT copy :=
#' @import ggplot2
#' @import cowplot
#'
#' @export
visualize_within_groups <- function(original_data,
                                    coverage_normalized,
                                    fully_normalized,
                                    group_names,
                                    sample_size = 100000,
                                    save_plot = FALSE,
                                    output_dir = ".",
                                    plot_width = 15,
                                    plot_height = 12,
                                    show_correlation_heatmaps = FALSE) {

  # Add initial error checking
  if (!all(sapply(list(original_data, coverage_normalized, fully_normalized), is.list))) {
    stop("All input data must be lists")
  }

  if(plot_width <= 0 || plot_height <= 0) {
    stop("Plot dimensions must be positive numbers")
  }

  # Sample data function definition
  sample_data <- function(data_list, stage) {
    set.seed(42)
    n_sites <- nrow(data_list[[1]][[1]])
    sample_idx <- sample(n_sites, min(sample_size, n_sites))

    dt <- rbindlist(
      lapply(names(data_list), function(group) {
        rbindlist(
          lapply(seq_along(data_list[[group]]), function(i) {
            data.table(
              group = group,
              replicate = paste0("Rep", i),
              coverage = data_list[[group]][[i]]$cov[sample_idx],
              rate = data_list[[group]][[i]]$rate[sample_idx],
              site_idx = sample_idx,
              stage = stage
            )
          })
        )
      })
    )

    return(dt)
  }

  # Prepare data
  cat("Preparing data...\n")
  plot_data <- rbindlist(list(
    sample_data(original_data, "Original"),
    sample_data(coverage_normalized, "Coverage Normalized"),
    sample_data(fully_normalized, "Fully Normalized")
  ))

  # Set factor levels for stages
  plot_data[, stage := factor(stage,
                              levels = c("Original",
                                         "Coverage Normalized",
                                         "Fully Normalized"))]

  # Print data summary
  print("Data summary:")
  print(plot_data[, .N, by = .(group, stage, replicate)])

  # Create basic plots
  cat("Creating basic plots...\n")
  p1 <- ggplot(plot_data, aes(x = replicate, y = log10(coverage), fill = replicate)) +
    geom_boxplot() +
    facet_grid(stage ~ group) +
    labs(title = "Coverage Distribution by Replicate",
         y = "log10(Coverage)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2 <- ggplot(plot_data, aes(x = rate, color = replicate)) +
    geom_density() +
    facet_grid(stage ~ group) +
    labs(title = "Methylation Rate Distribution by Replicate") +
    theme_minimal()

  # Create correlation plots
  cat("Creating correlation plots...\n")
  p3 <- p4 <- ggplot() + theme_minimal() + ggtitle("No correlation analysis available")  # Default empty plots

  tryCatch({
    cor_analysis <- analyze_correlation_changes(plot_data, n_bins = 10)
    if (!is.null(cor_analysis$plots)) {
      p3 <- cor_analysis$plots$correlations
      p4 <- cor_analysis$plots$distribution
    }
  }, error = function(e) {
    cat("Error in correlation analysis:", e$message, "\n")
    p3 <- p4 <- ggplot() + theme_minimal() +
      ggtitle("Error in correlation analysis")
  })

  # Create correlation heatmaps and combine plots
  if(show_correlation_heatmaps) {
    tryCatch({
      cat("Creating correlation heatmaps...\n")
      orig_plots <- .create_cor_plot(plot_data, "Original")
      norm_plots <- .create_cor_plot(plot_data, "Fully Normalized")

      # Close any open graphics devices
      while (dev.cur() > 1) dev.off()

      if(save_plot) {
        pdf_file <- file.path(output_dir, "within_group_normalization.pdf")
        pdf(pdf_file, width = plot_width, height = plot_height)

        # Basic plots and correlation analysis
        print(cowplot::plot_grid(p1, p2, p3, p4,
                                 ncol = 2,
                                 labels = c("A", "B", "C", "D")))

        # Correlation heatmaps
        for(p in orig_plots) print(p)
        for(p in norm_plots) print(p)

        dev.off()

        if(file.exists(pdf_file)) {
          cat(sprintf("Plots saved to %s\n", pdf_file))
        } else {
          warning("PDF file was not created successfully")
        }
      }

      # For display and return
      combined_plot <- cowplot::plot_grid(p1, p2, p3, p4,
                                          ncol = 2,
                                          labels = c("A", "B", "C", "D"))

      # Store heatmap plots
      p5 <- orig_plots
      p6 <- norm_plots

    }, error = function(e) {
      cat("Error in creating correlation heatmaps:", e$message, "\n")
      cat("Falling back to basic plots...\n")
      combined_plot <- cowplot::plot_grid(p1, p2, p3, p4,
                                          ncol = 2,
                                          labels = c("A", "B", "C", "D"))
    })
  } else {
    # Create basic combined plot
    combined_plot <- cowplot::plot_grid(p1, p2, p3, p4,
                                        ncol = 2,
                                        labels = c("A", "B", "C", "D"))

    # Save if requested
    if(save_plot) {
      filename <- file.path(output_dir, "within_group_normalization.pdf")
      ggsave(filename, combined_plot, width = plot_width, height = plot_height)
      cat(sprintf("Plots saved to %s\n", filename))
    }
  }

  return(list(
    plot_data = plot_data,
    combined_plot = combined_plot,
    individual_plots = if(show_correlation_heatmaps && exists("p5") && exists("p6")) {
      list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6)
    } else {
      list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
    }
  ))
}


#' Internal function to create correlation heatmaps
#' @keywords internal
.create_cor_plot <- function(data, stage_name) {
  # Initial diagnostics
  cat(sprintf("\nCreating correlation plot for stage: %s\n", stage_name))
  cat("Number of groups:", length(unique(data$group)), "\n")
  cat("Groups:", paste(unique(data$group), collapse=", "), "\n")

  # Filter data for specified stage first
  data <- data[stage == stage_name]

  # Process each group separately
  plot_list <- lapply(unique(data$group), function(g) {
    cat(sprintf("Processing group %s...\n", g))

    # Get data for this group
    group_data <- data[group == g]

    # Create wide format matrix
    wide_data <- dcast(group_data, site_idx ~ replicate, value.var = "rate")
    wide_data[, site_idx := NULL]

    # Calculate correlation
    cor_mat <- cor(wide_data)

    # Convert to long format
    cor_dt <- as.data.table(cor_mat, keep.rownames = "rn")
    cor_long <- melt(cor_dt,
                     id.vars = "rn",
                     variable.name = "variable",
                     value.name = "value")

    # Create individual plot for this group
    p <- ggplot(cor_long, aes(x = rn, y = variable, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(limits = c(-1, 1),
                           low = "blue", mid = "white", high = "red",
                           midpoint = 0) +
      geom_text(aes(label = sprintf("%.2f", value)),
                color = ifelse(cor_long$value > 0.5, "black", "white")) +
      labs(title = paste(g, "-", stage_name),
           x = "Replicate", y = "Replicate") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    return(p)
  })

  # Return list of plots
  return(plot_list)
}


#' Analyze and visualize rate distributions and correlations
#' @param plot_data Data table containing methylation data
#' @param n_bins Number of bins for rate analysis (default = 10 for deciles)
#' @param show_diagnostics Logical to print detailed bin information
#' @return List containing analysis results and plots
#' @keywords internal
analyze_correlation_changes <- function(plot_data, n_bins = 10, show_diagnostics = TRUE) {
  # Print initial diagnostics
  cat("Analyzing correlations...\n")
  cat("Number of rows in plot_data:", nrow(plot_data), "\n")
  cat("Groups present:", paste(unique(plot_data$group), collapse = ", "), "\n")
  cat("Stages present:", paste(levels(plot_data$stage), collapse = ", "), "\n")

  # Ensure factors
  plot_data[, group := factor(group)]
  plot_data[, stage := factor(stage, levels = c("Original",
                                                "Coverage Normalized",
                                                "Fully Normalized"))]

  # Calculate rate breaks and create bins
  probs <- seq(0, 1, length.out = n_bins + 1)
  rate_breaks <- quantile(plot_data$rate, probs = probs, na.rm = TRUE)

  # Create bin labels
  bin_labels <- if(n_bins == 10) {
    paste0("D", 1:10)  # Deciles
  } else if(n_bins == 5) {
    paste0("Q", 1:5)   # Quintiles
  } else {
    paste0("Bin", 1:n_bins)
  }

  # Print bin information if requested
  if(show_diagnostics) {
    cat("\nRate distribution by bin:\n")
    for(i in 1:n_bins) {
      bin_data <- plot_data[rate >= rate_breaks[i] & rate <= rate_breaks[i+1]]
      cat(sprintf("%s: %.4f to %.4f (n=%d, %.1f%% of data)\n",
                  bin_labels[i],
                  rate_breaks[i],
                  rate_breaks[i+1],
                  nrow(bin_data),
                  100 * nrow(bin_data) / nrow(plot_data)))
    }
  }

  # Assign bins to data
  plot_data[, rate_bin := cut(rate,
                              breaks = rate_breaks,
                              labels = bin_labels,
                              include.lowest = TRUE)]

  # Create distribution visualization
  p_dist <- ggplot(plot_data, aes(x = rate)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    geom_vline(xintercept = rate_breaks,
               color = "red",
               linetype = "dashed") +
    geom_text(data = data.frame(
      x = (rate_breaks[-1] + rate_breaks[-length(rate_breaks)])/2,
      y = 0,
      label = bin_labels
    ),
    aes(x = x, y = y, label = label),
    vjust = -0.5,
    angle = 90) +
    facet_wrap(~stage) +
    labs(title = "Rate Distribution with Bin Breaks",
         x = "Methylation Rate",
         y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Calculate correlations within each bin
  bin_cors <- plot_data[!is.na(rate_bin), {
    if(show_diagnostics) {
      cat(sprintf("\nProcessing group %s, stage %s, bin %s (n=%d)\n",
                  unique(group), unique(stage), unique(rate_bin), .N))
    }

    wide_data <- dcast(.SD, site_idx ~ replicate, value.var = "rate")
    wide_data[, site_idx := NULL]

    if(ncol(wide_data) < 2) return(NULL)

    cor_mat <- cor(wide_data, use = "pairwise.complete.obs")
    cor_vals <- cor_mat[upper.tri(cor_mat)]

    if(length(cor_vals) == 0) return(NULL)

    .(mean_cor = mean(cor_vals, na.rm = TRUE),
      sd_cor = sd(cor_vals, na.rm = TRUE),
      n_sites = .N,
      min_rate = min(rate),
      max_rate = max(rate))
  }, by = .(group, stage, rate_bin)]

  # Create correlation plots
  p1 <- ggplot(bin_cors, aes(x = rate_bin, y = mean_cor, color = stage)) +
    geom_point(size = 3) +
    geom_line(aes(group = stage)) +
    geom_errorbar(aes(ymin = pmax(mean_cor - sd_cor, 0),
                      ymax = pmin(mean_cor + sd_cor, 1)),
                  width = 0.2) +
    facet_wrap(~group) +
    labs(title = paste("Correlation by Rate", if(n_bins == 10) "Decile" else "Bin"),
         x = if(n_bins == 10) "Decile" else "Bin",
         y = "Mean Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Return results
  return(list(
    summary = bin_cors,
    rate_breaks = rate_breaks,
    plots = list(
      distribution = p_dist,
      correlations = p1
    )
  ))
}
