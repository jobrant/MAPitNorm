#' @title Internal Visualization Utility Functions
#' @description This file contains internal utility functions for visualization.
#' @name visualization-utils
NULL

utils::globalVariables(c(".data"))

# Data Preparation Functions ----------------------------------------------

#' Prepare and sample methylation data for visualization
#'
#' @param data_list List of methylation data at different stages
#' @param stages Vector of stage names
#' @param sample_size Number of sites to sample (NULL for all)
#' @param seed Random seed for reproducibility
#' @return Data frame for plotting
#' @keywords internal
.sample_data_for_visualization <- function(data_list, stages, sample_size = NULL, seed = 42) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Initialize empty data frame to store results
  all_data <- data.frame()

  # Process each stage
  for (i in seq_along(data_list)) {
    stage_name <- stages[i]
    stage_data <- data_list[[i]]

    # Check if we're dealing with grouped data (list of lists)
    is_grouped <- all(sapply(stage_data, is.list))

    if (is_grouped) {
      # Process grouped data
      for (group_name in names(stage_data)) {
        group_data <- stage_data[[group_name]]

        for (sample_name in names(group_data)) {
          sample_df <- group_data[[sample_name]]

          # Ensure we have rate column
          if (!("rate" %in% names(sample_df))) {
            if ("mc" %in% names(sample_df) && "cov" %in% names(sample_df)) {
              sample_df$rate <- sample_df$mc / sample_df$cov
            } else {
              warning(paste("Stage", stage_name, "Group", group_name, "Sample",
                            sample_name, "doesn't have rate column"))
              next
            }
          }

          # Extract needed columns and add metadata
          temp_df <- data.frame(
            stage = stage_name,
            group = group_name,
            sample = sample_name,
            rate = sample_df$rate,
            cov = if ("cov" %in% names(sample_df)) sample_df$cov else NA,
            row_id = seq_along(sample_df$rate)
          )

          # Add to result
          all_data <- rbind(all_data, temp_df)
        }
      }
    } else {
      # Process flat data (list of data frames)
      for (sample_name in names(stage_data)) {
        sample_df <- stage_data[[sample_name]]

        # Ensure we have rate column
        if (!("rate" %in% names(sample_df))) {
          if ("mc" %in% names(sample_df) && "cov" %in% names(sample_df)) {
            sample_df$rate <- sample_df$mc / sample_df$cov
          } else {
            warning(paste("Stage", stage_name, "Sample", sample_name,
                          "doesn't have rate column"))
            next
          }
        }

        # Extract needed columns and add metadata
        temp_df <- data.frame(
          stage = stage_name,
          group = "unknown",
          sample = sample_name,
          rate = sample_df$rate,
          cov = if ("cov" %in% names(sample_df)) sample_df$cov else NA,
          row_id = seq_along(sample_df$rate)
        )

        # Add to result
        all_data <- rbind(all_data, temp_df)
      }
    }
  }

  # Sample if requested and we have data
  if (!is.null(sample_size) && nrow(all_data) > 0 && sample_size < nrow(all_data)) {
    # Sample by unique combination of row_id
    unique_ids <- unique(all_data$row_id)
    sampled_ids <- sample(unique_ids, min(sample_size, length(unique_ids)))
    all_data <- all_data[all_data$row_id %in% sampled_ids, ]
  }

  # Make stage a factor with correct order
  all_data$stage <- factor(all_data$stage, levels = stages)

  return(all_data)
}


# Individual Plot Creation Functions --------------------------------------

#' Create density plot of methylation rates
#' @param data Prepared data frame
#' @param stages Vector of stage names
#' @param theme ggplot2 theme to apply
#' @param output_dir Output directory
#' @return List with plot and report file path
#' @keywords internal
.create_density_plot <- function(data, stages, theme = theme_minimal(), output_dir) {
  # Ensure stages is a factor with correct levels for ordering
  data$stage <- factor(data$stage, levels = stages)

  # Create plot with facets by stage
  p <- ggplot(data, aes(x = rate, color = sample)) +
    geom_density(alpha = 0.7) +
    facet_wrap(~stage) +  # Each stage gets its own panel
    labs(title = "Methylation Rate Distributions",
         x = "Methylation Rate",
         y = "Density") +
    theme

  report_file <- .save_plot_report(p, "density", output_dir, length(stages))

  return(list(plot = p, report_file = report_file))
}


#' Create QQ plot for methylation rate comparison
#'
#' @param data Prepared data frame with methylation information
#' @param stages Vector of stage names
#' @param theme ggplot theme to apply
#' @param output_dir Output directory
#' @return List with plot and report file path
#' @keywords internal
.create_qq_plot <- function(data, stages, theme = theme_minimal(), output_dir) {
  # Try to create hexbin plots for all sample pairs
  hexbin_plots <- tryCatch({
    .create_qq_hexbin_plot(data, stages, bins = 50, theme)
  }, error = function(e) {
    message("Error creating hexbin plots: ", e$message)
    return(NULL)
  })

  if (is.null(hexbin_plots)) {
    p <- .create_simple_qq_plot(data, stages, theme)
    report_file <- .save_plot_report(p, "qq", output_dir, 1)
    return(list(plot = p, report_file = report_file))
  }

  report_file <- file.path(output_dir, "qq_plots_report.pdf")
  tryCatch({
    pdf(report_file, width = 10, height = 8)
    grid::grid.newpage()
    grid::grid.text("QQ Plots Report",
                    x = 0.5, y = 0.8,
                    gp = grid::gpar(fontsize = 24, fontface = "bold"))
    grid::grid.text(paste("Total plots:", length(hexbin_plots)),
                    x = 0.5, y = 0.6,
                    gp = grid::gpar(fontsize = 16))
    for (i in seq_along(hexbin_plots)) {
      print(hexbin_plots[[i]])
    }
    dev.off()
    message("QQ plots saved to ", report_file)
  }, error = function(e) {
    message("Error saving QQ report: ", e$message)
    report_file <- NULL
  })

  max_plots <- min(4, length(hexbin_plots))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- patchwork::wrap_plots(hexbin_plots[1:max_plots], ncol = 2)
  } else {
    p <- hexbin_plots[[1]]
  }
  return(list(plot = p, report_file = report_file))
}


#' Create hexbin QQ plot for sample comparison
#'
#' @param data Prepared data frame with methylation information
#' @param stages Vector of stage names
#' @param bins Number of hexagonal bins
#' @param theme ggplot theme to apply
#' @return ggplot object
#' @keywords internal
.create_qq_hexbin_plot <- function(data, stages, bins = 50, theme = theme_minimal()) {
  # Get all possible sample pairs per group and stage
  plot_list <- list()

  # For each stage, create QQ plots for all sample pairs within each group
  for (stage_name in levels(data$stage)) {
    stage_data <- subset(data, stage == stage_name)

    # Process each group
    for (group_name in unique(stage_data$group)) {
      group_data <- subset(stage_data, group == group_name)
      samples <- unique(group_data$sample)

      # Skip if fewer than 2 samples
      if (length(samples) < 2) {
        message(sprintf("Skipping QQ plot for group %s in stage %s: fewer than 2 samples",
                        group_name, stage_name))
        next
      }

      # Create all pairwise combinations
      pairs <- combn(samples, 2, simplify = FALSE)

      # For each pair, create a QQ plot
      for (pair in pairs) {
        pair_name <- paste(pair, collapse = " vs ")

        # Get data for each sample in the pair
        sample1_data <- subset(group_data, sample == pair[1])
        sample2_data <- subset(group_data, sample == pair[2])

        # Ensure we have valid data for both samples
        if (nrow(sample1_data) == 0 || nrow(sample2_data) == 0) {
          message(sprintf("Missing data for %s in %s, %s",
                          pair_name, group_name, stage_name))
          next
        }

        # Create a data frame with matching row_ids
        common_rows <- intersect(sample1_data$row_id, sample2_data$row_id)

        if (length(common_rows) < 10) {
          message(sprintf("Insufficient matching data points for %s in %s, %s",
                          pair_name, group_name, stage_name))
          next
        }

        # Create a data frame with the rates from both samples
        plot_df <- data.frame(
          x = sample1_data$rate[match(common_rows, sample1_data$row_id)],
          y = sample2_data$rate[match(common_rows, sample2_data$row_id)]
        )

        # Remove any rows with NA
        plot_df <- na.omit(plot_df)

        # Check again if we have enough data
        if (nrow(plot_df) < 10) {
          message(sprintf("Insufficient valid data points for %s in %s, %s",
                          pair_name, group_name, stage_name))
          next
        }

        # Check for constant values which cause hexbin to fail
        x_range <- range(plot_df$x, na.rm = TRUE)
        y_range <- range(plot_df$y, na.rm = TRUE)

        if (diff(x_range) <= 0 || diff(y_range) <= 0) {
          message(sprintf("Zero range detected in data for %s in %s, %s",
                          pair_name, group_name, stage_name))
          next
        }

        # Create the plot with explicit mapping
        p <- ggplot(plot_df, aes(x = x, y = y)) +
          geom_hex(bins = bins) +
          geom_abline(color = "red", linetype = "dashed") +
          labs(title = paste(group_name, "-", stage_name),
               subtitle = pair_name,
               x = pair[1],
               y = pair[2]) +
          ggplot2::scale_fill_viridis_c(option = "magma") +
          theme

        # Add to list
        plot_id <- paste(stage_name, group_name, paste(pair, collapse = "_"), sep = "_")
        plot_list[[plot_id]] <- p
      }
    }
  }

  # Check if we created any plots
  if (length(plot_list) == 0) {
    message("No QQ plots could be created with the provided data")
    return(NULL)
  }

  return(plot_list)
}


#' Create simple QQ plot as fallback
#' @param data Prepared data frame
#' @param stages Vector of stage names
#' @param theme ggplot theme
#' @return ggplot object
#' @keywords internal
.create_simple_qq_plot <- function(data, stages, theme = theme_minimal()) {
  # Calculate average rates by group and stage
  avg_data <- aggregate(rate ~ group + stage, data, mean)

  # Create a simple comparison plot between stages
  if (length(stages) >= 2) {
    # Compare first two stages
    stage1 <- stages[1]
    stage2 <- stages[2]

    comp_data <- reshape2::dcast(avg_data, group ~ stage, value.var = "rate")

    p <- ggplot(comp_data, aes(x = .data[[stage1]], y = .data[[stage2]])) +
      geom_point(aes(color = group), size = 3) +
      geom_abline(color = "red", linetype = "dashed") +
      labs(title = paste("Comparison:", stage1, "vs", stage2),
           x = paste("Rate in", stage1),
           y = paste("Rate in", stage2)) +
      theme

    return(p)
  } else {
    # If only one stage, show rates by group
    p <- ggplot(avg_data, aes(x = group, y = rate, fill = group)) +
      geom_bar(stat = "identity") +
      facet_wrap(~stage) +
      labs(title = "Average Methylation Rate by Group",
           x = "Group",
           y = "Average Rate") +
      theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    return(p)
  }
}

# Determine Layout --------------------------------------------------------

#' Determine optimal layout for a set of plots
#'
#' @param plot_types Character vector of plot types
#' @return A layout specification (ncol, nrow)
#' @keywords internal
.get_optimal_layout <- function(plot_types) {
  n_plots <- length(plot_types)

  # Special cases for specific plot combinations
  if ("density" %in% plot_types && "qq" %in% plot_types && n_plots == 2) {
    # Density and QQ plots look better side by side
    return(list(ncol = 2, nrow = 1))
  }

  if ("boxplot" %in% plot_types && n_plots > 2) {
    # Boxplots often need more horizontal space
    return(list(ncol = 1, nrow = n_plots))
  }

  # Default layouts based on number of plots
  if (n_plots <= 2) {
    return(list(ncol = n_plots, nrow = 1))
  } else if (n_plots <= 4) {
    return(list(ncol = 2, nrow = ceiling(n_plots / 2)))
  } else {
    # Try to make grid as square as possible
    ncol <- round(sqrt(n_plots))
    nrow <- ceiling(n_plots / ncol)
    return(list(ncol = ncol, nrow = nrow))
  }
}


# Saving Functions --------------------------------------------------------
#' Save plot to PDF report
#' @param plot ggplot object
#' @param plot_type Plot type
#' @param output_dir Output directory
#' @param n_facets Number of facets (for layout)
#' @return File path of saved report
#' @keywords internal
.save_plot_report <- function(plot, plot_type, output_dir, n_facets) {
  report_file <- file.path(output_dir, paste0(plot_type, "_report.pdf"))
  tryCatch({
    pdf(report_file, width = 10, height = 8)
    grid::grid.newpage()
    grid::grid.text(paste(toupper(substr(plot_type, 1, 1)), substr(plot_type, 2, nchar(plot_type)), "Plot Report"),
                    x = 0.5, y = 0.8,
                    gp = grid::gpar(fontsize = 24, fontface = "bold"))
    grid::grid.text(paste("Generated on:", Sys.Date()),
                    x = 0.5, y = 0.6,
                    gp = grid::gpar(fontsize = 16))
    print(plot)
    dev.off()
    message(plot_type, " report saved to ", report_file)
    return(report_file)
  }, error = function(e) {
    message("Error saving ", plot_type, " report: ", e$message)
    return(NULL)
  })
}

