#' Create Methylation Normalization Visualization Report
#'
#' Generates visualization reports to compare methylation data across normalization stages.
#' Supports density, QQ, MA, boxplot, and correlation plots, with individual reports per
#' plot type and a combined report for multiple plot types.
#'
#' @param data_list List containing methylation data at different normalization stages
#' @param stages Character vector of stage names corresponding to data_list elements
#' @param plots Character vector of plot types to include
#' @param sample_size Number of sites to randomly sample for visualization (NULL for all)
#' @param output_file Optional file path to save visualization report
#' @param output_dir Directory to save additional plot reports (default is same as output_file).
#' @param theme ggplot2 theme to use apply (defaults to theme_minimual()).
#' @return A list containing:
#'   \item{plots}{List of ggplot objects}
#'   \item{combined}{Combined plot (if output_file specified)}
#'   \item{data}{Sampled data used for plotting}
#'   \item{report_files}{List of paths to individual plot reports}
#' @export
visualize_normalization <- function(data_list,
                                    stages = names(data_list),
                                    plots = c("density", "qq"),
                                    sample_size = 100000,
                                    output_file = NULL,
                                    output_dir = NULL,
                                    theme = theme_minimal()) {

  # Validation
  if (!is.list(data_list)) {
    stop("data_list must be a list of methylation data")
  }

  if (length(data_list) == 0) {
    stop("data_list is empty")
  }

  if (is.null(stages) || length(stages) == 0) {
    stages <- if (!is.null(names(data_list)) && all(names(data_list) != "")) {
      names(data_list)
    } else {
      paste("Stage", seq_along(data_list))
    }
  }

  if (length(stages) != length(data_list)) {
    stop("Length of stages must match length of data_list")
  }


  # Determine output directory
  if (is.null(output_dir)) {
    if (!is.null(output_file)) {
      output_dir <- dirname(output_file)
    } else {
      output_dir <- getwd()
    }
  }

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  # Fix any double slashes in paths
  if (!is.null(output_dir)) {
    output_dir <- gsub("//", "/", output_dir)
    output_dir <- gsub("\\\\\\\\", "\\\\", output_dir)
  }

  # Sample and prepare data
  sampled_data <- tryCatch({
    .sample_data_for_visualization(data_list, stages, sample_size)
  }, error = function(e) {
    message("Error in data preparation: ", e$message)
    return(NULL)
  })

  # Stop if data preparation failed
  if (is.null(sampled_data)) {
    stop("Failed to prepare data for visualization. Check that your data structure is correct.")
  }

  # Create requested plots
  plot_objects <- list()
  report_files <- list()

  # Create plots for each requested type
  for (plot_type in plots) {
    plot_result <- switch(
      plot_type,
      density = .create_density_plot(sampled_data, stages, theme, output_dir),
      qq = .create_qq_plot(sampled_data, stages, theme, output_dir)
      # ma = .create_ma_plot(sampled_data, stages, theme, output_dir),
      # boxplot = .create_boxplot(sampled_data, stages, theme, output_dir),
      # correlation = .create_correlation_plot(sampled_data, stages, theme, output_dir)
    )

    if (!is.null(plot_result$plot)) {
      plot_objects[[plot_type]] <- plot_result$plot
      if (!is.null(plot_result$report_file)) {
        report_files[[plot_type]] <- plot_result$report_file
      }
    } else {
      message("Failed to create ", plot_type, " plot")
    }
  }

  # Create combined report if multiple plots and output_file is specified
  combined_report <- NULL
  if (length(plot_objects) > 0 && !is.null(output_file)) {
    tryCatch({
      layout <- .get_optimal_layout(names(plot_objects))
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined_report <- patchwork::wrap_plots(plot_objects, ncol = layout$ncol, nrow = layout$nrow)
      } else {
        combined_report <- plot_objects[[1]]
        message("patchwork not available; returning first plot")
      }

      # Save combined report to PDF
      pdf(output_file, width = 12, height = 10)
      grid::grid.newpage()
      grid::grid.text("Methylation Normalization Report",
                      x = 0.5, y = 0.8,
                      gp = grid::gpar(fontsize = 24, fontface = "bold"))
      grid::grid.text(paste("Plots:", paste(names(plot_objects), collapse = ", ")),
                      x = 0.5, y = 0.6,
                      gp = grid::gpar(fontsize = 16))
      grid::grid.text(paste("Generated on:", Sys.Date()),
                      x = 0.5, y = 0.4,
                      gp = grid::gpar(fontsize = 14))
      print(combined_report)
      dev.off()
      message("Combined report saved to ", output_file)
    }, error = function(e) {
      message("Error saving combined report: ", e$message)
    })
  }

  return(list(
    plots = plot_objects,
    combined = combined_report,
    data = sampled_data,
    report_files = report_files
  ))
}

