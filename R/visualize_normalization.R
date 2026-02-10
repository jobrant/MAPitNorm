#' Create Methylation Normalization Visualization Report
#'
#' Generates a multi-page PDF report comparing methylation data across
#' normalization stages. This is a convenience wrapper around the
#' individual plot functions (\code{\link{plot_density}},
#' \code{\link{plot_qq}}, \code{\link{plot_ma}}, \code{\link{plot_pca}}).
#'
#' @param data_list Named list where each element contains methylation data
#'   at a different normalization stage (e.g.,
#'   \code{list(raw = raw_data, normalized = norm_data)}).
#' @param stages Character vector of stage names corresponding to
#'   \code{data_list} elements. Default uses \code{names(data_list)}.
#' @param plots Character vector of plot types to include. Options are
#'   \code{"density"}, \code{"qq"}, \code{"ma"}, and \code{"pca"}.
#'   Default is \code{c("density", "qq")}.
#' @param sample_size Number of sites to randomly sample for visualization.
#'   Default is 100000. Use \code{NULL} for all sites.
#' @param output_file File path to save the combined PDF report.
#'   If \code{NULL}, no file is saved.
#' @param output_dir Directory for individual plot reports.
#'   Defaults to the directory of \code{output_file}, or the working
#'   directory if \code{output_file} is also \code{NULL}.
#' @param theme ggplot2 theme to apply. Default is
#'   \code{\link[ggplot2]{theme_minimal}}.
#'
#' @return A list containing:
#'   \item{plots}{Named list of ggplot objects}
#'   \item{combined}{Combined patchwork plot (if output_file specified)}
#'   \item{report_files}{Named list of paths to individual PDF reports}
#'
#' @examples
#' \dontrun{
#' report <- visualize_normalization(
#'   data_list = list(
#'     "Before" = split_data,
#'     "After"  = normalized_data
#'   ),
#'   plots = c("density", "ma", "pca"),
#'   output_file = "normalization_report.pdf"
#' )
#' }
#'
#' @seealso \code{\link{plot_density}}, \code{\link{plot_qq}},
#'   \code{\link{plot_ma}}, \code{\link{plot_pca}}
#'
#' @export
visualize_normalization <- function(data_list,
                                    stages = names(data_list),
                                    plots = c("density", "qq"),
                                    sample_size = 100000,
                                    output_file = NULL,
                                    output_dir = NULL,
                                    theme = ggplot2::theme_minimal()) {
  
  # --- Validation ---
  if (!is.list(data_list) || length(data_list) == 0) {
    stop("'data_list' must be a non-empty list of methylation data.")
  }
  
  valid_plot_types <- c("density", "qq", "ma", "pca")
  invalid <- setdiff(plots, valid_plot_types)
  if (length(invalid) > 0) {
    stop("Invalid plot type(s): ", paste(invalid, collapse = ", "),
         ". Choose from: ", paste(valid_plot_types, collapse = ", "))
  }
  
  if (is.null(stages) || length(stages) == 0) {
    stages <- if (!is.null(names(data_list)) && all(nchar(names(data_list)) > 0)) {
      names(data_list)
    } else {
      paste("Stage", seq_along(data_list))
    }
  }
  
  if (length(stages) != length(data_list)) {
    stop("Length of 'stages' must match length of 'data_list'.")
  }
  
  # --- Output directory ---
  if (is.null(output_dir)) {
    output_dir <- if (!is.null(output_file)) dirname(output_file) else getwd()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  
  # --- Generate plots ---
  plot_objects <- list()
  report_files <- list()
  
  for (plot_type in plots) {
    p <- tryCatch({
      switch(
        plot_type,
        density = plot_density(data_list, stages, sample_size, theme),
        qq      = plot_qq(data_list, stages, sample_size, theme = theme),
        ma      = plot_ma(data_list, stages, sample_size, theme = theme),
        pca     = plot_pca(data_list, stages, theme = theme)
      )
    }, error = function(e) {
      message("Failed to create ", plot_type, " plot: ", e$message)
      NULL
    })
    
    if (!is.null(p)) {
      plot_objects[[plot_type]] <- p
      rf <- .save_plot_report(p, plot_type, output_dir, length(stages))
      if (!is.null(rf)) report_files[[plot_type]] <- rf
    }
  }
  
  if (length(plot_objects) == 0) {
    warning("No plots were successfully created.")
    return(list(plots = list(), combined = NULL, report_files = list()))
  }
  
  # --- Combined report ---
  combined_report <- NULL
  if (!is.null(output_file) && length(plot_objects) > 0) {
    tryCatch({
      layout <- .get_optimal_layout(names(plot_objects))
      
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined_report <- patchwork::wrap_plots(
          plot_objects, ncol = layout$ncol, nrow = layout$nrow
        )
      } else {
        combined_report <- plot_objects[[1]]
        message("patchwork not available; combined report contains first plot only.")
      }
      
      grDevices::pdf(output_file, width = 12, height = 10)
      grid::grid.newpage()
      grid::grid.text("Methylation Normalization Report",
                      x = 0.5, y = 0.8,
                      gp = grid::gpar(fontsize = 24, fontface = "bold"))
      grid::grid.text(
        paste("Plots:", paste(names(plot_objects), collapse = ", ")),
        x = 0.5, y = 0.6,
        gp = grid::gpar(fontsize = 16))
      grid::grid.text(
        paste("Generated on:", Sys.Date()),
        x = 0.5, y = 0.4,
        gp = grid::gpar(fontsize = 14))
      print(combined_report)
      grDevices::dev.off()
      message("Combined report saved to ", output_file)
    }, error = function(e) {
      message("Error saving combined report: ", e$message)
    })
  }
  
  return(list(
    plots        = plot_objects,
    combined     = combined_report,
    report_files = report_files
  ))
}