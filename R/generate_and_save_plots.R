# File: R/generate_and_save_plots.R

#' Generate and Save Plots
#'
#' Generates and saves plots comparing normalized and unnormalized data.
#'
#' @param norm_data Normalized data prepared for plotting.
#'
#' @param unnorm_data Unnormalized data prepared for plotting.
#'
#' @param measure Measurement type (e.g., "rate", "cov", "mc").
#'
#' @param ylabel Label for the y-axis.
#'
#' @param filename Output filename for the saved plots.
#'
#' @import ggplot2
#'
#' @import gridExtra
#'
#' @export
#'
generate_and_save_plots <- function(norm_data, unnorm_data, measure, ylabel, filename) {
  norm_plots <- create_plots(norm_data, measure, ylabel, normalized = TRUE)
  unnorm_plots <- create_plots(unnorm_data, measure, ylabel, normalized = FALSE)

  pdf(filename, width = 15, height = 10)
  print(grid.arrange(
    unnorm_plots[[1]], unnorm_plots[[2]],
    norm_plots[[1]], norm_plots[[2]],
    ncol = 2
  ))
  dev.off()
}
