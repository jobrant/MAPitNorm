# File: R/prepare_plot_data.R

#' Prepare Data for Plotting
#'
#' Prepare data frames for plotting normalization results.
#'
#' @param groups List of group replicate data frames.
#'
#' @param measure The measurement to prepare data for (e.g., "rate", "cov", "mc").
#'
#' @return A data frame ready for plotting.
#'
#' @import dplyr
#'
#' @import purrr
#'
#' @export
#'
prepare_data <- function(groups, measure) {
  map_dfr(seq_along(groups), function(i) {
    group_name <- paste0("M", i)
    map_dfr(groups[[i]], function(df) {
      data.frame(
        group = group_name,
        value = df[[measure]],
        rate = df$rate
      )
    })
  })
}
