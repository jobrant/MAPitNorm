#' Complete Methylation Data Normalization
#'
#' Performs a complete normalization workflow for methylation data,
#' handling both coverage and methylation rates.
#'
#' @param data_list List of methylation data by group
#' @param group_names Optional character vector of group names. If NULL (default),
#'   will be extracted from the metadata attached to data_list.
#' @param do_coverage_norm Logical, whether to perform coverage normalization
#' @param normalize_rates Logical, whether to perform rate normalization
#' @param coverage_between_groups Logical, normalize coverage between groups
#' @param rate_within_groups Logical, normalize rates within groups
#' @param rate_between_groups Logical, normalize rates between groups
#' @param within_alpha Alpha parameter for within-group rate preservation
#' @param between_alpha Alpha parameter for between-group rate preservation
#' @param min_coverage Minimum coverage threshold for including sites
#' @param max_quantiles Maximum number of quantiles to use for rate normalization
#' @param sites_per_quantile Target number of sites per quantile for rate normalization
#' @param diagnostics Logical, whether to output diagnostic information
#'
#' @details
#' This function provides a unified interface to the complete methylation normalization
#' workflow. It combines coverage normalization and methylation rate normalization in
#' a single call with sensible defaults.
#'
#' For advanced users the individual normalization functions
#' \code{\link{normalize_coverage}} and \code{\link{normalize_methylation_rates}}
#' can be called directly.
#'
#' @return Normalized methylation data list
#'
#' @examples
#' \dontrun{
#' # Basic usage with defaults
#' normalized_data <- normalize_methylation_data(
#'   data_list = split_data,
#'   group_names = group_names
#' )
#'
#' # Preserve biological differences between groups
#' normalized_data <- normalize_methylation_data(
#'   data_list = split_data,
#'   group_names = group_names,
#'   coverage_between_groups = FALSE,
#'   rate_between_groups = TRUE,
#'   between_alpha = 0.7
#' )
#' }
#'
#' @seealso
#' \code{\link{normalize_coverage}} for coverage normalization only
#' \code{\link{normalize_methylation_rates}} for rate normalization only
#'
#' @export
#'
normalize_methylation_data <- function(data_list,
                                       group_names = NULL,
                                       do_coverage_norm = TRUE,
                                       normalize_rates = TRUE,
                                       coverage_between_groups = FALSE,
                                       rate_within_groups = TRUE,
                                       rate_between_groups = TRUE,
                                       within_alpha = 0.3,
                                       between_alpha = 0.5,
                                       min_coverage = 5,
                                       max_quantiles = 50,
                                       sites_per_quantile = 1000,
                                       diagnostics = FALSE) {

  group_names <- extract_group_names(data_list, diagnostics)

  # Step 0: Filter by minimum coverage if specified
  if (min_coverage > 0) {
    data_list <- lapply(data_list, function(group) {
      lapply(group, function(df) {
        df[df$cov >= min_coverage, ]
      })
    })

    if (diagnostics) {
      cat("Filtered sites with coverage <", min_coverage, "\n")
    }
  }

  # Step 1: Coverage normalization
  if (do_coverage_norm) {
    data_list <- normalize_coverage(data_list,
                                    group_names = group_names,
                                    between_groups = coverage_between_groups,
                                    diagnostics = diagnostics)
  }

  # Step 2: Rate normalization
  if (normalize_rates) {
    data_list <- normalize_methylation_rates(data_list,
                                             group_names = group_names,
                                             within_groups = rate_within_groups,
                                             between_groups = rate_between_groups,
                                             within_alpha = within_alpha,
                                             between_alpha = between_alpha,
                                             max_quantiles = max_quantiles,
                                             sites_per_quantile = sites_per_quantile,
                                             diagnostics = diagnostics)
  }

  return(data_list)
}
