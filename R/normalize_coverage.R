# File: R/normalize_coverage.R

#' Normalize Coverage Data
#'
#' Normalizes methylation coverage data across replicates and groups by scaling
#' coverage values to the mean total coverage across all samples.
#'
#' @param data_list A list of lists containing methylation data frames. Each inner list
#'        represents a group and contains one or more replicate data frames.
#' @param group_names A character vector of group names.
#'
#' @return A list of lists with the same structure as the input, containing normalized
#'         coverage data.
#'
#' @details Each data frame in the input list should contain columns:
#'          \itemize{
#'            \item cov: coverage values
#'            \item mc: methylation counts
#'          }
#'          The function calculates scaling factors based on total coverage and
#'          normalizes both coverage and methylation counts accordingly.
#'
#' @examples
#' \dontrun{
#' split_data <- split_by_groups(gch_samples, group_names)
#' normalized_coverage <- normalize_coverage(split_data, group_names)
#' }
#' @import data.table
#' @export

normalize_coverage <- function(data_list, group_names) {
  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (!all(sapply(data_list, is.list))) stop("data_list must be a list of lists")
  if (!all(sapply(unlist(data_list, recursive = FALSE),
                  function(df) all(c("cov", "mc") %in% names(df))))) {
    stop("All data frames must contain 'cov' and 'mc' columns")
  }

  # Get all data frames in a flat list
  all_dfs <- unlist(data_list, recursive = FALSE)

  # Calculate scaling factors
  totals <- sapply(all_dfs, function(df) sum(df$cov))
  avg_total <- mean(totals)
  sf_total <- totals / avg_total

  print("Scaling factors:")
  print(sf_total)

  # Normalize each group
  normalized_groups <- lapply(names(data_list), function(group_name) {
    group_data <- data_list[[group_name]]
    group_sf <- sf_total[names(group_data)]

    # Print diagnostics
    print(paste("Processing group:", group_name))
    print("Group scaling factors:")
    print(group_sf)

    normalize_group(group_data, group_sf)
  })

  names(normalized_groups) <- names(data_list)
  return(normalized_groups)
}
