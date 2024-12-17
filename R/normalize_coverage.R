# File: R/normalize_coverage.R

#' Normalizes methylation coverage data within groups (replicates) and optionally between groups.
#'
#' @param data_list A list of lists containing methylation data frames. Each inner list
#'        represents a group and contains one or more replicate data frames.
#' @param group_names A character vector of group names.
#' @param between_groups Logical indicating whether to normalize between groups.
#'        Default is FALSE to preserve biological differences between groups.
#'
#' @details When between_groups = FALSE, each group is normalized independently,
#'          preserving potential biological differences in total accessibility between groups.
#'          When between_groups = TRUE, all samples are normalized to the same mean coverage,
#'          which might mask global accessibility changes between conditions.
#'
#' @return A list of lists with the same structure as the input, containing normalized
#'         coverage data.
#'
#' @examples
#' \dontrun{
#' split_data <- split_by_groups(gch_samples, group_names)
#' normalized_coverage <- normalize_coverage(split_data, group_names)
#' }
#' @import data.table
#' @export

normalize_coverage <- function(data_list, group_names, between_groups = FALSE) {
  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (!all(sapply(data_list, is.list))) stop("data_list must be a list of lists")
  if (!all(sapply(unlist(data_list, recursive = FALSE),
                  function(df) all(c("cov", "mc") %in% names(df))))) {
    stop("All data frames must contain 'cov' and 'mc' columns")
  }

  # Create a named vector to store scaling factors
  sf_total <- numeric()

  if(between_groups) {
    # Normalize across all samples
    for(group in names(data_list)) {
      for(sample_name in names(data_list[[group]])) {
        sf_total[sample_name] <- sum(data_list[[group]][[sample_name]]$cov)
      }
    }
    avg_total <- mean(sf_total)
    sf_total <- sf_total / avg_total

    print("Performing between-group normalization")
  } else {
    # Normalize within each group separately
    for(group in names(data_list)) {
      group_totals <- sapply(data_list[[group]], function(df) sum(df$cov))
      group_mean <- mean(group_totals)
      group_sf <- group_totals / group_mean
      sf_total[names(data_list[[group]])] <- group_sf
    }

    print("Performing within-group normalization only")
  }


  # Print diagnostics
  print("Scaling factors:")
  print(sf_total)

  # Normalize each group
  normalized_groups <- lapply(names(data_list), function(group_name) {
    group_data <- data_list[[group_name]]

    # Get scaling factors for this group using names
    group_sf <- sf_total[names(group_data)]

    print(paste("Processing group:", group_name))
    print("Group scaling factors:")
    print(group_sf)

    normalize_group(group_data, group_sf)
  })

  names(normalized_groups) <- names(data_list)
  return(normalized_groups)
}
