#' Normalize Methylation Rates
#'
#' Normalizes methylation rates using quantile-based normalization to account for
#' non-linear enzyme efficiency differences.
#'
#' @param data_list List of normalized coverage data
#' @param group_names Optional character vector of group names. If NULL (default),
#'   will be extracted from the metadata attached to data_list.
#' @param within_groups Logical indicating whether to normalize within groups (default = TRUE)
#' @param between_groups Logical indicating whether to normalize between groups (default = TRUE)
#' @param within_alpha Numeric between 0 and 1 controlling structure preservation within groups (default = 0.3)
#' @param between_alpha Numeric between 0 and 1 controlling structure preservation between groups (default = 0.5)
#' @param sites_per_quantile Target number of sites per quantile (default = 1000).
#'        For large datasets (>100k sites), this provides a good balance between
#'        normalization granularity and statistical power.
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param diagnostics Logical indicating whether to print diagnostic information
#'
#' @details This function performs a two-step normalization:
#'          1. Within-group normalization to handle technical variation between replicates
#'          2. Between-group normalization to adjust for global efficiency differences
#'          while preserving biological differences. The higher between_alpha preserves
#'          more of the biological differences between groups.
#'
#' @return List of lists containing normalized methylation data
#'
#' @importFrom data.table setDT copy :=
#' @importFrom stats quantile
#' @importFrom purrr map2
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
normalize_methylation_rates <- function(data_list,
                                        group_names = NULL,
                                        within_groups = TRUE,
                                        between_groups = TRUE,
                                        within_alpha = 0.3,
                                        between_alpha = 0.5,
                                        sites_per_quantile = 1000,
                                        max_quantiles = 50,
                                        diagnostics = TRUE) {

  # Extract group names from metadata if not provided
  if (is.null(group_names)) {
    # Check for metadata attribute
    if (is.list(data_list) && !is.null(attr(data_list, "sample_metadata"))) {
      sample_metadata <- attr(data_list, "sample_metadata")
      group_names <- unique(sample_metadata$group_id)
      if (diagnostics) {
        message("Using group names extracted from metadata: ",
                paste(group_names, collapse = ", "))
      }
    } else if (is.list(data_list) && all(names(data_list) != "")) {
      # If no metadata but data_list has names, use those as group names
      group_names <- names(data_list)
      if (diagnostics) {
        message("Using list names as group names: ",
                paste(group_names, collapse = ", "))
      }
    } else {
      stop("No group_names provided and unable to extract from metadata")
    }
  }

  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (!all(sapply(data_list, is.list))) stop("data_list must be a list of lists")

  result <- data_list

  # Step 1: Within-group normalization
  if(within_groups) {
    if(diagnostics) cat("\n=== Performing within-group normalization ===\n")
    result <- lapply(names(result), function(group) {
      if(diagnostics) cat(sprintf("\nProcessing group: %s\n", group))
      .normalize_methylation_within_set(
        result[[group]],
        alpha = within_alpha,
        sites_per_quantile = sites_per_quantile,
        max_quantiles = max_quantiles,
        diagnostics = diagnostics
      )
    })
    names(result) <- names(data_list)
  }

  # Step 2: Between-group normalization
  if(between_groups) {
    if(diagnostics) cat("\n=== Performing between-group normalization ===\n")
    # Flatten data for between-group normalization
    all_samples <- unlist(result, recursive = FALSE)

    # Perform gentler between-group normalization
    normalized_data <- .normalize_methylation_within_set(
      all_samples,
      alpha = between_alpha,  # Higher alpha to preserve more group differences
      sites_per_quantile = sites_per_quantile * 2,  # Fewer quantiles for between-group
      max_quantiles = max_quantiles,
      diagnostics = diagnostics
    )

    # Reconstruct the group structure
    result <- list()
    current_idx <- 1
    for(group in names(data_list)) {
      n_samples <- length(data_list[[group]])
      result[[group]] <- normalized_data[current_idx:(current_idx + n_samples - 1)]
      current_idx <- current_idx + n_samples
    }
  }

  names(result) <- names(data_list)
  return(result)
}
