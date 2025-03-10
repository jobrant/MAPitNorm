#' Sample methylation data for visualization
#'
#' @param original_data Original data list
#' @param within_norm Within-group normalized data
#' @param full_norm Fully normalized data
#' @param n_sites Number of sites to sample (default = 100000)
#' @param seed Random seed for reproducibility (default = 42)
#'
#' @return List containing sampled versions of input data maintaining the same structure
#' @export
sample_for_visualization <- function(original_data, within_norm, full_norm,
                                     n_sites = 100000, seed = 42) {
  # Set seed for reproducibility
  set.seed(seed)

  # Get total number of sites from first sample of first group
  total_sites <- nrow(original_data[[1]][[1]])

  # Sample indices
  sample_idx <- sort(sample(total_sites, min(n_sites, total_sites)))

  # Function to subset data list
  subset_data <- function(data_list) {
    lapply(data_list, function(group) {
      lapply(group, function(df) {
        df[sample_idx, ]
      })
    })
  }

  # Create sampled versions of each dataset
  list(
    original = subset_data(original_data),
    within_norm = subset_data(within_norm),
    full_norm = subset_data(full_norm)
  )
}
