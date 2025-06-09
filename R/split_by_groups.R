#' Split data into Groups and Replicates
#'
#' Splits data list into group specific lists so all replicates per group are together.
#' Works with data loaded using load_data().
#'
#' @param data_list List containing methylation data as returned by load_data()..
#'
#' @param group_names Optional vector of group names to include. If NULL (default),
#'   all groups found in the sample metadata will be used.
#'
#' @return A list containing samples organized by group,
#' with each group containing its replicates.
#'
#' @export
#'

split_by_groups <- function(data_list, group_names = NULL) {
  # Check inputs
  if (!is.list(data_list)) {
    stop("data_list must be a list")
  }

  # Get sample metadata
  sample_metadata <- attr(data_list, "sample_metadata")
  if (is.null(sample_metadata)) {
    stop("Input doesn't have sample metadata. Make sure it's from load_data().")
  }

  # Extract available groups
  available_groups <- unique(sample_metadata$group_id)

  # If group not provided, use all available groups
  if (is.null(group_names)) {
    group_names <- available_groups
    message("Using all available groups: ", paste(group_names, collapse = ","))
  } else {
    # Validate requested groups
    if (!is.character(group_names)) {
      stop("group_names must be a character vector")
    }

    # Check if all requested groups exist in the data
    missing_groups <- setdiff(group_names, available_groups)
    if (length(missing_groups) > 0) {
      warning("The following groups are not found in the data: ",
              paste(missing_groups, collapse = ", "))
    }

    # Filter to only include existing groups
    group_names <- intersect(group_names, available_groups)
    if (length(group_names) == 0) {
      stop("None of the requested groups were found in the data.")
    }
  }

  result_list <- list()

  # For each group, collect all its samples
  for (group in group_names) {
    # Get sample IDs for this group
    group_sample_ids <- sample_metadata$sample_id[sample_metadata$group_id == group]

    #Extract corresponding samples from the data_list
    group_samples <- data_list[group_sample_ids]

    if (length(group_samples) == 0) {
      warning("No Samples found for group: ", group)
    } else {
      result_list[[group]] <- group_samples
    }
  }

  # Preserve sample metadata at the top level
  attr(result_list, "sample_metadata") <- sample_metadata

  return(result_list)

}

