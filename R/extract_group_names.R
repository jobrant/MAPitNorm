#' Extract group names from data
#' @param data_list A list of methylation data
#' @param diagnostics Whether to print diagnostic messages
#' @return Character vector of group names
#' @keywords internal
extract_group_names <- function(data_list, diagnostics = FALSE) {
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
    stop("Unable to extract group names from data")
  }
  return(group_names)
}


## To use add this bit to appropriate functions
normalize_coverage <- function(data_list, group_names = NULL, between_groups = FALSE,
                               diagnostics = FALSE) {
  # Extract group names if not provided
  if (is.null(group_names)) {
    group_names <- extract_group_names(data_list, diagnostics)
  }

  # Rest of function...
}
