# File: R/split_by_groups.R

#' Split data into Groups and Replicates
#'
#' Splits data list into group specific lists so all replicates per group are together
#'
#' @param data_list List containing GCH or HCG methylation data.
#'
#' @param group_names Vector of group names to split by.
#'
#' @return A list containing GCH and HCG samples with replicates by group.
#'
#' @export
#'

split_by_groups <- function(data_list, group_names) {
  # Check inputs
  if (!is.list(data_list)) stop("data_list must be a list")
  if (is.null(names(data_list))) stop("data_list must have names")
  if (!is.character(group_names)) stop("group_names must be a character vector")

  # Create a pattern that matches any of the group names
  name_groups <- sapply(names(data_list), function(x) {
    # Remove 'allc_' prefix first
    x_clean <- sub("^allc_", "", x)
    # Then check for group matches
    matches <- group_names[sapply(group_names, function(g) startsWith(x_clean, g))]
    if(length(matches) == 0) {
      stop("Sample name '", x, "' doesn't match any group names")
    } else {
      matches[which.max(nchar(matches))]
    }
  })

  # Check if any samples were found for each group
  for(group in group_names) {
    if(sum(name_groups == group) == 0) {
      warning("No samples found for group: ", group)
    }
  }

  # Split the list
  split(data_list, factor(name_groups, levels = group_names))
}
