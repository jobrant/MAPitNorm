# File: R/find_shared_sites.R

#' Find Shared Sites
#'
#' Identify shared sites across all samples within GCH or HCG data.
#'
#' @param samples List of sample data frames.
#'
#'
#' @return A vector of shared site identifiers.
#'
#' @export
#'
find_shared_sites <- function(samples) {
  Reduce(intersect, lapply(samples, function(x) x$uniqueID))
}
