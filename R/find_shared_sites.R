#' Find Shared Sites
#'
#' Identify shared sites across all samples within GCH or HCG data loaded using load_data().
#'
#' @param allc A list returned by load_data(), containing a `data` element with a list
#'   of sample data.tables, each with a `uniqueID` column.
#'
#' @return A vector of shared site identifiers (uniqueID values).
#'
#' @export
find_shared_sites <- function(allc) {
  # Input validation
  if (!is.list(allc) || !"data" %in% names(allc)) {
    stop("Input must be a list with a 'data' element, as returned by load_data().")
  }
  if (length(allc$data) == 0) {
    stop("No samples found in allc$data.")
  }

  # Extract uniqueID from each sample's data.table
  Reduce(intersect, lapply(allc$data, function(x) x$uniqueID))
}
