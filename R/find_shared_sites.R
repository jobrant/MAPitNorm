#' Find and Filter Shared Sites
#'
#' Identify shared sites across all samples and optionally filters data to only include these sites.
#'
#' @param allc A list of samples returned by load_data(), with each element
#' being a data.table containing a `uniqueID` column.
#' @param filter Logical, whether to filter the data to only include shared sites (default = TRUE).
#'  If FALSE, only returns the vector of shared site IDs.
#' @param quiet Logical, whether to suppress progress messages (default = FALSE).
#'
#' @return If filter=TRUE, returns the filtered list with only shared sites.
#'  If filter=FALSE, returns a vector of shared site identifiers (uniqueID values).
#'
#' @export
find_shared_sites <- function(allc, filter = TRUE, quiet = FALSE) {
  # Input validation
  if (!is.list(allc) || length(allc) == 0) {
    stop("Input must be a non-empty list of samples, as returned by load_data().")
  }

  # Check if the first element is a data.table with uniqueID column
  if (!("data.table" %in% class(allc[[1]])) || !("uniqueID" %in% colnames(allc[[1]]))) {
    stop("Each element in the list must be a data.table with a 'uniqueID' column.")
  }

  # Create an index for faster filtering later
  if (filter && requireNamespace("data.table", quietly = TRUE)) {
    # Make sure all samples have uniqueID as key for faster operations
    for (i in seq_along(allc)) {
      if (!data.table::haskey(allc[[i]]) ||
          !identical(data.table::key(allc[[i]]), "uniqueID")) {
        data.table::setkeyv(allc[[i]], "uniqueID")
        if (!quiet) message("Setting key on sample ", i, " for faster filtering")
      }
    }
  }

  if (!quiet) message("Finding shared sites across ", length(allc), " samples...")

  if (length(allc) > 20) {
    # For many samples, more efficient to use a frequency table approach
    # Extract all uniqueIDs
    all_ids <- unlist(lapply(allc, function(x) x$uniqueID))

    # Count frequency of each ID
    id_table <- table(all_ids)

    # Find IDs that appear in all samples
    shared_sites <- names(id_table)[id_table == length(allc)]

  } else {
    # For fewer samples, Reduce approach is fine
    shared_sites <- Reduce(intersect, lapply(allc, function(x) x$uniqueID))
  }

  if (!quiet) message("Found ", length(shared_sites), " shared sites.")

  # If filter=FALSE, return only the shared site IDs
  if (!filter) {
    return(shared_sites)
  }

  # Filter each sample to include only shared sites
  if (!quiet) {
    message("Filtering samples to include only shared sites...")
    pb <- utils::txtProgressBar(min = 0, max = length(allc), style = 3)
  }

  # Fast filtering using data.table
  filtered_data <- vector("list", length(allc))
  names(filtered_data) <- names(allc)

  for (i in seq_along(allc)) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      # Fast filtering with data.table
      filtered_data[[i]] <- allc[[i]][shared_sites]
    } else {
      # Fallback to base R subsetting
      filtered_data[[i]] <- allc[[i]][allc[[i]]$uniqueID %in% shared_sites, ]
    }

    if (!quiet) utils::setTxtProgressBar(pb, i)
  }

  if (!quiet) {
    close(pb)
    message("Filtering complete.")
  }

  # Preserve attributes
  attributes(filtered_data) <- attributes(allc)

  return(filtered_data)
}

