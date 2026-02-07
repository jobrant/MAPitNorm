#' Convert SummarizedExperiment back to a nested list of data tables
#'
#' @description
#' This function performs the inverse operation of \code{create_RSE}. It extracts
#' genomic coordinates and assay data, then reconstructs a hierarchical list
#' based on the groups provided in the \code{Group_SampleName} format.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param groups Character vector. The names of specific groups to extract.
#' If \code{NULL} (default), all samples are extracted.
#'
#' @return A doubly-nested named list of \code{data.table} objects.
#'
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import data.table
#' @importFrom stringr str_split
#' @export

convertSE_to_list <- function(rse, groups = NULL, verbose = TRUE) {

  # Validate input object
  if (!is(rse, "SummarizedExperiment") && !is(rse, "RangedSummarizedExperiment")) {
    stop("Input 'rse' must be a SummarizedExperiment or RangedSummarizedExperiment object.")
  }

  # Step 1: Parse Group and Sample names from Row Names
  # We use rownames(colData) which are the experiment's column names
  all_col_ids <- rownames(SummarizedExperiment::colData(rse))

  # Parse the "Sample_Sample" format
  name_split <- stringr::str_split(all_col_ids, pattern = "_", n = 2, simplify = TRUE)
  parsed_samples <- name_split[, 1]
  parsed_reps    <- name_split[, 2]

  # Handle Group Filtering
  if (!is.null(groups)) {
    if (verbose) message("Filtering for specified groups...")
    keep_idx <- which(parsed_samples %in% groups)

    if (length(keep_idx) == 0) {
      stop("None of the specified groups were found in the column IDs.")
    }

    # Subset SE or RSE and our parsed vectors
    rse <- rse[, keep_idx]

    # CRITICAL: Drop unused factor levels in colData (e.g., removing Group B labels)
    SummarizedExperiment::colData(rse) <- droplevels(SummarizedExperiment::colData(rse))

    parsed_samples <- parsed_samples[keep_idx]
    parsed_reps    <- parsed_reps[keep_idx]
    all_col_ids    <- all_col_ids[keep_idx]
  }



  # Step 2: Extracting genomic metadata from rowRanges..."
  if (verbose) message("Extracting genomic metadata from rowRanges...")
  # Extract the 'spine' (chr, pos, strand, etc.)
  # mcols() extracts the metadata columns like 'site' and 'uniqueID'
  coords <- as.data.table(SummarizedExperiment::rowRanges(rse))

  # Ensure column names match your original requirements
  # GRanges 'seqnames' becomes 'chr', 'start' becomes 'pos'
  setnames(coords, old = c("seqnames", "start"), new = c("chr", "pos"))

  # Drop IRanges specific columns like 'end' and 'width' if the width is 1 throughout
  if(length(unique(coords$width))==1){
    coords <- coords[, .(chr, pos, strand, site, uniqueID)]
  }

  # Extract assay data
  if (verbose) message("Extracting assay data...")
  cov_arr <- SummarizedExperiment::assay(rse, "coverage")
  rate_arr <- SummarizedExperiment::assay(rse, "rate")

  # Create a flat list of data.tables
  flat_list <- lapply(seq_len(ncol(rse)), function(i) {
    # Combine coordinates with measurements for this specific column
    dt <- cbind(coords,
                cov = cov_arr[, i],
                rate = rate_arr[, i])

    # Optional: Calculate mc (methylated counts) if you need the original full set
    dt[, mc := cov * rate]

    # Return as a data.table
    return(data.table::as.data.table(dt))
  })

  names(flat_list) <- colnames(rse)

  if (verbose) message ("Nesting into Group/Sample hierarchy...")

  # Split by Group (Tier 1)
  unique_samples <- unique(parsed_samples)

  lol <- lapply(unique_samples, function(s_name) {
    sample_idx <- which(parsed_samples == s_name)
    sample_reps <- flat_list[sample_idx]
    names(sample_reps) <- parsed_reps[sample_idx]
    return(sample_reps)
  })

  names(lol) <- unique_samples

  if (verbose) message("Conversion complete")

  return(lol)
}
