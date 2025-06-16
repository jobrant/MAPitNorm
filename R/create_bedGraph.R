#' Converts normalized data into bedGraph files
#'
#' @param data Input data, can be either:
#'        1. A list of lists containing normalized samples
#'        2. A single data frame containing one sample
#'        3. A character string specifying path to input file
#' @param type Character string for type of value to include in bedGraph.
#'        Should be one of "cov" or "meth" (for coverage bedGraph or methylation rate bedGraph, respectively.)
#' @param all Logical. Should function be performed on all files in object. Default is TRUE.
#'        Ignored if input is a single sample or file.
#' @param out Directory to write output bedGraph files to. Default is current working directory.
#' @param group_name Name of the group (optional)
#' @param sample_name Name to use for output file (optional)
#' @param use_cpp Logical. Should function be performed on all files in object.
#'        Default is TRUE. Will fall back to R implementation if C++ is not available.
#' @details
#' Converts methylation data to bedGraph format. Can handle multiple input types:
#' - If given a normalized data object (list of lists), processes all samples by default
#' - If given a single sample data frame, processes just that sample
#' - If given a file path, loads and processes that file
#'
#' @return Invisibly returns list of created bedGraph file paths
#' @export
create_bedgraph <- function(data,
                            type = c("meth", "cov"),
                            all = TRUE,
                            out = ".",
                            group_name = NULL,
                            sample_name = NULL,
                            use_cpp = TRUE) {

  # Match and check type argument
  type <- match.arg(type)

  # Check and create output directory if needed
  out <- normalizePath(out, mustWork = FALSE)
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE)
  }

  # Determine C++ implementation availability
  has_rcpp <- requireNamespace("Rcpp", quietly = TRUE) &&
    is.function(get("createBedgraphCpp", envir = asNamespace("MAPitNorm")))

  use_cpp_impl <- use_cpp && has_rcpp

  if (use_cpp && !has_rcpp) {
    message("Rcpp or createBedgraphCpp not available. Using R implementation.")
    message("For faster processing, ensure Rcpp is installed and createBedgraphCpp is exported.")
    use_cpp_impl <- FALSE
  }

  # Store implementation choice
  .options <- list(use_cpp_impl = use_cpp_impl)

    # Handle different input types
  if (is.character(data)) {
    # Input is a file path
    if (!file.exists(data)) {
      stop("Input file does not exist: ", data)
    }
    sample_name <- if(is.null(sample_name)) basename(data) else sample_name
    return(process_single_file(data, type, out, sample_name, .options))

  } else if (is.data.frame(data)) {
    # Input is a single sample
    return(process_single_sample(sample_df = data,
                                 type = type,
                                 out = out,
                                 sample_name = sample_name,
                                 .options = .options))

  } else if (is.list(data)) {
    # Check if it's a nested list (list of lists)
    if (all) {
      if (!all(sapply(data, is.list))) {
        # Single group from a list of lists
        group_name <- deparse(substitute(data))  # Get the variable name
        group_name <- sub(".*\\$", "", group_name)  # Extract just the group name
        return(process_single_sample(data, type, out,
                                     group_name = group_name,
                                     sample_name = sample_name,
                                     .options = .options))
      }
      # Full list of lists
      return(process_all_samples(data, type, out, .options))
    } else {
      stop("When 'all=FALSE', input should be a single sample or file path")
    }
  } else {
    stop("Invalid input type. Must be either a file path, data frame, or list of samples")
  }
}
