#' Load and Preprocess Methylation Data
#'
#' Load methylation data from files based on a sample sheet.
#'
#' @param dir_path Path to the directory containing the data files.
#' @param sample_sheet Path to the sample sheet file or a data.frame/data.table containing
#'   sample information with columns: group_id, replicate, sample_id, file_name.
#' @param type Character string specifying which type of files to load.
#'   If NULL (default), it will use file names from the sample sheet as is.
#' @param groups Character vector of group_ids to include. If NULL (default),
#'   all groups will be loaded.
#' @param cores Number of cores to use for parallel processing. Default is 1
#'   (no parallel processing).
#' @param use_cpp Logical, whether to use the C++ implementation for faster loading.
#'   Default is TRUE. Will fall back to R implementation if C++ is not available.
#' @param single_file Optional path to a single methylation file to load directly
#'   When provided, the sample_sheet parameter is ignored.
#' @return A list of data.tables containing processed methylation data, or a
#'    single data.table if single_file is provided.
#'
#' @export
#' @useDynLib MAPitNorm, .registration = TRUE
load_data <- function(dir_path, sample_sheet, type = NULL, groups = NULL,
                      cores = 1, use_cpp = TRUE, single_file = NULL) {

  # Handle single file case
  if (!is.null(single_file)) {
    if (!missing(sample_sheet)) {
      message("Ignoring sample_sheet parameter since single_file is provided")
    }

    # If single_file is a path relative to dir_path, construct the full path
    if (missing(dir_path) || is.null(dir_path)) {
      file_path <- single_file
    } else {
      file_path <- file.path(dir_path, single_file)
    }

    return(load_single_file(file_path, use_cpp = use_cpp))
  }

  # If we reach here, we're loading multiple files

  if (is.null(sample_sheet)) {
    stop("sample_sheet must be provided when not using single_file")
  }

  # Load sample sheet
  sample_data <- .load_sample_sheet(sample_sheet)

  # Filter by groups if specified
  if (!is.null(groups)) {
    sample_data <- .filter_by_groups(sample_data, groups)
  }

  # Filter by type if specified
  if (!is.null(type)) {
    sample_data <- .filter_by_type(sample_data, type)
  }

  # Construct file paths
  files.list <- file.path(dir_path, sample_data$file_name)
  names(files.list) <- sample_data$sample_id

  # Check file existence
  files.list <- .check_files_exist(files.list)
  if (length(files.list) < 2) {
    stop("At least two valid files are required for analysis.")
  }

  # Determine whether to use C++ implementation
  has_rcpp <- requireNamespace("Rcpp", quietly = TRUE) &&
    exists("readMethylationFiles", mode = "function")

  use_cpp_impl <- use_cpp && has_rcpp

  if (use_cpp && !has_rcpp) {
    message("Rcpp package not available. Using R implementation.")
    message("For faster processing, install the Rcpp package.")
    use_cpp_impl <- FALSE
  }

  message(sprintf("Loading %d files using %s implementation...",
                  length(files.list),
                  ifelse(use_cpp_impl, "C++", "R")))

  start_time <- Sys.time()

  # Load data using appropriate implementation
  if (use_cpp_impl) {
    all_samples <- .load_data_cpp(files.list, cores)
  } else {
    all_samples <- .load_data_r(files.list, cores)
  }

  # Post-process and attach metadata
  all_samples <- .post_process_samples(all_samples, sample_data, files.list)

  end_time <- Sys.time()
  elapsed <- end_time - start_time
  message(sprintf("Total loading time: %.2f seconds", as.numeric(elapsed, units="secs")))

  return(all_samples)
}


