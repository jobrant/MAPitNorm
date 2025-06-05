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
#'
#' @return A list of data.tables containing processed methylation data.
#'
#' @importFrom data.table fread setnames := is.data.table as.data.table rbindlist
#' @importFrom methods is
#'
#' @export
load_data <- function(dir_path, sample_sheet, type = NULL, groups = NULL,
                      cores = 1, use_cpp = TRUE) {

  # Load sample sheet
  sample_data <- load_sample_sheet(sample_sheet)

  # Filter by groups if specified
  if (!is.null(groups)) {
    sample_data <- filter_by_groups(sample_data, groups)
  }

  # Filter by type if specified
  if (!is.null(type)) {
    sample_data <- filter_by_type(sample_data, type)
  }

  # Construct file paths
  files.list <- file.path(dir_path, sample_data$file_name)
  names(files.list) <- sample_data$sample_id

  # Check file existence
  files.list <- check_files_exist(files.list)
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
    all_samples <- load_data_cpp(files.list, cores)
  } else {
    all_samples <- load_data_r(files.list, cores)
  }

  # Post-process and attach metadata
  all_samples <- post_process_samples(all_samples, sample_data, files.list)

  end_time <- Sys.time()
  elapsed <- end_time - start_time
  message(sprintf("Total loading time: %.2f seconds", as.numeric(elapsed, units="secs")))

  return(all_samples)
}

# Internal function to load using C++ implementation
load_data_cpp <- function(files.list, cores) {
  # Call the C++ function
  raw_data <- readMethylationFiles(files.list)
  names(raw_data) <- names(files.list)

  # Convert to data.tables and add calculated columns
  all_samples <- lapply(seq_along(raw_data), function(i) {
    df <- as.data.table(raw_data[[i]])
    df[, rate := ifelse(cov > 0, mc/cov, 0)]
    df[, uniqueID := paste(chr, pos, site, sep="_")]
    return(df)
  })

  names(all_samples) <- names(raw_data)
  return(all_samples)
}

# Internal function as fallback using pure R
load_data_r <- function(files.list, cores) {
  # Import helper functions using ::
  # This ensures the package works even if these packages aren't loaded
  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    # Use parallel processing
    all_samples <- load_data_parallel(files.list, cores)
  } else {
    # Use sequential processing with progress bar
    all_samples <- load_data_sequential(files.list)
  }

  return(all_samples)
}

