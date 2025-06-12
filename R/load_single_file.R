#' Load a single methylation data file
#'
#' Load a single methylation data file for visualization or other purposes.
#'
#' @param file_path Path to the methylation data file.
#' @param use_cpp Logical, whether to use the C++ implementation for faster loading.
#'   Default is TRUE. Will fall back to R implementation if C++ is not available.
#' @return A data.table containing the methylation data.
#' @export
load_single_file <- function(file_path, use_cpp = TRUE) {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  # Check if C++ implementation is available
  has_cpp <- use_cpp && requireNamespace("Rcpp", quietly = TRUE)

  message(sprintf("Loading file using %s implementation: %s",
                  ifelse(has_cpp, "C++", "R"),
                  basename(file_path)))

  start_time <- Sys.time()

  # Load data using appropriate implementation
  if (has_cpp) {
    # Use the C++ function
    df <- readMethylationFile(file_path)
    df <- data.table::as.data.table(df)
    df[, rate := ifelse(cov > 0, mc/cov, 0)]
    df[, uniqueID := paste(chr, pos, site, sep="_")]
  } else {
    df <- .load_single_file_r(file_path)
  }

  end_time <- Sys.time()
  elapsed <- end_time - start_time
  message(sprintf("Loading time: %.2f seconds", as.numeric(elapsed, units="secs")))

  return(df)
}
