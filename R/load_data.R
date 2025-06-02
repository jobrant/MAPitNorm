#' Load and Preprocess Data
#'
#' Load methylation data from files and preprocess them.
#' Files are expected to be named groupName_replicateName_site.tsv
#' Site is either GCH or HCG.
#'
#' @param dir_path Path to the directory containing the data files.
#' @param type Character string specifying which type of files to load.
#' @param cores Number of cores to use for parallel processing. Default is 1
#'   (no parallel processing).
#'
#' @return A list of data.tables containing processed methylation data.
#'
#' @importFrom stringr str_sub
#' @importFrom progress progress_bar
#' @importFrom data.table fread setnames :=
#' @importFrom parallel makeCluster stopCluster parLapply detectCores clusterExport
#'
#' @examples
#' \dontrun{
#' # Load GCH files sequentially
#' gch_samples <- load_data("path/to/files", type = "GCH")
#'
#' # Load HCG files using parallel processing
#' hcg_samples <- load_data("path/to/files", type = "HCG", cores = 4)
#' }
#'
#' @export
#'

load_data <- function(dir_path, sample_key_path, type = c("GCH", "HCG"), cores = 1) {
  # Match arg
  type <- match.arg(type)

  # Validate sample key path
  if (missing(sample_key_path) || !file.exists(sample_key_path)) {
    stop("Sample key file path must be provided and exist.")
  }

  # Read sample key (assuming tab-separated with headers)
  sample_key <- fread(sample_key_path)
  required_cols <- c("sample", "Time_point", "file_path")
  if (!all(required_cols %in% colnames(sample_key))) {
    stop("Sample key must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Validate file paths and filter for methylation type
  sample_key[, file_path_full := {
    # Check if file_path is absolute or relative
    is_absolute <- startsWith(file_path, "/") | startsWith(file_path, "\\") | grepl("^[A-Za-z]:", file_path)
    full_path <- ifelse(is_absolute, file_path, file.path(dir_path, file_path))
    # Normalize path separators
    normalizePath(full_path, mustWork = FALSE)
  }]

  # Filter for files containing the methylation type and verify existence
  sample_key <- sample_key[str_detect(basename(file_path_full), type) & file.exists(file_path_full)]
  if (nrow(sample_key) == 0) {
    stop(sprintf("No %s files found matching sample key entries.", type))
  }
  if (nrow(sample_key) < 2) {
    stop(sprintf("At least two %s files with valid paths are required.", type))
  }

  # Warn about non-existent or non-matching files
  invalid_files <- sample_key[!file.exists(file_path_full) | !str_detect(basename(file_path_full), type)]
  if (nrow(invalid_files) > 0) {
    warning("Some sample key entries have invalid or non-matching file paths: ",
            paste(invalid_files$sample, invalid_files$Time_point, sep = "_", collapse = ", "))
  }

  # Define columns needed for processing
  cols_needed <- c("chr", "pos", "strand", "site", "mc", "cov")

  # Load files
  files_to_load <- sample_key$file_path_full
  if (cores > 1) {
    cl <- makeCluster(min(cores, detectCores() - 1))
    on.exit(stopCluster(cl))
    clusterExport(cl, varlist = "cols_needed", envir = environment())
    clusterEvalQ(cl, { library(data.table) })
    all_data <- parLapply(cl, files_to_load, function(file) {
      df <- fread(file, select = 1:6, showProgress = FALSE)
      setnames(df, cols_needed)
      df[, `:=`(rate = mc / cov, uniqueID = paste(chr, pos, site, sep = "_"))]
      return(df)
    })
  } else {
    pb <- progress_bar$new(
      format = "  Loading [:bar] :percent in :elapsed",
      total = length(files_to_load),
      clear = FALSE
    )
    all_data <- lapply(files_to_load, function(file) {
      pb$tick()
      df <- fread(file, select = 1:6, showProgress = FALSE)
      setnames(df, cols_needed)
      df[, `:=`(rate = mc / cov, uniqueID = paste(chr, pos, site, sep = "_"))]
      return(df)
    })
  }

  # Set names for the data list
  sample_names <- paste(sample_key$sample, sample_key$Time_point, sep = "_")
  names(all_data) <- sample_names

  # Prepare metadata
  metadata <- copy(sample_key)
  metadata[, sample_name := sample_names]
  metadata[, file_path_full := NULL] # Remove temporary column

  # Return data and metadata
  return(list(data = all_data, metadata = metadata))
}
