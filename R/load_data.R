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

load_data <- function(dir_path, type = c("GCH", "HCG"), cores = 1) {
  # Match arg
  type <- match.arg(type)

  # List only files of requested type
  files.list <- list.files(
    path = dir_path,
    pattern = paste0(".*", type, ".*\\.tsv"),
    full.names = T
    )

  # Validation step
  if(length(files.list) == 0) {
    stop(sprintf("No %s methylation call files found in the directory.", type))
  }
  if(length(files.list) < 2) {
    stop(sprintf("At least two %s files are required for analysis.", type))
  }

  # Name files
  file.names <- basename(tools::file_path_sans_ext(files.list, compression = T))
  names(files.list) <- file.names

  message(sprintf("Loading %d %s files...", length(files.list), type))

  # Define columns we need
  cols_needed <- c("chr", "pos", "strand", "site", "mc", "cov")

  message(sprintf("Starting processing with %d cores", cores))
  start_time <- Sys.time()

  # Set up parallel or sequential processing
  if(cores > 1) {
    # Parallel processing
    tryCatch({
      message("Cluster created at: ", Sys.time())
      cl <- parallel::makeCluster(min(cores, parallel::detectCores() -1))
      on.exit(parallel::stopCluster(cl), add = T)

      parallel::clusterExport(cl = cl, varlist = "cols_needed", envir = environment())
      parallel::clusterEvalQ(cl, {
        library(data.table)
        NULL
      })

      all_samples <- parallel::parLapply(cl, files.list, function(file) {
        tryCatch({
          df <- fread(
            file,
            select = 1:6,
            showProgress = FALSE
          )
          setnames(df, cols_needed)
          df[, ":="(
            rate = mc/cov,
            uniqueID = paste(chr, pos, site, sep="_")
          )]
          return(df)
        }, error = function(e) {
          message("Error processing file: ", file, "\n", e$message)
          return(NULL)
        })
      })

    }, error = function(e) {
      message("Error in parallel processing: ", e$message)
      message("Falling back to sequential processing...")
      return(NULL)
    })

    # If parallel processing failed, fall back to sequential
    if(is.null(all_samples)) {
      cores <- 1
    }
  }

  if(cores == 1) {
    pb <- progress::progress_bar$new(
      format = "  Loading [:bar] :percent in :elapsed",
      total = length(files.list),
      clear = FALSE
    )

    all_samples <- lapply(files.list, function(file) {
      pb$tick()
      df <- data.table::fread(
        file,
        select = 1:6,
        showProgress = FALSE
      )
      data.table::setnames(df, cols_needed)
      df[, ":="(
        rate = mc/cov,
        uniqueID = paste(chr, pos, site, sep="_")
      )]
      return(df)
    })
  }

  # Remove any NULL entries from failed processing
  all_samples <- all_samples[!sapply(all_samples, is.null)]
  names(all_samples) <- file.names[!sapply(all_samples, is.null)]

  gc(verbose = FALSE)
  return(all_samples)
}
