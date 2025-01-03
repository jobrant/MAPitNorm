# File: R/load_data.R

#' Load and Preprocess Data
#'
#' Load methylation data from files, preprocess, and split into GCH and
#' HCG samples. Files are expected to be named groupName_replicateName_site.tsv
#' Site is either GCH or HCG.
#'
#' @param dir_path Path to the directory containing the data files.
#'
#' @return A list containing GCH and HCG samples as data frames.
#'
#' @importFrom purrr map
#' @importFrom stringr str_sub
#' @importFrom progress progress_bar
#' @importFrom data.table fread
#'
#' @export
#'

load_data <- function(dir_path) {
  # List all TSV files in the directory
  files.list <- list.files(path = dir_path, pattern = "tsv", full.names = T)

  # Check if files.list is empty
  if(length(files.list) == 0) {
    stop("No methylation call files found in the directory.")
  }

  # Filter based on file type
  gch_files <- files.list[grepl("GCH", files.list)]
  hcg_files <- files.list[grepl("HCG", files.list)]

  # Count the number of GCH and HCG files
  num_gch_files <- length(gch_files)
  num_hcg_files <- length(hcg_files)

  # Print the file count information
  message(paste("Number of GCH files: ", num_gch_files))
  message(paste("Number of HCG files: ", num_hcg_files))

  # Check if at least two GCH files are present
  if(num_gch_files < 2) {
    stop("At least two GCH files are required for analysis.")
  }

  files.names <- basename(str_sub(string = files.list, end = -8))
  names(files.list) <- files.names

  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  Loading [:bar] :percent in :elapsed",
    total = length(files.list), clear = FALSE
  )

  # Load and preprocess each file, updating the progress bar
  all_samples <- map(files.list, function(file) {
    pb$tick()  # Update progress bar
    df <- fread(file)
    colnames(df) <- c("chr", "pos", "strand", "site", "mc", "cov", "stat")
    df <- df[, .("chr", "pos", "strand", "site", "mc", "cov")]
    df$rate <- df$mc / df$cov
    df$uniqueID <- paste(df$chr, df$pos, df$site, sep = "_")
    return(df)
  })

  gch_samples <- all_samples[grepl(pattern = "GCH", x = names(all_samples))]
  hcg_samples <- all_samples[grepl(pattern = "HCG", x = names(all_samples))]

  list(GCH = gch_samples, HCG = hcg_samples)
}
