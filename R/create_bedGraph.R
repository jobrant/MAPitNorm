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
#'
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
                            sample_name = NULL) {

  # Match and check type argument
  type <- match.arg(type)

  # Check and create output directory if needed
  out <- normalizePath(out, mustWork = FALSE)
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE)
  }

  # Handle different input types
  if (is.character(data)) {
    # Input is a file path
    if (!file.exists(data)) {
      stop("Input file does not exist: ", data)
    }
    sample_name <- if(is.null(sample_name)) basename(data) else sample_name
    return(process_single_file(data, type, out, sample_name))

  } else if (is.data.frame(data)) {
    # Input is a single sample
    return(process_single_sample(sample_df = data,
                                 type = type,
                                 out = out,
                                 sample_name = sample_name))

  } else if (is.list(data)) {
    # Check if it's a nested list (list of lists)
    if (all) {
      if (!all(sapply(data, is.list))) {
        # Single group from a list of lists
        group_name <- deparse(substitute(data))  # Get the variable name
        group_name <- sub(".*\\$", "", group_name)  # Extract just the group name
        return(process_single_sample(data, type, out,
                                     group_name = group_name,
                                     sample_name = sample_name))
      }
      # Full list of lists
      return(process_all_samples(data, type, out))
    } else {
      stop("When 'all=FALSE', input should be a single sample or file path")
    }
  } else {
    stop("Invalid input type. Must be either a file path, data frame, or list of samples")
  }
}


#' Process a single methylation file into bedGraph format
#' @param file_path Path to methylation file
#' @param type Type of bedGraph to create ("meth" or "cov")
#' @param out Output directory path
#' @return Path to created bedGraph file
#' @keywords internal
process_single_file <- function(file_path, type, out) {
  # Load and process file
  data <- read_methylation_file(file_path)  # You'll need to implement this
  process_single_sample(data, type, out)
}


#' Process a single sample data frame into bedGraph format
#' @param sample_df Data frame containing methylation data
#' @param type Type of bedGraph to create ("meth" or "cov")
#' @param out Output directory path
#' @param sample_name Name of the sample (optional)
#'
#' @importFrom data.table fwrite
#' @return Path to created bedGraph file
#' @keywords internal
process_single_sample <- function(sample_df, type, out, group_name = NULL, sample_name = NULL) {
  # Convert to bedGraph format
  bg_data <- data.table::data.table(
    chr = paste0("chr", sample_df$chr),
    start = sample_df$pos - 1,  # BED format is 0-based
    end = sample_df$pos,
    value = if(type == "meth") sample_df$rate else sample_df$cov
  )

  # Create output filename
  if(is.null(sample_name)) {
    file_prefix <- "sample"
    sample_name <- basename(tempfile())
  } else {
    file_prefix <- if(is.null(group_name)) sample_name else paste(group_name, sample_name, sep="_")
  }

  out_file <- file.path(out, paste0(file_prefix, "_", type, ".bedGraph"))

  data.table::fwrite(bg_data, out_file,
              sep = "\t", quote = FALSE, col.names = FALSE)

  return(out_file)
}


#' Process all samples in a data list into bedGraph format
#' @param data_list List of lists containing methylation data
#' @param type Type of bedGraph to create ("rate" or "cov")
#' @param out Output directory path
#' @return Vector of paths to created bedGraph files
#' @keywords internal
process_all_samples <- function(data_list, type, out) {
  out_files <- lapply(names(data_list), function(group) {
    lapply(names(data_list[[group]]), function(sample) {
      process_single_sample(
        sample_df = data_list[[group]][[sample]],
        type = type,
        out = out,
        group_name = group,
        sample_name = sample
      )
    })
  })

  return(unlist(out_files))
}



#' Reads single methylation data file
#' @param file_path Path to methylation file
#'
#' @return data.table containing methylation data
#'
#' @importFrom progress progress_bar
#' @importFrom data.table fread setnames :=
#'
#' @keywords internal
read_methylation_file <- function(file_path) {
  cols_needed <- c("chr", "pos", "strand", "site", "mc", "cov")
  pb <- progress::progress_bar$new(
    format = "  Loading [:bar] :percent in :elapsed",
    total = length(1),
    clear = FALSE
  )
  pb$tick()
  df <- data.table::fread(
    file_path,
    select = 1:6,
    showProgress = FALSE
  )
  data.table::setnames(df, cols_needed)
  df[, ":="(
    rate = mc/cov,
    uniqueID = paste(chr, pos, site, sep="_")
  )]
  return(df)
}


