#' @param use_cpp Logical, whether to use the C++ implementation for faster loading.
#'   Default is TRUE. Will fall back to R implementation if C++ is not available.

#' Internal function to load sample sheet
#'
#' @param sample_sheet Path to sample sheet file or data.frame/data.table
#' @return A data.table containing the sample information
#' @keywords internal
#' @noRd
load_sample_sheet <- function(sample_sheet) {
  if (is.character(sample_sheet) && length(sample_sheet) == 1) {
    message("Loading sample sheet from: ", sample_sheet)
    sample_data <- data.table::fread(sample_sheet)
  } else if (is.data.frame(sample_sheet)) {
    # Convert to data.table if it's a data.frame
    if (!data.table::is.data.table(sample_sheet)) {
      sample_data <- data.table::as.data.table(sample_sheet)
    } else {
      sample_data <- sample_sheet
    }
  } else {
    stop("sample_sheet must be either a file path or a data.frame/data.table")
  }

  # Validate sample sheet format
  required_cols <- c("group_id", "replicate", "sample_id", "file_name")
  missing_cols <- setdiff(required_cols, colnames(sample_data))
  if (length(missing_cols) > 0) {
    stop("Sample sheet is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  return(sample_data)
}


#' Internal function to filter sample data by groups
#'
#' @param sample_data A data.table containing sample information
#' @param groups Vector of group IDs to include
#' @return Filtered sample_data
#' @keywords internal
#' @noRd
filter_by_groups <- function(sample_data, groups) {
  if (!is.null(groups)) {
    # Check if requested groups exist in the sample sheet
    missing_groups <- setdiff(groups, unique(sample_data$group_id))
    if (length(missing_groups) > 0) {
      warning("The following requested groups are not in the sample sheet: ",
              paste(missing_groups, collapse = ", "))
    }

    # Filter the sample data
    original_samples <- nrow(sample_data)
    sample_data <- sample_data[group_id %in% groups]

    if (nrow(sample_data) == 0) {
      stop("No samples found for the requested groups: ", paste(groups, collapse = ", "))
    }

    message(sprintf("Filtered to %d samples from %d groups",
                    nrow(sample_data), length(intersect(groups, unique(sample_data$group_id)))))
  }

  return(sample_data)
}


#' Internal function to filter sample data by type
#'
#' @param sample_data A data.table containing sample information
#' @param type character string of methylation type, e.g. GCH or HCG
#' @return Filtered sample_data
#' @keywords internal
#' @noRd
filter_by_type <- function(sample_data, type) {
  if (!is.null(type)) {
    original_samples <- nrow(sample_data)
    sample_data <- sample_data[grep(type, file_name)]
    if (nrow(sample_data) == 0) {
      stop(sprintf("No %s files found in sample sheet.", type))
    }
    message(sprintf("Filtered sample sheet to %d %s files", nrow(sample_data), type))
  }

  return(sample_data)
}

#' Internal function to validate input files
#'
#' @param files.list list or character vector of file names
#' @return list or character vector of file names
#' @keywords internal
#' @noRd
check_files_exist <- function(files.list) {
  existing_files <- file.exists(files.list)
  if (!all(existing_files)) {
    warning(sprintf("%d files from sample sheet not found in directory",
                    sum(!existing_files)))
    message("Missing files: ", paste(files.list[!existing_files], collapse=", "))
    files.list <- files.list[existing_files]
  }

  # Validation step
  if (length(files.list) == 0) {
    stop("No files found in the directory matching the sample sheet.")
  }
  if (length(files.list) < 2) {
    stop("At least two files are required for analysis.")
  }

  message(sprintf("Loading %d files...", length(files.list)))

  return(files.list)
}

#' Internal function to load a single file using R
#' @return data.table of methylation data
#' @keywords internal
#' @noRd
#'
load_single_file_r <- function(file) {
  tryCatch({
    # Use base R to read the gzipped file
    con <- gzfile(file, "r")
    data_lines <- readLines(con)
    close(con)

    # Split the tab-delimited lines
    split_lines <- strsplit(data_lines, "\t")

    # Extract the columns we need
    n_rows <- length(split_lines)

    # Preallocate vectors for efficiency
    chr <- character(n_rows)
    pos <- integer(n_rows)
    strand <- character(n_rows)
    site <- character(n_rows)
    mc <- integer(n_rows)
    cov <- integer(n_rows)

    # Extract data column by column
    for (i in 1:n_rows) {
      line <- split_lines[[i]]
      chr[i] <- line[1]
      pos[i] <- as.integer(line[2])
      strand[i] <- line[3]
      site[i] <- line[4]
      mc[i] <- as.integer(line[5])
      cov[i] <- as.integer(line[6])
    }

    # Create the data.table
    df <- data.table::data.table(
      chr = chr,
      pos = pos,
      strand = strand,
      site = site,
      mc = mc,
      cov = cov
    )

    # Safely calculate rate
    df[, rate := ifelse(cov > 0, mc/cov, 0)]
    df[, uniqueID := paste(chr, pos, site, sep="_")]

    return(df)
  }, error = function(e) {
    message("Error processing file: ", file)
    message("Error message: ", e$message)
    return(NULL)
  })
}

#' Internal function to sequentially load files using R implementation
#' @return data.table of methylation data
#' @keywords internal
#' @noRd
#'
load_data_sequential <- function(files.list) {
  # Create progress bar
  pb <- progress::progress_bar$new(
    format = "  Loading [:bar] :percent in :elapsed",
    total = length(files.list),
    clear = FALSE
  )

  all_samples <- lapply(files.list, function(file) {
    pb$tick()
    load_single_file_r(file)
  })

  names(all_samples) <- names(files.list)
  return(all_samples)
}

#' Internal function to sequentially load files in parallel using R implementation
#' @return data.table of methylation data
#' @keywords internal
#' @noRd
#'
load_data_parallel <- function(files.list, cores) {
  all_samples <- NULL

  if (cores > 1) {
    # Parallel processing
    tryCatch({
      message("Cluster created at: ", Sys.time())
      cl <- parallel::makeCluster(min(cores, parallel::detectCores() - 1))
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Export the loading function
      parallel::clusterExport(cl = cl, varlist = c("load_file"), envir = environment())

      parallel::clusterEvalQ(cl, {
        library(data.table)
        NULL
      })

      all_samples <- parallel::parLapply(cl, files.list, load_single_file_r)
      names(all_samples) <- names(files.list)

    }, error = function(e) {
      message("Error in parallel processing: ", e$message)
      message("Falling back to sequential processing...")
      return(NULL)
    })
  }

  # If parallel processing failed, fall back to sequential
  if (is.null(all_samples)) {
    all_samples <- load_data_sequential(files.list)
  }

  return(all_samples)
}

#' Internal function post-process loaded samples
#' @return data.table, or list of data.tables, of methylation data
#' @keywords internal
#' @noRd
#'
post_process_samples <- function(all_samples, sample_data, files.list) {
  # Common post-processing for both implementations
  valid_samples <- !sapply(all_samples, is.null)
  all_samples <- all_samples[valid_samples]

  if (length(all_samples) == 0) {
    warning("No files were successfully loaded!")
    attr(all_samples, "sample_metadata") <- sample_data[0,]
    return(all_samples)
  }

  valid_sample_ids <- names(files.list)[valid_samples]
  filtered_sample_data <- sample_data[match(valid_sample_ids, sample_data$sample_id)]
  attr(all_samples, "sample_metadata") <- filtered_sample_data

  return(all_samples)
}
