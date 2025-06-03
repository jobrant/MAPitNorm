#' Load and Preprocess Data
#'
#' Load methylation data from files based on a sample sheet.
#'
#' @param dir_path Path to the directory containing the data files.
#' @param sample_sheet Path to the sample sheet file or a data.frame/data.table
#' containing sample information with columns: group_id, replicate, sample_id, file_name.
#' @param type Character string specifying which type of files to load. If NULL (default),
#' it will use file names from the sample sheet as is.
#' @param groups Character vector of group_ids to include. If NULL (default), all
#' groups will be loaded.
#' @param cores Number of cores to use for parallel processing. Default is 1
#'   (no parallel processing).
#'
#' @return A list of data.tables containing processed methylation data.
#'
#' @importFrom stringr str_sub
#' @importFrom progress progress_bar
#' @importFrom data.table fread setnames := is.data.table as.data.table
#' @importFrom parallel makeCluster stopCluster parLapply detectCores clusterExport
#' @importFrom utils read.delim
#'
#' @examples
#' \dontrun{
#' # Load all files based on sample sheet
#' all_samples <- load_data("path/to/files", "path/to/sample_sheet.tsv")
#'
#' # Load only samples from specific groups
#' group1_samples <- load_data("path/to/files", "path/to/sample_sheet.tsv",
#'                             groups = c("M1", "M2"))
#'
#' # Combine group filtering with type filtering
#' group1_hcg <- load_data("path/to/files", "path/to/sample_sheet.tsv",
#'                         type = "HCG", groups = "M1")
#' }
#'
#' @export
#'

load_data <- function(dir_path, sample_sheet, type = NULL, groups = NULL, cores = 1) {
  # Load the sample sheet if it's a file path
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

  # Filter by groups if specified
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

  # Filter for specific type if provided
  if (!is.null(type)) {
    original_samples <- nrow(sample_data)
    sample_data <- sample_data[grep(type, file_name)]
    if (nrow(sample_data) == 0) {
      stop(sprintf("No %s files found in sample sheet.", type))
    }
    message(sprintf("Filtered sample sheet to %d %s files", nrow(sample_data), type))
  }

  # Construct file paths
  files.list <- file.path(dir_path, sample_data$file_name)
  names(files.list) <- sample_data$sample_id

  # Check if files exist
  existing_files <- file.exists(files.list)
  if (!all(existing_files)) {
    warning(sprintf("%d files from sample sheet not found in directory",
                    sum(!existing_files)))
    message("Missing files: ", paste(files.list[!existing_files], collapse=", "))
    sample_data <- sample_data[existing_files]
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

  # Define columns we need - note we know there's no header from inspection
  cols_needed <- c("chr", "pos", "strand", "site", "mc", "cov", "context")

  message(sprintf("Starting processing with %d cores", cores))
  start_time <- Sys.time()

  # Define an optimized loading function based on what we now know about the files
  load_file <- function(file) {
    tryCatch({
      # Use base R to read the gzipped file since fread has issues
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
      # This avoids creating a large intermediate matrix
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

  all_samples <- list()

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

      all_samples <- parallel::parLapply(cl, files.list, load_file)

    }, error = function(e) {
      message("Error in parallel processing: ", e$message)
      message("Falling back to sequential processing...")
      return(NULL)
    })

    # If parallel processing failed, fall back to sequential
    if (is.null(all_samples)) {
      cores <- 1
    }
  }

  if (cores == 1) {
    pb <- progress::progress_bar$new(
      format = "  Loading [:bar] :percent in :elapsed",
      total = length(files.list),
      clear = FALSE
    )

    all_samples <- lapply(files.list, function(file) {
      pb$tick()
      load_file(file)
    })
  }

  # Remove any NULL entries from failed processing
  valid_samples <- !sapply(all_samples, is.null)
  all_samples <- all_samples[valid_samples]

  if (length(all_samples) == 0) {
    warning("No files were successfully loaded!")
    # Still attach an empty metadata
    attr(all_samples, "sample_metadata") <- sample_data[0,]
    return(all_samples)
  }

  # Get the sample IDs that were successfully loaded
  valid_sample_ids <- names(files.list)[valid_samples]

  # Filter the sample data to match successfully loaded samples
  filtered_sample_data <- sample_data[match(valid_sample_ids, sample_data$sample_id)]

  # Attach sample metadata to the returned list as an attribute
  attr(all_samples, "sample_metadata") <- filtered_sample_data

  # Summarize loaded files by group
  group_counts <- table(filtered_sample_data$group_id)
  message("Files loaded by group:")
  for (grp in names(group_counts)) {
    message("  Group ", grp, ": ", group_counts[grp], " files")
  }

  message("Successfully loaded ", length(all_samples), " out of ", length(files.list), " files")

  # Calculate total elapsed time
  end_time <- Sys.time()
  message("Total time: ", format(end_time - start_time))

  gc(verbose = FALSE)
  return(all_samples)
}
