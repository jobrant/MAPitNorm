#' @title Internal Utility Functions
#' @description
#' This file contains internal utility functions that are not exported

# This is a special note about parameters that only exist in the main function
# and not in these internal helpers, but are documented for clarity:
# - use_cpp is a parameter in load_data() that controls whether to use C++ implementation


# Data Loading Utilities --------------------------------------------------

#' Internal function to load sample sheet
#' @param sample_sheet Path to sample sheet file or data.frame/data.table
#' @return A data.table containing the sample information
#' @keywords internal
#' @noRd
.load_sample_sheet <- function(sample_sheet) {
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
#' @param sample_data A data.table containing sample information
#' @param groups Vector of group IDs to include
#' @return Filtered sample_data
#' @keywords internal
#' @noRd
.filter_by_groups <- function(sample_data, groups) {
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
#' @param sample_data A data.table containing sample information
#' @param type character string of methylation type, e.g. GCH or HCG
#' @return Filtered sample_data
#' @keywords internal
#' @noRd
.filter_by_type <- function(sample_data, type) {
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
#' @param files.list list or character vector of file names
#' @return list or character vector of file names
#' @keywords internal
#' @noRd
.check_files_exist <- function(files.list) {
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

# Internal function to load using C++ implementation
.load_data_cpp <- function(files.list, cores) {
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
.load_data_r <- function(files.list, cores) {
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

#' Internal function to load a single file using R
#' @return data.table of methylation data
#' @keywords internal
#' @noRd
.load_single_file_r <- function(file) {
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
.load_data_sequential <- function(files.list) {
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
.load_data_parallel <- function(files.list, cores) {
  all_samples <- NULL

  if (cores > 1) {
    # Parallel processing
    tryCatch({
      message("Cluster created at: ", Sys.time())
      cl <- parallel::makeCluster(min(cores, parallel::detectCores() - 1))
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Export the loading function
      parallel::clusterExport(cl = cl, varlist = c("load_single_file_r"), envir = environment())

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
.post_process_samples <- function(all_samples, sample_data, files.list) {
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


# Coverage and Methylation Normalization Utilities ------------------------

#' Internal function to normalize a group of samples
#' @param replicates List of replicate data frames
#' @param sf Vector of scaling factors
#' @return List of normalized data frames
#' @keywords internal
.normalize_group <- function(replicates, sf) {
  purrr::map2(replicates, sf, function(df, s) {
    if(is.na(s)) stop("Scaling factor is NA")
    dt <- data.table::copy(df)
    data.table::setDT(dt)
    data.table::setDT(dt)[, c("cov", "mc", "rate") := list(cov/s,
                                                           mc/s,
                                                           mc/cov
    )]
    return(dt)
  })
}

#' Internal function to normalize methylation rates within a set of samples
#' @param sample_list List of replicate data frames
#' @param alpha Numeric between 0 and 1 controlling structure preservation (default = 0.3)
#' @param sites_per_quantile Target number of sites per quantile (default = 1000)
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param diagnostics Logical indicating whether to print diagnostic information
#' @keywords internal
.normalize_methylation_within_set <- function(sample_list,
                                             alpha = 0.3,
                                             sites_per_quantile = 1000,
                                             max_quantiles = 50,
                                             diagnostics = TRUE) {

  if(diagnostics) cat("Calculating average rates...\n")

  # Calculate average methylation rates
  avg_rates <- rowMeans(do.call(cbind, lapply(sample_list, function(df) df$rate)))

  # Dynamic quantile calculation
  n_sites <- length(avg_rates)
  n_quantiles <- as.integer(max(5, min(n_sites/sites_per_quantile, max_quantiles)))
  if(diagnostics) cat(sprintf("Using %d quantiles for %d sites\n", n_quantiles, n_sites))

  # Create quantiles
  probs <- seq(0, 1, length.out = n_quantiles + 1)
  breaks <- unique(quantile(avg_rates, prob = probs, names = FALSE))

  if(length(breaks) < 3) {
    warning("Too few unique values for quantile binning. Reducing number of quantiles.")
    n_quantiles <- length(unique(avg_rates)) - 1
    if(n_quantiles < 2) stop("Not enough unique values for quantile normalization")
    probs <- seq(0, 1, length.out = n_quantiles + 1)
    breaks <- unique(quantile(avg_rates, prob = probs, names = FALSE))
  }

  rate_quantiles <- cut(avg_rates, breaks = breaks, include.lowest = TRUE)

  # Check that quantiles are valid and match data size
  if (length(rate_quantiles) != n_sites) {
    stop(sprintf("Quantile vector length (%d) doesn't match site count (%d)",
                 length(rate_quantiles), n_sites))
  }

  # Initialize result list and progress bar
  result <- sample_list
  if(diagnostics) {
    cat("Normalizing by quantile...\n")
    pb <- txtProgressBar(min = 0, max = length(levels(rate_quantiles)), style = 3)
  }

  # Normalize each quantile
  for(i in seq_along(levels(rate_quantiles))) {
    quant <- levels(rate_quantiles)[i]
    sites_in_quantile <- which(rate_quantiles == quant)

    # Check for empty quantiles
    if(length(sites_in_quantile) == 0) {
      if(diagnostics) cat(sprintf("Warning: Quantile %s is empty, skipping\n", quant))
      if(diagnostics) setTxtProgressBar(pb, i)
      next
    }

    # Validate site indices
    max_sites <- min(sapply(result, nrow))
    if(any(sites_in_quantile > max_sites)) {
      problem_sites <- sites_in_quantile[sites_in_quantile > max_sites]
      if(diagnostics) cat(sprintf("Warning: Removing %d out-of-bounds sites\n",
                                  length(problem_sites)))
      sites_in_quantile <- sites_in_quantile[sites_in_quantile <= max_sites]

      # Skip if no valid sites remain
      if(length(sites_in_quantile) == 0) {
        if(diagnostics) setTxtProgressBar(pb, i)
        next
      }
    }

    # Calculate average rate for this quantile
    quant_rates <- sapply(result, function(df) mean(df$rate[sites_in_quantile]))
    avg_quant_rate <- mean(quant_rates)

    # Calculate and apply scaling factors with structure preservation
    result <- map2(result, quant_rates, function(df, sample_rate) {
      dt <- data.table::copy(df)
      data.table::setDT(dt)

      # Safety check for this specific data table
      if(max(sites_in_quantile) > nrow(dt)) {
        valid_sites <- sites_in_quantile[sites_in_quantile <= nrow(dt)]
        if(length(valid_sites) == 0) {
          return(dt)  # Return unchanged if no valid sites
        }
        sites_in_quantile <- valid_sites
      }

      # Get current rates for these sites
      current_rates <- dt$rate[sites_in_quantile]

      # Preserve structure while normalizing to average
      centered <- current_rates - mean(current_rates)
      new_rates <- avg_quant_rate + alpha * centered

      # Update methylation counts and rates
      dt[sites_in_quantile, `:=`(
        mc = cov * new_rates,
        rate = new_rates
      )]

      return(dt)
    })

    if(diagnostics) setTxtProgressBar(pb, i)
  }

  if(diagnostics) {
    close(pb)
    cat("\nNormalization complete!\n")
  }

  return(result)
}

