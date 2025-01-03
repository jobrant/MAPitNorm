#' Normalizes methylation coverage data within groups (replicates) and optionally between groups.
#'
#' @param data_list A list of lists containing methylation data frames. Each inner list
#'        represents a group and contains one or more replicate data frames. Each data frame
#'        must contain columns 'cov' and 'mc'.
#' @param group_names A character vector of group names corresponding to the names in data_list.
#' @param between_groups Logical indicating whether to normalize between groups.
#'        Default is FALSE to preserve biological differences between groups.
#' @param diagnostics Logical indicating whether to print diagnostic information.
#'        Default is FALSE.
#'
#' @details When between_groups = FALSE, each group is normalized independently,
#'          preserving potential biological differences in total accessibility between groups.
#'          When between_groups = TRUE, all samples are normalized to the same mean coverage,
#'          which might mask global accessibility changes between conditions.
#'
#' @return A list of lists with the same structure as the input, containing normalized
#'         coverage data. Each data frame in the output will contain normalized versions
#'         of 'cov' and 'mc' columns, plus a calculated 'rate' column (mc/cov).
#'
#' @examples
#' \dontrun{
#' split_data <- split_by_groups(gch_samples, group_names)
#' normalized_coverage <- normalize_coverage(split_data, group_names)
#' }
#' @importFrom data.table setDT copy := setDTthreads
#' @importFrom purrr map2 transpose
#' @importFrom parallel detectCores
#'
#' @export

normalize_coverage <- function(data_list, group_names, between_groups = FALSE,  diagnostics = FALSE) {
  data.table::setDTthreads(threads = parallel::detectCores())  # Optional: for better performance

  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (!all(sapply(data_list, is.list))) stop("data_list must be a list of lists")
  if (!all(sapply(unlist(data_list, recursive = FALSE),
                  function(df) all(c("cov", "mc") %in% names(df))))) {
    stop("All data frames must contain 'cov' and 'mc' columns")
  }

  # Create a named vector to store scaling factors
  sf_total <- numeric()

  if(between_groups) {
    # Normalize across all samples
    for(group in names(data_list)) {
      for(sample_name in names(data_list[[group]])) {
        sf_total[sample_name] <- sum(data_list[[group]][[sample_name]]$cov)
      }
    }
    avg_total <- mean(sf_total)
    sf_total <- sf_total / avg_total

    print("Performing between-group normalization")
  } else {
    # Normalize within each group separately
    for(group in names(data_list)) {
      group_totals <- sapply(data_list[[group]], function(df) sum(df$cov))
      group_mean <- mean(group_totals)
      group_sf <- group_totals / group_mean
      sf_total[names(data_list[[group]])] <- group_sf
    }

    print("Performing within-group normalization only")
  }

  # Print diagnostics
  print("Scaling factors:")
  print(sf_total)

  # Normalize each group
  normalized_groups <- lapply(names(data_list), function(group_name) {
    group_data <- data_list[[group_name]]

    # Get scaling factors for this group using names
    group_sf <- sf_total[names(group_data)]

    print(paste("Processing group:", group_name))
    print("Group scaling factors:")
    print(group_sf)

    .normalize_group(group_data, group_sf)
  })

  # Add diagnostic check
  if(diagnostics) {
    # Sample a few values before normalization
    sample_before <- lapply(data_list, function(group) {
      lapply(group, function(df) {
        head(df$cov, 5)  # First 5 coverage values
      })
    })

    # After normalization
    sample_after <- lapply(normalized_groups, function(group) {
      lapply(group, function(df) {
        head(df$cov, 5)  # First 5 coverage values
      })
    })

    cat("\nSample values before normalization:\n")
    print(sample_before)
    cat("\nSample values after normalization:\n")
    print(sample_after)
  }

  names(normalized_groups) <- names(data_list)
  return(normalized_groups)
}

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
