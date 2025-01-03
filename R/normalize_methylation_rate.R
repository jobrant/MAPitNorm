#' Normalize Methylation Rates
#'
#' Normalizes methylation rates using quantile-based normalization to account for
#' non-linear enzyme efficiency differences.
#'
#' @param data_list List of normalized coverage data
#' @param group_names Character vector of group names
#' @param within_groups Logical indicating whether to normalize within groups (default = TRUE)
#' @param between_groups Logical indicating whether to normalize between groups (default = TRUE)
#' @param within_alpha Numeric between 0 and 1 controlling structure preservation within groups (default = 0.3)
#' @param between_alpha Numeric between 0 and 1 controlling structure preservation between groups (default = 0.5)
#' @param sites_per_quantile Target number of sites per quantile (default = 1000).
#'        For large datasets (>100k sites), this provides a good balance between
#'        normalization granularity and statistical power.
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param diagnostics Logical indicating whether to print diagnostic information
#'
#' @details This function performs a two-step normalization:
#'          1. Within-group normalization to handle technical variation between replicates
#'          2. Between-group normalization to adjust for global efficiency differences
#'          while preserving biological differences. The higher between_alpha preserves
#'          more of the biological differences between groups.
#'
#' @return List of lists containing normalized methylation data
#'
#' @importFrom data.table setDT copy :=
#' @importFrom stats quantile
#' @importFrom purrr map2
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
normalize_methylation_rates <- function(data_list,
                                        group_names,
                                        within_groups = TRUE,
                                        between_groups = TRUE,
                                        within_alpha = 0.3,
                                        between_alpha = 0.5,
                                        sites_per_quantile = 1000,
                                        max_quantiles = 50,
                                        diagnostics = TRUE) {

  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (!all(sapply(data_list, is.list))) stop("data_list must be a list of lists")

  result <- data_list

  # Step 1: Within-group normalization
  if(within_groups) {
    if(diagnostics) cat("\n=== Performing within-group normalization ===\n")
    result <- lapply(names(result), function(group) {
      if(diagnostics) cat(sprintf("\nProcessing group: %s\n", group))
      normalize_methylation_within_set(
        result[[group]],
        alpha = within_alpha,
        sites_per_quantile = sites_per_quantile,
        max_quantiles = max_quantiles,
        diagnostics = diagnostics
      )
    })
    names(result) <- names(data_list)
  }

  # Step 2: Between-group normalization
  if(between_groups) {
    if(diagnostics) cat("\n=== Performing between-group normalization ===\n")
    # Flatten data for between-group normalization
    all_samples <- unlist(result, recursive = FALSE)

    # Perform gentler between-group normalization
    normalized_data <- normalize_methylation_within_set(
      all_samples,
      alpha = between_alpha,  # Higher alpha to preserve more group differences
      sites_per_quantile = sites_per_quantile * 2,  # Fewer quantiles for between-group
      max_quantiles = max_quantiles,
      diagnostics = diagnostics
    )

    # Reconstruct the group structure
    result <- list()
    current_idx <- 1
    for(group in names(data_list)) {
      n_samples <- length(data_list[[group]])
      result[[group]] <- normalized_data[current_idx:(current_idx + n_samples - 1)]
      current_idx <- current_idx + n_samples
    }
  }

  names(result) <- names(data_list)
  return(result)
}


#' Internal function to normalize methylation rates within a set of samples
#' @param sample_list List of replicate data frames
#' @param alpha Numeric between 0 and 1 controlling structure preservation (default = 0.3)
#' @param sites_per_quantile Target number of sites per quantile (default = 1000)
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param diagnostics Logical indicating whether to print diagnostic information
#'
#' @keywords internal

normalize_methylation_within_set <- function(sample_list,
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

    if(length(sites_in_quantile) > 0) {
      # Calculate average rate for this quantile
      quant_rates <- sapply(result, function(df) mean(df$rate[sites_in_quantile]))
      avg_quant_rate <- mean(quant_rates)

      # Calculate and apply scaling factors with structure preservation
      result <- map2(result, quant_rates, function(df, sample_rate) {
        dt <- data.table::copy(df)
        data.table::setDT(dt)

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
    }

    if(diagnostics) setTxtProgressBar(pb, i)
  }

  if(diagnostics) {
    close(pb)
    cat("\nNormalization complete!\n")
  }

  return(result)
}
