#' Normalize Methylation Rates
#'
#' Normalizes methylation rates using quantile-based normalization to account for
#' non-linear enzyme efficiency differences.
#'
#' @param data_list List of normalized coverage data
#' @param group_names Character vector of group names
#' @param alpha Numeric between 0 and 1 controlling structure preservation (default = 0.3)
#' @param use_smoothing Logical indicating whether to apply smoothing (default = TRUE)
#' @param min_quantiles Minimum number of quantiles to use (default = 5)
#' @param sites_per_quantile Target number of sites per quantile (default = 20)
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param between_groups Logical indicating whether to normalize between groups
#' @param diagnostics Logical indicating whether to print diagnostic information
#'
#' @details This function normalizes methylation rates within quantile bins to account
#' for non-linear enzyme efficiency differences. When between_groups = FALSE, normalization
#' is performed within each group separately to preserve biological differences between groups.
#'
#' @return List of lists containing normalized methylation data
#'
#' @importFrom data.table setDT copy :=
#' @importFrom purrr map2
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
normalize_methylation_rates <- function(data_list,
                                        group_names,
                                        alpha = 0.3,
                                        use_smoothing = TRUE,
                                        min_quantiles = 5,
                                        sites_per_quantile = 20,
                                        max_quantiles = 50,
                                        between_groups = FALSE,
                                        diagnostics = TRUE) {

  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (!all(sapply(data_list, is.list))) stop("data_list must be a list of lists")

  if(between_groups) {
    if(diagnostics) cat("Performing between-group normalization\n")
    all_samples <- unlist(data_list, recursive = FALSE)
    normalized_data <- normalize_methylation_within_set(
      all_samples,
      alpha = alpha,
      use_smoothing = use_smoothing,
      min_quantiles = min_quantiles,
      sites_per_quantile = sites_per_quantile,
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
  } else {
    if(diagnostics) cat("Performing within-group normalization\n")
    result <- lapply(data_list, function(group) {
      normalize_methylation_within_set(
        group,
        alpha = alpha,
        use_smoothing = use_smoothing,
        min_quantiles = min_quantiles,
        sites_per_quantile = sites_per_quantile,
        max_quantiles = max_quantiles,
        diagnostics = diagnostics
      )
    })
  }

  names(result) <- names(data_list)
  return(result)
}


#' Internal function to normalize methylation rates within a set of samples
#' @param sample_list List of replicate data frames
#' @param alpha Numeric between 0 and 1 controlling structure preservation (default = 0.3)
#' @param use_smoothing Logical indicating whether to apply smoothing (default = TRUE)
#' @param min_quantiles Minimum number of quantiles to use (default = 5)
#' @param sites_per_quantile Target number of sites per quantile (default = 20)
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param diagnostics Logical indicating whether to print diagnostic information
#' @importFrom stats ksmooth quantile
#' @keywords internal

normalize_methylation_within_set <- function(sample_list,
                                             alpha = 0.3,
                                             use_smoothing = TRUE,
                                             min_quantiles = 5,
                                             sites_per_quantile = 20,
                                             max_quantiles = 50,
                                             diagnostics = TRUE) {

  if(diagnostics) cat("Calculating average rates...\n")
  # Calculate average methylation rates
  avg_rates <- rowMeans(do.call(cbind, lapply(sample_list, function(df) df$rate)))

  # Dynamic quantile calculation
  n_sites <- length(avg_rates)
  n_quantiles <- as.integer(max(min_quantiles, min(n_sites/sites_per_quantile, max_quantiles)))
  if(diagnostics) cat(sprintf("Using %d quantiles for %d sites\n", n_quantiles, n_sites))

  # Create quantiles
  if(diagnostics) cat("Creating quantile bins...\n")
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

      # Apply hybrid normalization to each sample
      result <- map2(result, quant_rates, function(df, sample_rate) {
        dt <- data.table::copy(df)
        data.table::setDT(dt)

        # Get current rates for these sites
        current_rates <- dt$rate[sites_in_quantile]

        if(use_smoothing && length(sites_in_quantile) > 3) {  # Only smooth if enough points
          # Apply smoothing with explicit output points
          bw <- diff(range(sites_in_quantile))/n_quantiles
          smoothed <- ksmooth(sites_in_quantile,
                              current_rates,
                              bandwidth = bw,
                              x.points = sites_in_quantile)$y

          # Combine smoothing with structure preservation
          centered <- current_rates - smoothed
          new_rates <- smoothed + alpha * centered
        } else {
          # Just preserve structure
          centered <- current_rates - mean(current_rates)
          new_rates <- avg_quant_rate + alpha * centered
        }

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
