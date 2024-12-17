#' Normalize Methylation Rates
#'
#' Normalizes methylation rates using quantile-based normalization to account for
#' non-linear enzyme efficiency differences.
#'
#' @param data_list List of normalized coverage data (output from normalize_coverage_data)
#' @param group_names Character vector of group names
#' @param n_quantiles Number of quantiles to use (default = 10 for deciles)
#' @param between_groups Logical indicating whether to normalize between groups
#' @param method Character indicating normalization method ("quantile" or "loess")
#'
#' @return List of lists containing normalized methylation data
#'
#' @import data.table
#' @import purrr
#'
#' @export
#'
normalize_methylation_rates <- function(data_list,
                                        group_names,
                                        n_quantiles = 10,
                                        between_groups = FALSE,
                                        method = "quantile") {

  if(between_groups) {
    # Flatten list for between-group normalization
    all_samples <- unlist(data_list, recursive = FALSE)
    normalized_data <- normalize_methylation_within_set(all_samples,
                                                        n_quantiles = n_quantiles,
                                                        method = method)

    # Reconstruct the group structure
    result <- list()
    current_idx <- 1
    for(group in names(data_list)) {
      n_samples <- length(data_list[[group]])
      result[[group]] <- normalized_data[current_idx:(current_idx + n_samples - 1)]
      current_idx <- current_idx + n_samples
    }

  } else {
    # Normalize within each group separately
    result <- lapply(data_list, function(group) {
      normalize_methylation_within_set(group,
                                       n_quantiles = n_quantiles,
                                       method = method)
    })
  }

  return(result)
}

#' Internal function to normalize methylation rates within a set of samples
#'
#' @keywords internal
normalize_methylation_within_set <- function(sample_list, n_quantiles, method) {
  # Calculate average methylation rates
  avg_rates <- rowMeans(do.call(cbind, lapply(sample_list, function(df) df$rate)))

  if(method == "quantile") {
    # Create quantiles
    probs <- seq(0, 1, length.out = n_quantiles + 1)
    rate_quantiles <- cut(avg_rates,
                          quantile(avg_rates, prob = probs, names = FALSE),
                          include.lowest = TRUE)

    # Initialize result list
    result <- sample_list

    # Normalize each quantile
    for(quant in levels(rate_quantiles)) {
      sites_in_quantile <- which(rate_quantiles == quant)

      # Sum methylation counts for these sites
      quant_sums <- sapply(result, function(df) sum(df$mc[sites_in_quantile]))
      avg_quant_sum <- mean(quant_sums)

      # Calculate and apply scaling factors
      sf_quant <- quant_sums / avg_quant_sum

      result <- map2(result, sf_quant, function(df, sf) {
        dt <- data.table::copy(df)
        dt[sites_in_quantile, mc := mc / sf]
        dt[sites_in_quantile, rate := mc / cov]
        dt
      })
    }

  } else if(method == "loess") {
    # Alternative: LOESS normalization
    # This would implement a smoothed normalization; need to talk to Mike and Rhonda about this
    stop("LOESS normalization not yet implemented")
  }

  return(result)
}
