#' Normalize Methylation Rates with Method Comparison
#'
#' @param data_list List of normalized coverage data
#' @param group_names Character vector of group names
#' @param n_quantiles Number of quantiles for quantile normalization
#' @param between_groups Logical for between-group normalization
#' @param compare_methods Logical to generate comparison plots
#' @param loess_span Numeric value for LOESS smoothing parameter (0.2-0.8)
#' @param sample_size Number of sites to sample for method comparison
#' @param testing Logical indicating whether to use sampling for actual normalization
#'
#' @import data.table
#' @import ggplot2
#' @import purrr
#'
#' @export
#'
normalize_methylation_rates_compare <- function(data_list,
                                                group_names,
                                                n_quantiles = 10,
                                                compare_methods = TRUE,
                                                loess_span = 0.2,
                                                sample_size = 100000,
                                                testing = FALSE) {

  # Function to sample data
  sample_data <- function(data, size) {
    # Get consistent sample indices
    set.seed(42)  # for reproducibility
    sample_idx <- sample(nrow(data[[1]][[1]]), min(size, nrow(data[[1]][[1]])))

    # Apply sampling to all groups and replicates
    lapply(data, function(group) {
      lapply(group, function(df) df[sample_idx,])
    })
  }

  # Sample data if testing or comparing methods
  if(testing || compare_methods) {
    cat(sprintf("Sampling %d sites for analysis...\n", sample_size))
    working_data <- sample_data(data_list, sample_size)
  } else {
    working_data <- data_list
  }

  if(compare_methods) {
    cat("Comparing normalization methods...\n")

    # Store results for both methods
    results_quantile <- normalize_methylation_within_set(working_data[[1]],
                                                         n_quantiles = n_quantiles,
                                                         method = "quantile",
                                                         loess_span = loess_span)

    results_loess <- normalize_methylation_within_set(working_data[[1]],
                                                      method = "loess",
                                                      loess_span = loess_span)

    # Generate comparison plots
    compare_normalization_methods(working_data[[1]][[1]],
                                  results_quantile[[1]],
                                  results_loess[[1]])
  }

  # Return list with both normalizations if comparing, otherwise just quantile
  if(compare_methods) {
    return(list(
      quantile = results_quantile,
      loess = results_loess
    ))
  } else {
    return(normalize_methylation_within_set(working_data[[1]],
                                            n_quantiles = n_quantiles,
                                            method = "quantile"))
  }
}

#' Internal function for normalization methods
#' @keywords internal
normalize_methylation_within_set <- function(sample_list,
                                             n_quantiles = 10,
                                             method = "quantile",
                                             loess_span = 0.2) {

  # Calculate average rates
  cat("Calculating average rates...\n")
  avg_rates <- rowMeans(do.call(cbind, lapply(sample_list, function(df) df$rate)))

  if(method == "quantile") {
    cat("Performing quantile normalization...\n")
    probs <- seq(0, 1, length.out = n_quantiles + 1)
    rate_quantiles <- cut(avg_rates,
                          quantile(avg_rates, prob = probs, names = FALSE),
                          include.lowest = TRUE)

    result <- sample_list

    # Add progress bar for quantiles
    pb <- txtProgressBar(min = 0, max = length(levels(rate_quantiles)), style = 3)

    for(i in seq_along(levels(rate_quantiles))) {
      quant <- levels(rate_quantiles)[i]
      sites_in_quantile <- which(rate_quantiles == quant)
      quant_sums <- sapply(result, function(df) sum(df$mc[sites_in_quantile]))
      avg_quant_sum <- mean(quant_sums)
      sf_quant <- quant_sums / avg_quant_sum

      result <- map2(result, sf_quant, function(df, sf) {
        dt <- data.table::copy(df)
        dt[sites_in_quantile, mc := mc / sf]
        dt[sites_in_quantile, rate := mc / cov]
        dt
      })

      setTxtProgressBar(pb, i)
    }

    close(pb)
    cat("\nQuantile normalization complete!\n")

  } else if(method == "loess") {
    cat("Performing LOESS normalization...\n")
    pb <- txtProgressBar(min = 0, max = length(sample_list), style = 3)

    result <- lapply(seq_along(sample_list), function(i) {
      dt <- data.table::copy(sample_list[[i]])

      loess_fit <- loess(dt$rate ~ avg_rates, span = loess_span)
      predicted_rates <- predict(loess_fit)

      dt[, c("mc", "rate") := list(cov * predicted_rates, predicted_rates)]

      setTxtProgressBar(pb, i)
      return(dt)
    })

    close(pb)
    cat("\nLOESS normalization complete!\n")
  }

  return(result)
}

#' Compare Normalization Methods
#'
#' @param original_data Original data.table
#' @param quantile_norm Quantile normalized data
#' @param loess_norm LOESS normalized data
#'
#' @import ggplot2
#' @export
compare_normalization_methods <- function(original_data, quantile_norm, loess_norm) {
  # Prepare data for plotting
  plot_data <- data.frame(
    original = original_data$rate,
    quantile = quantile_norm$rate,
    loess = loess_norm$rate
  )

  # Scatter plot
  p1 <- ggplot(plot_data) +
    geom_point(aes(x = original, y = quantile, color = "Quantile"), alpha = 0.5) +
    geom_point(aes(x = original, y = loess, color = "LOESS"), alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(x = "Original Methylation Rate",
         y = "Normalized Rate",
         title = "Comparison of Normalization Methods",
         color = "Method") +
    theme_minimal()

  # Density plot
  p2 <- ggplot(plot_data) +
    geom_density(aes(x = original, color = "Original"), alpha = 0.5) +
    geom_density(aes(x = quantile, color = "Quantile"), alpha = 0.5) +
    geom_density(aes(x = loess, color = "LOESS"), alpha = 0.5) +
    labs(x = "Methylation Rate",
         y = "Density",
         title = "Distribution of Methylation Rates",
         color = "Method") +
    theme_minimal()

  # Print plots side by side
  gridExtra::grid.arrange(p1, p2, ncol = 2)

  # Print summary statistics
  cat("\nSummary Statistics:\n")
  print(summary(plot_data))

  # Print correlation coefficients
  cat("\nCorrelation with original rates:\n")
  cat("Quantile:", cor(plot_data$original, plot_data$quantile), "\n")
  cat("LOESS:", cor(plot_data$original, plot_data$loess), "\n")
}
