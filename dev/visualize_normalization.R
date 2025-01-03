#' Visualize Normalization Effects
#'
#' Creates diagnostic plots to visualize the effects of both coverage and methylation rate normalization.
#'
#' @param original_data Original data list before normalization
#' @param coverage_normalized Data after coverage normalization
#' @param fully_normalized Data after both coverage and methylation rate normalization
#' @param group_names Character vector of group names
#' @param n_quantiles Number of quantiles for visualization (default = 10)
#' @param save_plot Logical indicating whether to save plots to file (default = FALSE)
#' @param output_dir Directory to save plots (default = current directory)
#' @param sample_size Number of sites to sample for visualization (default = 100000)
#' @importFrom ggplot2 ggplot aes geom_density facet_wrap scale_x_log10 labs
#'             theme_minimal geom_line geom_boxplot element_text theme
#' @importFrom data.table rbindlist data.table setDT :=
#' @importFrom gridExtra grid.arrange
#' @importFrom stats quantile
#' @importFrom grDevices pdf dev.off
#'
#' @export
#'
visualize_normalization <- function(original_data,
                                    coverage_normalized,
                                    fully_normalized,
                                    group_names,
                                    n_quantiles = 10,
                                    sample_size = 100000,
                                    save_plot = FALSE,
                                    output_dir = ".") {

  cat("Preparing data for visualization...\n")

  # Function to sample data efficiently
  sample_data <- function(data_list, stage, size) {
    # Use same indices for all samples for consistency
    set.seed(42)
    sample_idx <- sample(nrow(data_list[[1]][[1]]), size)

    # Convert to data.table and bind
    dt <- rbindlist(
      lapply(names(data_list), function(group) {
        rbindlist(
          lapply(data_list[[group]], function(df) {
            data.table(
              group = group,
              coverage = df$cov[sample_idx],
              methylation = df$mc[sample_idx],
              rate = df$rate[sample_idx],
              stage = stage
            )
          })
        )
      })
    )
    return(dt)
  }

  cat("Sampling data...\n")
  # Combine sampled data from all stages
  plot_data <- rbindlist(list(
    sample_data(original_data, "Original", sample_size),
    sample_data(coverage_normalized, "Coverage Normalized", sample_size),
    sample_data(fully_normalized, "Fully Normalized", sample_size)
  ))

  cat("Creating plots...\n")

  # 1. Coverage Distribution Plot
  p1 <- ggplot(plot_data, aes(x = coverage, color = group)) +
    geom_density() +
    facet_wrap(~stage) +
    scale_x_log10() +
    labs(title = "Coverage Distribution",
         x = "Coverage (log scale)",
         y = "Density") +
    theme_minimal()

  # 2. Methylation Rate Distribution
  p2 <- ggplot(plot_data, aes(x = rate, color = group)) +
    geom_density() +
    facet_wrap(~stage) +
    labs(title = "Methylation Rate Distribution",
         x = "Methylation Rate",
         y = "Density") +
    theme_minimal()

  # 3. Quantile Plot
  cat("Calculating quantiles...\n")
  quantile_data <- plot_data[, .(
    quantiles = list(quantile(rate, probs = seq(0, 1, length.out = n_quantiles)))
  ), by = .(group, stage)]

  quantile_data <- quantile_data[, .(
    quantile = seq_along(unlist(quantiles)),
    value = unlist(quantiles)
  ), by = .(group, stage)]

  p3 <- ggplot(quantile_data, aes(x = quantile, y = value, color = group)) +
    geom_line() +
    facet_wrap(~stage) +
    labs(title = "Quantile Plot of Methylation Rates",
         x = "Quantile",
         y = "Methylation Rate") +
    theme_minimal()

  # 4. Box plots of rates by group
  p4 <- ggplot(plot_data, aes(x = group, y = rate, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~stage) +
    labs(title = "Distribution of Methylation Rates by Group",
         x = "Group",
         y = "Methylation Rate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  cat("Combining plots...\n")
  combined_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

  if(save_plot) {
    filename <- file.path(output_dir, "normalization_diagnostics.pdf")
    ggsave(filename, combined_plot, width = 15, height = 12)
    cat(sprintf("Plots saved to %s\n", filename))
  }

  return(combined_plot)
}
