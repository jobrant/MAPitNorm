# File: R/create_plots.R

#' Create Plots
#'
#' Create plots for methylation data.
#'
#' @param data Prepared data for plotting.
#'
#' @param measure Measurement type (e.g., "rate", "cov", "mc").
#'
#' @param ylabel Label for the y-axis.
#'
#' @param normalized Boolean indicating if the data is normalized.
#'
#' @return List of ggplot objects.
#'
#' @import ggplot2
#'
#' @import dplyr
#'
#' @import tidyr
#'
#' @export
#'
create_plots <- function(data, measure, ylabel, normalized = TRUE) {
  data$Decile <- cut(data$rate, breaks = quantile(data$rate, probs = seq(0, 1, 0.1)),
                     include.lowest = TRUE, labels = FALSE)

  summary_data <- data %>%
    group_by(group, Decile) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

  diff_data <- summary_data %>%
    pivot_wider(names_from = group, values_from = mean_value) %>%
    pivot_longer(cols = -c(Decile, M1), names_to = "group", values_to = "value") %>%
    mutate(diff = M1 - value)

  norm_status <- if(normalized) "Normalized" else "Un-normalized"

  p1 <- ggplot(summary_data, aes(x = Decile, y = mean_value, color = group)) +
    geom_line() +
    geom_point() +
    labs(x = "Decile of Average Methylation Rate", y = ylabel,
         title = paste(norm_status, ylabel)) +
    theme_minimal()

  p2 <- ggplot(diff_data, aes(x = Decile, y = diff, color = group)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Decile of Average Methylation Rate", y = paste(ylabel, "Difference"),
         title = paste(norm_status, ylabel, "Difference")) +
    theme_minimal()

  list(p1, p2)
}
