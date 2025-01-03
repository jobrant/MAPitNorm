#' Create data frame for plotting
#' @param data_list List of methylation data
#' @param stage Label for normalization stage
#' @return Data frame with methylation rates and metadata
#' @keywords internal
create_plot_df <- function(data_list, stage) {
  df <- do.call(rbind, lapply(names(data_list), function(group) {
    do.call(rbind, lapply(names(data_list[[group]]), function(sample) {
      data.frame(
        rate = data_list[[group]][[sample]]$rate,
        group = group,
        sample = sample,
        stage = stage
      )
    }))
  }))

  # Make stage a factor with ordered levels
  df$stage <- factor(df$stage,
                     levels = c("Original",
                                "Within-group normalized",
                                "Fully normalized"))
  return(df)
}

#' Visualize Normalization Effects
#'
#' @param original_data Original data list before normalization
#' @param within_norm Data after within-group normalization
#' @param full_norm Data after both normalizations
#' @param output_file Optional file path to save plots
#' @return A grid arrangement of four diagnostic plots
#' @importFrom ggplot2 ggplot aes geom_density geom_point geom_boxplot facet_wrap
#' @importFrom ggplot2 theme_minimal labs geom_abline geom_hline theme element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices dev.off
#' @export
visualize_normalization <- function(original_data, within_norm, full_norm, output_file = NULL) {
  p1 <- plot_rate_distributions(original_data, within_norm, full_norm)
  p2 <- plot_qq_comparison(original_data, within_norm, full_norm)
  p3 <- plot_ma(original_data, within_norm, full_norm)
  p4 <- plot_boxplots(original_data, within_norm, full_norm)

  combined <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 12, height = 10)
  }

  return(combined)
}

#' Function to visualize density plots of rate distribution
#'
#' @param original_data Original data list
#' @param within_norm Data after within-group normalization
#' @param full_norm Data after both normalizations
#' @return A ggplot object showing density distributions
#' @export
plot_rate_distributions <- function(original_data, within_norm, full_norm) {
  # Combine data with labels
  plot_data <- rbind(
    create_plot_df(original_data, "Original"),
    create_plot_df(within_norm, "Within-group normalized"),
    create_plot_df(full_norm, "Fully normalized")
  )

  ggplot(plot_data, aes(x = rate, color = sample)) +
    geom_density() +
    facet_wrap(~stage) +
    theme_minimal() +
    labs(title = "Distribution of methylation rates",
         x = "Methylation rate",
         y = "Density")
}


#' Create enhanced QQ comparison plots with density information
#' @param original_data Original data list
#' @param within_norm Within-group normalized data
#' @param full_norm Fully normalized data
#' @param color_scheme Color scheme for density plot ("magma", "viridis", "plasma", "inferno")
#' @param bins Number of hexagonal bins (default = 50)
#' @param sample_size Optional sampling size (NULL for full dataset)
#' @param seed Random seed for reproducibility (default = 42,
#'        the answer to life, the universe, and everything)
#' @return A ggplot object showing sample comparisons with density information
#' @export
plot_qq_comparison <- function(original_data, within_norm, full_norm,
                               color_scheme = "magma",
                               bins = 50,
                               sample_size = NULL,
                               seed = 42) {
  # Helper function to create QQ data
  create_qq_data <- function(data_list, stage) {
    df <- do.call(rbind, lapply(names(data_list), function(group) {
      samples <- names(data_list[[group]])
      if(length(samples) >= 2) {
        data.frame(
          sample1_rate = data_list[[group]][[samples[1]]]$rate,
          sample2_rate = data_list[[group]][[samples[2]]]$rate,
          group = group,
          stage = stage
        )
      }
    }))
    df$stage <- factor(df$stage,
                       levels = c("Original",
                                  "Within-group normalized",
                                  "Fully normalized"))
    return(df)
  }

  # Optional sampling
  if(!is.null(sample_size)) {
    sampled <- sample_for_visualization(original_data, within_norm, full_norm,
                                        n_sites = sample_size, seed = seed)
    original_data <- sampled$original
    within_norm <- sampled$within_norm
    full_norm <- sampled$full_norm
  }

  plot_data <- rbind(
    create_qq_data(original_data, "Original"),
    create_qq_data(within_norm, "Within-group normalized"),
    create_qq_data(full_norm, "Fully normalized")
  )

  ggplot(plot_data, aes(x = sample1_rate, y = sample2_rate)) +
    geom_hex(bins = bins) +
    geom_abline(color = "red", linetype = "dashed") +
    facet_grid(group ~ stage) +
    scale_fill_viridis_c(
      option = color_scheme,
      name = "Number of\ndata points",
      labels = scales::comma,
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 3.5,     # Increased width
        barheight = 10,   # Increased height
        title.hjust = 0.5,
        label.position = "right",
        ticks.linewidth = 1
      )
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.key.height = unit(2, "cm"),  # Explicitly set legend height
      legend.key.width = unit(1, "cm"),   # Explicitly set legend width
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    ) +
    labs(title = "Sample-to-sample comparison",
         subtitle = paste0("Density shown using ", bins, " bins"),
         x = "Sample 1 methylation rate",
         y = "Sample 2 methylation rate")
}


#' Function to visualize differences vs average (MA plot)
#' @param original_data Original data list before normalization
#' @param within_norm Data after within-group normalization
#' @param full_norm Data after both normalizations
#' @return A ggplot object showing MA plots
#' @export
plot_ma <- function(original_data, within_norm, full_norm,
                    max_points = 100000) {
create_ma_data <- function(data_list, stage) {
  df <- do.call(rbind, lapply(names(data_list), function(group) {
    samples <- names(data_list[[group]])
    if(length(samples) >= 2) {
      rates1 <- data_list[[group]][[samples[1]]]$rate
      rates2 <- data_list[[group]][[samples[2]]]$rate
      data.frame(
        mean_rate = (rates1 + rates2) / 2,
        rate_diff = rates1 - rates2,
        group = group,
        stage = stage
      )
    }
  }))
  df$stage <- factor(df$stage,
                     levels = c("Original",
                                "Within-group normalized",
                                "Fully normalized"))
  return(df)
}

  plot_data <- rbind(
    create_ma_data(original_data, "Original"),
    create_ma_data(within_norm, "Within-group normalized"),
    create_ma_data(full_norm, "Fully normalized")
  )

  ggplot(plot_data, aes(x = mean_rate, y = rate_diff)) +
    geom_hex(bins = 50) +
    geom_hline(yintercept = 0, color = "red") +
    facet_grid(group ~ stage) +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    labs(title = "MA plot: Difference vs Average",
         x = "Mean rate",
         y = "Rate difference")
}

#' Function to visualize distribution changes as Box plots
#' @param original_data Original data list before normalization
#' @param within_norm Data after within-group normalization
#' @param full_norm Data after both normalizations
#' @return A ggplot object showing boxplots
#' @export
plot_boxplots <- function(original_data, within_norm, full_norm) {
  create_box_data <- function(data_list, stage) {
    df <- do.call(rbind, lapply(names(data_list), function(group) {
      samples <- names(data_list[[group]])
      if(length(samples) >= 2) {
        rates1 <- data_list[[group]][[samples[1]]]$rate
        rates2 <- data_list[[group]][[samples[2]]]$rate
        data.frame(
          mean_rate = (rates1 + rates2) / 2,
          rate_diff = rates1 - rates2,
          group = group,
          stage = stage
        )
      }
    }))
    df$stage <- factor(df$stage,
                       levels = c("Original",
                                  "Within-group normalized",
                                  "Fully normalized"))
    return(df)
  }
  plot_data <- rbind(
    create_plot_df(original_data, "Original"),
    create_plot_df(within_norm, "Within-group normalized"),
    create_plot_df(full_norm, "Fully normalized")
  )

  ggplot(plot_data, aes(x = sample, y = rate, fill = group)) +
    geom_boxplot() +
    facet_wrap(~stage) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Distribution of rates by sample",
         x = "Sample",
         y = "Rate")
}
