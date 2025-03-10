#' Create data frame for plotting
#' @param data_list List of methylation data
#' @param stage Label for normalization stage
#' @param group_name Name of the group (optional)
#' @return Data frame with methylation rates and metadata
#' @keywords internal
create_plot_df <- function(data_list, stage, group_name = NULL) {
  df <- do.call(rbind, lapply(names(data_list), function(sample) {
    data.frame(
      rate = data_list[[sample]]$rate,
      sample = sample,
      stage = stage,
      group = if(is.null(group_name)) "group" else group_name
    )
  }))

  # Make stage a factor with ordered levels
  df$stage <- factor(df$stage,
                     levels = c("Original",
                                "Within-group normalized",
                                "Fully normalized"))
  return(df)
}

#' Create plot data frame for groups
#' @param data_list List of methylation data
#' @param stage Label for normalization stage
#' @param group_names Vector of group names
#' @return Data frame with methylation rates and group information
#' @keywords internal
create_plot_df_groups <- function(data_list, stage, group_names) {
  do.call(rbind, lapply(seq_along(group_names), function(i) {
    group <- group_names[i]
    do.call(rbind, lapply(names(data_list[[i]]), function(sample) {
      data.frame(
        rate = data_list[[i]][[sample]]$rate,
        sample = sample,
        group = group,
        stage = stage
      )
    }))
  }))
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


#' Create helper function for single group QQ comparison
#' @param original_data Original data list
#' @param within_norm Within-group normalized data
#' @param full_norm Fully normalized data
#' @param color_scheme Color scheme for density plot ("magma", "viridis", "plasma", "inferno")
#' @param bins Number of hexagonal bins (default = 50)
#' @return A ggplot object showing sample comparisons with density information
#' @keywords internal
plot_qq_comparison_single_group <- function(original, normalized, group_name,
                                            color_scheme = "magma",
                                            bins = 50) {
  # Get all pairwise comparisons for this group
  create_qq_data_single <- function(data, stage) {
    samples <- names(data)
    pairs <- combn(samples, 2, simplify = FALSE)
    do.call(rbind, lapply(pairs, function(pair) {
      data.frame(
        sample1_rate = data[[pair[1]]]$rate,
        sample2_rate = data[[pair[2]]]$rate,
        pair = paste(pair[1], "vs", pair[2]),
        stage = stage
      )
    }))
  }

  plot_data <- rbind(
    create_qq_data_single(original, "Original"),
    create_qq_data_single(normalized, "Normalized")
  )

  # Set the order explicitly
  plot_data$stage <- factor(plot_data$stage,
                            levels = c("Original",
                                       "Within-group normalized"))

  ggplot(plot_data, aes(x = sample1_rate, y = sample2_rate)) +
    geom_hex(bins = bins) +
    geom_abline(color = "red", linetype = "dashed") +
    facet_grid(pair ~ stage) +
    scale_fill_viridis_c(
      option = color_scheme,
      name = "Number of\ndata points",
      labels = scales::comma,
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 1.5,
        barheight = 10
      )
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.height = unit(2, "cm"),
      legend.key.width = unit(1, "cm"),
      strip.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = paste("Sample comparisons -", group_name),
         x = "Sample 1 methylation rate",
         y = "Sample 2 methylation rate")
}

#' Create within-group normalization report
#' @param original_data Original data list
#' @param normalized_data After within-group normalization
#' @param group_names Vector of group names
#' @param output_file Optional file path to save report
#' @param color_scheme Color scheme for density plots
#' @param bins Number of hexagonal bins
#' @return A list containing both plot arrangements
#' @export
plot_within_group_report <- function(original_data, normalized_data,
                                     group_names,
                                     output_file = NULL,
                                     color_scheme = "magma",
                                     bins = 50) {

  # Create plots for each group
  plots <- lapply(seq_along(group_names), function(i) {
    group <- group_names[i]

    # QQ comparison plot
    p1 <- plot_qq_comparison_single_group(
      original_data[[i]],
      normalized_data[[i]],
      group_name = group,
      color_scheme = color_scheme,
      bins = bins
    )

    # Distribution plot
    plot_data <- rbind(
      create_plot_df(original_data[[i]], "Original", group),
      create_plot_df(normalized_data[[i]], "Within-group normalized", group)
    )

    # Set the order explicitly
    plot_data$stage <- factor(plot_data$stage,
                              levels = c("Original",
                                         "Within-group normalized"))

    p2 <- ggplot(plot_data, aes(x = rate, color = sample)) +
      geom_density() +
      facet_wrap(~stage) +
      theme_minimal() +
      labs(title = paste("Rate distributions -", group),
           x = "Methylation rate",
           y = "Density")

    list(qq = p1, dist = p2)
  })

  # Separate QQ and distribution plots
  qq_plots <- lapply(plots, function(x) x$qq)
  dist_plots <- lapply(plots, function(x) x$dist)

  # Create two separate arrangements
  qq_page <- gridExtra::grid.arrange(
    grobs = qq_plots,
    ncol = 2,
    top = "Within-group normalization: Sample comparisons"
  )

  dist_page <- gridExtra::grid.arrange(
    grobs = dist_plots,
    ncol = 2,
    top = "Within-group normalization: Rate distributions"
  )

  if (!is.null(output_file)) {
    # Save as multi-page PDF
    pdf(output_file, width = 15, height = 5 * ceiling(length(group_names)/2))
    grid::grid.draw(qq_page)
    grid::grid.newpage()
    grid::grid.draw(dist_page)
    dev.off()
  }

  # Return both arrangements
  return(list(qq_page = qq_page, dist_page = dist_page))
}

#' Create distribution data for group comparisons
#' @param data_list List of methylation data
#' @param stage Label for normalization stage
#' @param group_names Vector of group names
#' @return Data frame with methylation rates and group information
#' @keywords internal
create_group_dist_data <- function(data_list, stage, group_names) {
  do.call(rbind, lapply(seq_along(group_names), function(i) {
    group <- group_names[i]
    # Keep individual sample rates
    do.call(rbind, lapply(names(data_list[[i]]), function(sample) {
      data.frame(
        rate = data_list[[i]][[sample]]$rate,
        sample = sample,
        group = group,
        stage = stage
      )
    }))
  }))
}

#' Create between-group normalization report
#' @param within_norm After within-group normalization
#' @param full_norm After both normalizations
#' @param group_names Vector of group names
#' @param output_file Optional file path to save report
#' @param color_scheme Color scheme for density plots
#' @param bins Number of hexagonal bins
#' @return A grid arrangement of diagnostic plots
#' @export
plot_between_group_report <- function(within_norm, full_norm,
                                      group_names,
                                      output_file = NULL,
                                      color_scheme = "magma",
                                      bins = 50) {

  # Calculate group means for comparison
  calculate_group_means <- function(data_list, groups) {
    means_list <- lapply(seq_along(groups), function(i) {
      rates <- do.call(cbind, lapply(data_list[[i]], function(df) df$rate))
      rowMeans(rates)
    })
    names(means_list) <- groups
    means_list
  }

  # Create plots
  # 1. Group mean distributions
  plot_data_dist <- rbind(
    create_group_dist_data(within_norm, "Before between-group norm", group_names),
    create_group_dist_data(full_norm, "After between-group norm", group_names)
  )

  plot_data_dist$stage <- factor(plot_data_dist$stage,
                                 levels = c("Before between-group norm",
                                            "After between-group norm"))

  p1 <- ggplot(plot_data_dist, aes(x = rate, color = group)) +
    geom_density() +
    facet_wrap(~stage) +
    theme_minimal() +
    labs(title = "Distribution of methylation rates by group",
         x = "Methylation rate",
         y = "Density")

  # 2. Group means comparison
  means_before <- calculate_group_means(within_norm, group_names)
  means_after <- calculate_group_means(full_norm, group_names)

  plot_data_means <- data.frame(
    mean_before = unlist(means_before),
    mean_after = unlist(means_after),
    group = rep(group_names, each = length(means_before[[1]]))
  )

  p2 <- ggplot(plot_data_means, aes(x = mean_before, y = mean_after)) +
    geom_hex(bins = bins) +
    geom_abline(color = "red", linetype = "dashed") +
    scale_fill_viridis_c(option = color_scheme) +
    theme_minimal() +
    labs(title = "Group means before vs after normalization",
         x = "Mean rate before",
         y = "Mean rate after")

  # 3. Boxplots by group
  plot_data_box <- rbind(
    create_plot_df_groups(within_norm, "Before between-group norm", group_names),
    create_plot_df_groups(full_norm, "After between-group norm", group_names)
  )

  plot_data_box$stage <- factor(plot_data_box$stage,
                                levels = c("Before between-group norm",
                                           "After between-group norm"))

  p3 <- ggplot(plot_data_box, aes(x = group, y = rate, fill = group)) +
    geom_boxplot() +
    facet_wrap(~stage) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Distribution of rates by group",
         x = "Group",
         y = "Rate")

  # Combine plots
  combined <- gridExtra::grid.arrange(p1, p2, p3,
                                      ncol = 2,
                                      top = "Between-group normalization effects")

  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 15, height = 10)
  }

  return(combined)
}

