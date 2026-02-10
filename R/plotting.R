#' @title Visualization Functions for Methylation Data
#' @description
#' Functions for visualizing methylation data at various stages of the
#' normalization workflow. Includes density plots, QQ plots, MA plots,
#' and PCA/MDS plots.
#' @name plotting
NULL

utils::globalVariables(c(
  ".data", "PC1", "PC2", "Dim1", "Dim2", "M", "A",
  "sample1_rate", "sample2_rate"
))


# =========================================================================
# Exported Plot Functions
# =========================================================================

#' Density Plot of Methylation Rates
#'
#' Creates density plots of methylation rate distributions. When multiple
#' normalization stages are provided, each stage is shown in a separate
#' facet for easy comparison.
#'
#' @param meth_data A nested list of methylation data (group -> sample ->
#'   data.table), as returned by \code{\link{split_by_groups}} or the
#'   normalization functions. Can also be a named list of such objects
#'   for multi-stage comparison.
#' @param stages Optional character vector of stage names. Required when
#'   \code{meth_data} is a list of normalization stages. Ignored for
#'   single-stage data.
#' @param sample_size Number of sites to randomly sample for plotting.
#'   Default is 100000. Use \code{NULL} for all sites.
#' @param theme A ggplot2 theme. Default is \code{\link[ggplot2]{theme_minimal}}.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Single stage
#' plot_density(split_data)
#'
#' # Compare normalization stages
#' plot_density(
#'   list(raw = raw_data, normalized = norm_data),
#'   stages = c("raw", "normalized")
#' )
#' }
#'
#' @export
plot_density <- function(meth_data,
                         stages = NULL,
                         sample_size = 100000,
                         theme = ggplot2::theme_minimal()) {
  
  prepared <- .prepare_plot_data(meth_data, stages, sample_size)
  data <- prepared$data
  stages <- prepared$stages
  
  data$stage <- factor(data$stage, levels = stages)
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = rate, color = sample)) +
    ggplot2::geom_density(linewidth = 0.7) +
    ggplot2::facet_wrap(~stage) +
    ggplot2::labs(
      title = "Methylation Rate Distributions",
      x = "Methylation Rate",
      y = "Density"
    ) +
    theme
  
  return(p)
}


#' QQ Plot for Methylation Rate Comparison
#'
#' Creates quantile-quantile plots comparing methylation rates between
#' pairs of samples within each group. Uses hexagonal binning for large
#' datasets. When multiple stages are provided, plots are generated for
#' each stage.
#'
#' @inheritParams plot_density
#' @param bins Number of hexagonal bins for the hex plot. Default is 50.
#'
#' @return A ggplot object (or patchwork composite if multiple pairs exist).
#'
#' @examples
#' \dontrun{
#' plot_qq(split_data)
#' }
#'
#' @export
plot_qq <- function(meth_data,
                    stages = NULL,
                    sample_size = 100000,
                    bins = 50,
                    theme = ggplot2::theme_minimal()) {
  
  prepared <- .prepare_plot_data(meth_data, stages, sample_size)
  data <- prepared$data
  stages <- prepared$stages
  
  # Try hexbin first, fall back to simple
  hex_plots <- tryCatch(
    .create_qq_hexbin_plot(data, stages, bins, theme),
    error = function(e) {
      message("Hexbin QQ plot failed: ", e$message, ". Using simple QQ plot.")
      NULL
    }
  )
  
  if (is.null(hex_plots) || length(hex_plots) == 0) {
    return(.create_simple_qq_plot(data, stages, theme))
  }
  
  if (length(hex_plots) == 1) {
    return(hex_plots[[1]])
  }
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    max_plots <- min(4, length(hex_plots))
    return(patchwork::wrap_plots(hex_plots[seq_len(max_plots)], ncol = 2))
  }
  
  return(hex_plots[[1]])
}


#' MA Plot for Methylation Data
#'
#' Creates MA (M vs A) plots comparing methylation rates between pairs
#' of samples. The M-axis shows the difference in rates between two
#' samples, and the A-axis shows the average rate. Deviations from
#' M = 0 indicate systematic biases. When multiple stages are provided,
#' a faceted plot allows comparison of bias before and after normalization.
#'
#' @inheritParams plot_density
#' @param bins Number of hexagonal bins. Default is 60.
#' @param pair Optional character vector of length 2 specifying two sample
#'   names to compare. If \code{NULL} (default), the first two samples
#'   within the first group are used.
#'
#' @return A ggplot object (or patchwork composite for multiple pairs).
#'
#' @details
#' For each pair of samples, M and A are defined as:
#' \itemize{
#'   \item \code{M = rate_sample1 - rate_sample2} (difference)
#'   \item \code{A = (rate_sample1 + rate_sample2) / 2} (mean)
#' }
#'
#' A well-normalized dataset should show points centered around M = 0
#' across the full range of A.
#'
#' @examples
#' \dontrun{
#' # Compare first two samples automatically
#' plot_ma(list(raw = raw_data, normalized = norm_data),
#'         stages = c("raw", "normalized"))
#'
#' # Specify a sample pair
#' plot_ma(split_data, pair = c("WT_1", "WT_2"))
#' }
#'
#' @export
plot_ma <- function(meth_data,
                    stages = NULL,
                    sample_size = 100000,
                    bins = 60,
                    pair = NULL,
                    theme = ggplot2::theme_minimal()) {
  
  prepared <- .prepare_plot_data(meth_data, stages, sample_size)
  data <- prepared$data
  stages <- prepared$stages
  
  plot_list <- list()
  
  for (stage_name in stages) {
    stage_data <- data[data$stage == stage_name, ]
    
    # Determine which pair to plot
    if (!is.null(pair)) {
      s1 <- pair[1]
      s2 <- pair[2]
    } else {
      # Use first group, first two samples
      first_group <- unique(stage_data$group)[1]
      group_samples <- unique(stage_data$sample[stage_data$group == first_group])
      
      if (length(group_samples) < 2) {
        message("Fewer than 2 samples in group '", first_group,
                "' for stage '", stage_name, "'. Skipping.")
        next
      }
      s1 <- group_samples[1]
      s2 <- group_samples[2]
    }
    
    # Extract matched data
    d1 <- stage_data[stage_data$sample == s1, ]
    d2 <- stage_data[stage_data$sample == s2, ]
    common_rows <- intersect(d1$row_id, d2$row_id)
    
    if (length(common_rows) < 10) {
      message("Insufficient shared sites for MA plot: ", s1, " vs ", s2)
      next
    }
    
    r1 <- d1$rate[match(common_rows, d1$row_id)]
    r2 <- d2$rate[match(common_rows, d2$row_id)]
    
    ma_df <- data.frame(
      A = (r1 + r2) / 2,
      M = r1 - r2
    )
    ma_df <- stats::na.omit(ma_df)
    
    # Skip if no variation
    if (diff(range(ma_df$A)) <= 0 || diff(range(ma_df$M)) <= 0) {
      message("No variation in MA data for stage '", stage_name, "'. Skipping.")
      next
    }
    
    p <- ggplot2::ggplot(ma_df, ggplot2::aes(x = A, y = M)) +
      ggplot2::geom_hex(bins = bins) +
      ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                           color = "blue", se = FALSE, linewidth = 0.8) +
      ggplot2::scale_fill_viridis_c(option = "magma") +
      ggplot2::labs(
        title = stage_name,
        subtitle = paste(s1, "vs", s2),
        x = "A (Mean Rate)",
        y = "M (Rate Difference)"
      ) +
      ggplot2::coord_cartesian(ylim = c(
        max(-1, stats::quantile(ma_df$M, 0.001, na.rm = TRUE) * 1.2),
        min(1, stats::quantile(ma_df$M, 0.999, na.rm = TRUE) * 1.2)
      )) +
      theme
    
    plot_list[[stage_name]] <- p
  }
  
  if (length(plot_list) == 0) {
    stop("Could not create any MA plots with the provided data.")
  }
  
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  }
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(plot_list, nrow = 1))
  }
  
  return(plot_list[[1]])
}


#' PCA Plot of Methylation Data
#'
#' Creates a PCA (principal component analysis) or MDS (multidimensional
#' scaling) plot of methylation rate data. Samples are projected onto the
#' first two components and colored by group. This is useful for assessing
#' whether normalization improves replicate clustering while preserving
#' biological group separation.
#'
#' @inheritParams plot_density
#' @param method Character string: \code{"pca"} (default) or \code{"mds"}.
#' @param n_sites Number of most variable sites to use for dimensionality
#'   reduction. Default is 5000. Using too many sites can obscure structure;
#'   using too few reduces stability.
#' @param labels Logical. If \code{TRUE} (default), sample names are shown
#'   on the plot.
#' @param point_size Numeric size for plot points. Default is 3.
#'
#' @return A ggplot object. When multiple stages are provided, returns a
#'   faceted or patchwork composite plot.
#'
#' @details
#' The function selects the top \code{n_sites} most variable sites across
#' all samples, constructs a site-by-sample matrix, and performs PCA or
#' classical MDS. Sites with any missing values are excluded.
#'
#' For multi-stage comparison, each stage is reduced independently so that
#' the axes reflect the variance structure of that specific stage.
#'
#' @examples
#' \dontrun{
#' # Single stage PCA
#' plot_pca(split_data)
#'
#' # Compare stages
#' plot_pca(
#'   list(raw = raw_data, normalized = norm_data),
#'   stages = c("raw", "normalized")
#' )
#'
#' # Use MDS instead
#' plot_pca(split_data, method = "mds")
#' }
#'
#' @export
plot_pca <- function(meth_data,
                     stages = NULL,
                     method = c("pca", "mds"),
                     n_sites = 5000,
                     labels = TRUE,
                     point_size = 3,
                     theme = ggplot2::theme_minimal()) {
  
  method <- match.arg(method)
  
  # Detect whether input is multi-stage or single-stage
  stage_list <- .detect_stages(meth_data, stages)
  
  plot_list <- list()
  
  for (stage_name in names(stage_list)) {
    stage_data <- stage_list[[stage_name]]
    
    result <- tryCatch(
      .compute_reduction(stage_data, method, n_sites),
      error = function(e) {
        message("Dimensionality reduction failed for '", stage_name, "': ", e$message)
        NULL
      }
    )
    
    if (is.null(result)) next
    
    coords <- result$coords
    var_explained <- result$var_explained
    
    # Build axis labels
    if (method == "pca") {
      x_lab <- sprintf("PC1 (%.1f%%)", var_explained[1])
      y_lab <- sprintf("PC2 (%.1f%%)", var_explained[2])
    } else {
      x_lab <- "Dimension 1"
      y_lab <- "Dimension 2"
    }
    
    p <- ggplot2::ggplot(coords, ggplot2::aes(
      x = .data[["dim1"]], y = .data[["dim2"]],
      color = .data[["group"]]
    )) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::labs(
        title = stage_name,
        x = x_lab,
        y = y_lab,
        color = "Group"
      ) +
      theme
    
    if (labels) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data[["sample"]]),
        vjust = -0.8, hjust = 0.5, size = 3, show.legend = FALSE
      )
    }
    
    plot_list[[stage_name]] <- p
  }
  
  if (length(plot_list) == 0) {
    stop("Could not create any PCA/MDS plots with the provided data.")
  }
  
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  }
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(plot_list, nrow = 1))
  }
  
  return(plot_list[[1]])
}


# =========================================================================
# Internal Helper Functions
# =========================================================================

#' Detect whether input is single-stage or multi-stage and return a
#' named list of stage data
#' @return Named list where each element is a nested group/sample list
#' @keywords internal
#' @noRd
.detect_stages <- function(meth_data, stages = NULL) {
  # Check if the input is a list of stages (list of list of lists)
  # vs a single stage (list of lists of data.tables)
  if (.is_single_stage(meth_data)) {
    stage_name <- if (!is.null(stages)) stages[1] else "Data"
    return(stats::setNames(list(meth_data), stage_name))
  }
  
  # Multi-stage: meth_data is a list of nested lists
  if (is.null(stages)) {
    stages <- if (!is.null(names(meth_data))) names(meth_data) else
      paste("Stage", seq_along(meth_data))
  }
  
  if (length(stages) != length(meth_data)) {
    stop("Length of 'stages' must match length of 'meth_data'.")
  }
  
  stats::setNames(meth_data, stages)
}


#' Check if data is a single normalization stage
#' @return Logical
#' @keywords internal
#' @noRd
.is_single_stage <- function(meth_data) {
  if (!is.list(meth_data) || length(meth_data) == 0) return(FALSE)
  
  # A single stage is a list of lists where the innermost elements
  # are data.frames/data.tables (group -> sample -> dt)
  first <- meth_data[[1]]
  
  if (!is.list(first)) return(FALSE)
  
  # Check if the first element of the first group is a data.frame
  first_inner <- first[[1]]
  return(is.data.frame(first_inner) || data.table::is.data.table(first_inner))
}


#' Prepare data for plotting by flattening and optionally sampling
#'
#' Takes either a single-stage or multi-stage nested list and returns
#' a flat data.frame suitable for ggplot2.
#'
#' @param meth_data Nested list or list of nested lists
#' @param stages Stage names (can be NULL)
#' @param sample_size Number of sites to sample (NULL for all)
#' @return List with \code{data} (data.frame) and \code{stages} (character vector)
#' @keywords internal
#' @noRd
.prepare_plot_data <- function(meth_data, stages = NULL, sample_size = NULL) {
  stage_list <- .detect_stages(meth_data, stages)
  stages <- names(stage_list)
  
  all_data <- .flatten_stages(stage_list, stages)
  
  # Sample if requested
  if (!is.null(sample_size) && nrow(all_data) > 0 && sample_size < nrow(all_data)) {
    unique_ids <- unique(all_data$row_id)
    set.seed(42)
    sampled_ids <- sample(unique_ids, min(sample_size, length(unique_ids)))
    all_data <- all_data[all_data$row_id %in% sampled_ids, ]
  }
  
  all_data$stage <- factor(all_data$stage, levels = stages)
  
  list(data = all_data, stages = stages)
}


#' Flatten a named list of stages into a long data.frame
#' @return data.frame with columns: stage, group, sample, rate, cov, row_id
#' @keywords internal
#' @noRd
.flatten_stages <- function(stage_list, stages) {
  result_parts <- list()
  idx <- 1
  
  for (stage_name in stages) {
    stage_data <- stage_list[[stage_name]]
    
    for (group_name in names(stage_data)) {
      group_data <- stage_data[[group_name]]
      
      for (sample_name in names(group_data)) {
        df <- group_data[[sample_name]]
        
        # Ensure rate column exists
        if (!"rate" %in% names(df)) {
          if (all(c("mc", "cov") %in% names(df))) {
            df$rate <- ifelse(df$cov > 0, df$mc / df$cov, 0)
          } else {
            warning(sprintf("Skipping %s/%s/%s: no rate column.",
                            stage_name, group_name, sample_name))
            next
          }
        }
        
        result_parts[[idx]] <- data.frame(
          stage   = stage_name,
          group   = group_name,
          sample  = sample_name,
          rate    = df$rate,
          cov     = if ("cov" %in% names(df)) df$cov else NA_real_,
          row_id  = seq_along(df$rate),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
  }
  
  if (length(result_parts) == 0) {
    stop("No valid sample data found. Check your data structure.")
  }
  
  do.call(rbind, result_parts)
}


#' Compute PCA or MDS coordinates from nested methylation data
#'
#' @param stage_data Nested group/sample list for a single stage
#' @param method "pca" or "mds"
#' @param n_sites Number of most-variable sites to retain
#' @return List with \code{coords} (data.frame) and \code{var_explained} (numeric)
#' @keywords internal
#' @noRd
.compute_reduction <- function(stage_data, method, n_sites) {
  # Build a site-by-sample matrix
  sample_info <- list()
  rate_columns <- list()
  idx <- 1
  
  for (group_name in names(stage_data)) {
    group_data <- stage_data[[group_name]]
    for (sample_name in names(group_data)) {
      df <- group_data[[sample_name]]
      
      # Ensure rate exists
      if (!"rate" %in% names(df)) {
        if (all(c("mc", "cov") %in% names(df))) {
          df$rate <- ifelse(df$cov > 0, df$mc / df$cov, 0)
        } else {
          warning("Skipping sample ", sample_name, ": no rate column")
          next
        }
      }
      
      rate_columns[[idx]] <- df$rate
      sample_info[[idx]] <- data.frame(
        sample = sample_name,
        group  = group_name,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  if (length(rate_columns) < 3) {
    stop("At least 3 samples are required for PCA/MDS.")
  }
  
  # Verify all samples have the same number of sites
  site_counts <- vapply(rate_columns, length, integer(1))
  if (length(unique(site_counts)) != 1) {
    stop("All samples must have the same number of sites for PCA/MDS. ",
         "Run find_shared_sites() first.")
  }
  
  # Combine into matrix (sites x samples)
  rate_mat <- do.call(cbind, rate_columns)
  info_df <- do.call(rbind, sample_info)
  
  # Remove sites with any NA
  complete <- stats::complete.cases(rate_mat)
  rate_mat <- rate_mat[complete, , drop = FALSE]
  
  if (nrow(rate_mat) < 100) {
    stop("Fewer than 100 complete sites remaining. Check for missing data.")
  }
  
  # Select most variable sites
  site_vars <- apply(rate_mat, 1, stats::var)
  top_idx <- order(site_vars, decreasing = TRUE)[seq_len(min(n_sites, nrow(rate_mat)))]
  rate_mat <- rate_mat[top_idx, , drop = FALSE]
  
  var_explained <- c(NA, NA)
  
  if (method == "pca") {
    # PCA on transposed matrix (samples as rows)
    pca_result <- stats::prcomp(t(rate_mat), center = TRUE, scale. = FALSE)
    coords <- data.frame(
      dim1   = pca_result$x[, 1],
      dim2   = pca_result$x[, 2],
      sample = info_df$sample,
      group  = info_df$group,
      stringsAsFactors = FALSE
    )
    total_var <- sum(pca_result$sdev^2)
    var_explained <- (pca_result$sdev[1:2]^2 / total_var) * 100
  } else {
    # Classical MDS
    dist_mat <- stats::dist(t(rate_mat))
    mds_result <- stats::cmdscale(dist_mat, k = 2, eig = TRUE)
    coords <- data.frame(
      dim1   = mds_result$points[, 1],
      dim2   = mds_result$points[, 2],
      sample = info_df$sample,
      group  = info_df$group,
      stringsAsFactors = FALSE
    )
    if (!is.null(mds_result$eig)) {
      pos_eig <- mds_result$eig[mds_result$eig > 0]
      var_explained <- (pos_eig[1:2] / sum(pos_eig)) * 100
    }
  }
  
  list(coords = coords, var_explained = var_explained)
}


# --- QQ plot internals (kept from original) ------------------------------

#' Create hexbin QQ plots for all within-group sample pairs
#' @keywords internal
#' @noRd
.create_qq_hexbin_plot <- function(data, stages, bins = 50,
                                   theme = ggplot2::theme_minimal()) {
  plot_list <- list()
  
  for (stage_name in levels(data$stage)) {
    stage_data <- data[data$stage == stage_name, ]
    
    for (group_name in unique(stage_data$group)) {
      group_data <- stage_data[stage_data$group == group_name, ]
      samples <- unique(group_data$sample)
      
      if (length(samples) < 2) next
      
      pairs <- utils::combn(samples, 2, simplify = FALSE)
      
      for (pair in pairs) {
        d1 <- group_data[group_data$sample == pair[1], ]
        d2 <- group_data[group_data$sample == pair[2], ]
        
        common_rows <- intersect(d1$row_id, d2$row_id)
        if (length(common_rows) < 10) next
        
        plot_df <- data.frame(
          x = d1$rate[match(common_rows, d1$row_id)],
          y = d2$rate[match(common_rows, d2$row_id)]
        )
        plot_df <- stats::na.omit(plot_df)
        
        if (nrow(plot_df) < 10) next
        if (diff(range(plot_df$x)) <= 0 || diff(range(plot_df$y)) <= 0) next
        
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_hex(bins = bins) +
          ggplot2::geom_abline(color = "red", linetype = "dashed") +
          ggplot2::scale_fill_viridis_c(option = "magma") +
          ggplot2::labs(
            title = paste(group_name, "-", stage_name),
            subtitle = paste(pair, collapse = " vs "),
            x = pair[1],
            y = pair[2]
          ) +
          theme
        
        plot_id <- paste(stage_name, group_name,
                         paste(pair, collapse = "_"), sep = "_")
        plot_list[[plot_id]] <- p
      }
    }
  }
  
  if (length(plot_list) == 0) return(NULL)
  return(plot_list)
}


#' Simple QQ plot fallback when hexbin fails
#' @keywords internal
#' @noRd
.create_simple_qq_plot <- function(data, stages,
                                   theme = ggplot2::theme_minimal()) {
  avg_data <- stats::aggregate(rate ~ group + stage, data, mean)
  
  if (length(stages) >= 2) {
    stage1 <- stages[1]
    stage2 <- stages[2]
    comp_data <- reshape2::dcast(avg_data, group ~ stage, value.var = "rate")
    
    p <- ggplot2::ggplot(comp_data,
                         ggplot2::aes(x = .data[[stage1]], y = .data[[stage2]])) +
      ggplot2::geom_point(ggplot2::aes(color = group), size = 3) +
      ggplot2::geom_abline(color = "red", linetype = "dashed") +
      ggplot2::labs(
        title = paste("Comparison:", stage1, "vs", stage2),
        x = paste("Rate in", stage1),
        y = paste("Rate in", stage2)
      ) +
      theme
  } else {
    p <- ggplot2::ggplot(avg_data,
                         ggplot2::aes(x = group, y = rate, fill = group)) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~stage) +
      ggplot2::labs(
        title = "Average Methylation Rate by Group",
        x = "Group",
        y = "Average Rate"
      ) +
      theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}


# --- Report layout and saving helpers -----------------------------------

#' Determine optimal layout for a set of plots
#' @keywords internal
#' @noRd
.get_optimal_layout <- function(plot_types) {
  n <- length(plot_types)
  if (n <= 2) return(list(ncol = n, nrow = 1))
  if (n <= 4) return(list(ncol = 2, nrow = ceiling(n / 2)))
  ncol <- round(sqrt(n))
  list(ncol = ncol, nrow = ceiling(n / ncol))
}


#' Save a single plot to a PDF report with a title page
#' @keywords internal
#' @noRd
.save_plot_report <- function(plot, plot_type, output_dir, n_facets) {
  report_file <- file.path(output_dir, paste0(plot_type, "_report.pdf"))
  tryCatch({
    grDevices::pdf(report_file, width = 10, height = 8)
    grid::grid.newpage()
    grid::grid.text(
      paste0(toupper(substr(plot_type, 1, 1)),
             substr(plot_type, 2, nchar(plot_type)),
             " Plot Report"),
      x = 0.5, y = 0.8,
      gp = grid::gpar(fontsize = 24, fontface = "bold")
    )
    grid::grid.text(
      paste("Generated on:", Sys.Date()),
      x = 0.5, y = 0.6,
      gp = grid::gpar(fontsize = 16)
    )
    print(plot)
    grDevices::dev.off()
    message(plot_type, " report saved to ", report_file)
    return(report_file)
  }, error = function(e) {
    message("Error saving ", plot_type, " report: ", e$message)
    return(NULL)
  })
}