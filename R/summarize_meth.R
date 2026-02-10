#' Summarize Methylation Data Object
#'
#' Prints a concise overview of a methylation data object as returned by
#' \code{\link{load_data}}, \code{\link{split_by_groups}}, or the
#' normalization functions.
#'
#' @param meth_data A list of methylation data (flat sample list or
#'   nested group/replicate structure).
#' @param show_site_types Logical. If \code{TRUE} (default), shows a breakdown
#'   of methylation site types (e.g., GCH, HCG) when a \code{site} column
#'   is present.
#'
#' @return Invisibly returns a list of summary statistics, including
#'   \code{n_samples}, \code{n_groups}, \code{n_sites}, and
#'   \code{per_sample} details. Called primarily for its side effect of
#'   printing to the console.
#'
#' @examples
#' \dontrun{
#' all_samples <- load_data(dir_path = "data/", sample_sheet = "samples.csv")
#' summarize_meth(all_samples)
#'
#' split_data <- split_by_groups(all_samples)
#' summarize_meth(split_data)
#' }
#'
#' @export
summarize_meth <- function(meth_data, show_site_types = TRUE) {
  if (!is.list(meth_data) || length(meth_data) == 0) {
    stop("'meth_data' must be a non-empty list.")
  }

  # Detect structure: flat (sample list) vs nested (group -> sample)
  is_nested <- all(vapply(meth_data, is.list, logical(1)))

  # Collect metadata if available
  sample_metadata <- attr(meth_data, "sample_metadata")

  # --- Header ---
  cat("MAPitNorm Methylation Data Summary\n")
  cat(strrep("=", 40), "\n")

  if (is_nested) {
    .summarize_nested(meth_data, sample_metadata, show_site_types)
  } else {
    .summarize_flat(meth_data, sample_metadata, show_site_types)
  }
}


#' Summarize a flat sample list
#' @keywords internal
#' @noRd
.summarize_flat <- function(meth_data, sample_metadata, show_site_types) {
  n_samples <- length(meth_data)
  cat("Structure:  flat sample list\n")
  cat("Samples:    ", n_samples, "\n")

  if (!is.null(sample_metadata)) {
    groups <- unique(sample_metadata$group_id)
    cat("Groups:     ", paste(groups, collapse = ", "), "\n")
  }

  cat(strrep("-", 40), "\n")

  # Per-sample summary
  stats <- .collect_sample_stats(meth_data)
  .print_sample_table(stats)

  # Site types
  if (show_site_types) {
    .print_site_types(meth_data[[1]])
  }

  cat(strrep("=", 40), "\n")

  return(invisible(list(
    structure  = "flat",
    n_samples  = n_samples,
    n_groups   = if (!is.null(sample_metadata)) length(unique(sample_metadata$group_id)) else NA,
    per_sample = stats
  )))
}


#' Summarize a nested group/sample list
#' @keywords internal
#' @noRd
.summarize_nested <- function(meth_data, sample_metadata, show_site_types) {
  group_names <- names(meth_data)
  n_groups <- length(group_names)
  n_total <- sum(vapply(meth_data, length, integer(1)))

  cat("Structure:  nested (group/sample)\n")
  cat("Groups:     ", n_groups, "\n")
  cat("Samples:    ", n_total, " total\n\n")

  all_stats <- list()

  for (group in group_names) {
    group_data <- meth_data[[group]]
    n_reps <- length(group_data)
    cat(sprintf("  Group '%s' (%d replicate%s)\n",
                group, n_reps, ifelse(n_reps == 1, "", "s")))

    stats <- .collect_sample_stats(group_data)
    .print_sample_table(stats, indent = "    ")
    all_stats[[group]] <- stats
    cat("\n")
  }

  # Site types from first available sample
  if (show_site_types) {
    first_sample <- meth_data[[1]][[1]]
    .print_site_types(first_sample)
  }

  cat(strrep("=", 40), "\n")

  return(invisible(list(
    structure  = "nested",
    n_groups   = n_groups,
    n_samples  = n_total,
    groups     = group_names,
    per_group  = all_stats
  )))
}


#' Collect per-sample statistics
#' @return data.frame with one row per sample
#' @keywords internal
#' @noRd
.collect_sample_stats <- function(sample_list) {
  stats <- data.frame(
    sample       = names(sample_list),
    n_sites      = integer(length(sample_list)),
    mean_cov     = numeric(length(sample_list)),
    median_cov   = numeric(length(sample_list)),
    mean_rate    = numeric(length(sample_list)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(sample_list)) {
    df <- sample_list[[i]]

    stats$n_sites[i] <- nrow(df)

    if ("cov" %in% names(df)) {
      stats$mean_cov[i]   <- round(mean(df$cov, na.rm = TRUE), 1)
      stats$median_cov[i] <- round(stats::median(df$cov, na.rm = TRUE), 1)
    }

    if ("rate" %in% names(df)) {
      stats$mean_rate[i] <- round(mean(df$rate, na.rm = TRUE), 4)
    } else if (all(c("mc", "cov") %in% names(df))) {
      valid <- df$cov > 0
      stats$mean_rate[i] <- round(mean(df$mc[valid] / df$cov[valid], na.rm = TRUE), 4)
    }
  }

  return(stats)
}


#' Print a formatted sample statistics table
#' @keywords internal
#' @noRd
.print_sample_table <- function(stats, indent = "  ") {
  # Truncate long sample names
  display_names <- ifelse(nchar(stats$sample) > 20,
                          paste0(substr(stats$sample, 1, 17), "..."),
                          stats$sample)

  cat(sprintf("%s%-20s %10s %8s %8s %9s\n",
              indent, "Sample", "Sites", "MeanCov", "MedCov", "MeanRate"))
  cat(indent, strrep("-", 58), "\n", sep = "")

  for (i in seq_len(nrow(stats))) {
    cat(sprintf("%s%-20s %10s %8.1f %8.1f %9.4f\n",
                indent,
                display_names[i],
                format(stats$n_sites[i], big.mark = ","),
                stats$mean_cov[i],
                stats$median_cov[i],
                stats$mean_rate[i]))
  }
}


#' Print site type breakdown
#' @keywords internal
#' @noRd
.print_site_types <- function(df) {
  if (is.null(df) || !"site" %in% names(df)) return(invisible(NULL))

  site_counts <- table(df$site)
  if (length(site_counts) == 0) return(invisible(NULL))

  cat(strrep("-", 40), "\n")
  cat("Site types (from first sample):\n")
  for (st in names(site_counts)) {
    cat(sprintf("  %-8s %s sites\n", st, format(site_counts[st], big.mark = ",")))
  }
}
