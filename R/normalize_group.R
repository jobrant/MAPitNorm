# File: R/normalize_group.R

#' Normalize Group Data
#'
#' Internal function to normalize coverage and methylation counts for replicate groups.
#'
#' @param replicates List of replicate data frames.
#'
#' @param sf Scaling factors to apply.
#'
#' @return List of normalized replicate data frames.
#'
#' @import purrr
#'
#' @keywords internal
#'
normalize_group <- function(replicates, sf) {
  map2(replicates, sf[names(replicates)], function(df, s) {
    if(is.na(s)) stop("Scaling factor is NA")
    dt <- data.table::copy(df)
    dt[, c("cov", "mc") := list(cov/s, mc/s)]
    dt[, rate := mc/cov]
    dt
  })
}
