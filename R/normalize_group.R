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
  map2(replicates, sf, function(df, s) {

    # Print diagnostics
    print("Scaling factor being applied:")
    print(s)
    print("Sample of original values:")
    print(head(df$cov))
    print(head(df$mc))

    # Check for NA in scaling factor
    if(is.na(s)) stop("Scaling factor is NA")

    # Create a copy of the data.table
    df_new <- data.table::copy(df)

    # Modify the copy
    df_new[, `:=`(
      cov = cov / s,
      mc = mc / s,
      rate = mc / cov
    )]

    # Print sample of results
    print("Sample of normalized values:")
    print(head(df_new$cov))
    print(head(df_new$mc))
    print(head(df_new$rate))

    df_new
  })
}
