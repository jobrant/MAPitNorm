% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_methylation_rate.R
\name{normalize_methylation_rates}
\alias{normalize_methylation_rates}
\title{Normalize Methylation Rates}
\usage{
normalize_methylation_rates(
  data_list,
  group_names = NULL,
  within_groups = TRUE,
  between_groups = TRUE,
  within_alpha = 0.3,
  between_alpha = 0.5,
  sites_per_quantile = 1000,
  max_quantiles = 50,
  diagnostics = TRUE
)
}
\arguments{
\item{data_list}{List of normalized coverage data}

\item{group_names}{Optional character vector of group names. If NULL (default),
will be extracted from the metadata attached to data_list.}

\item{within_groups}{Logical indicating whether to normalize within groups (default = TRUE)}

\item{between_groups}{Logical indicating whether to normalize between groups (default = TRUE)}

\item{within_alpha}{Numeric between 0 and 1 controlling structure preservation within groups (default = 0.3)}

\item{between_alpha}{Numeric between 0 and 1 controlling structure preservation between groups (default = 0.5)}

\item{sites_per_quantile}{Target number of sites per quantile (default = 1000).
For large datasets (>100k sites), this provides a good balance between
normalization granularity and statistical power.}

\item{max_quantiles}{Maximum number of quantiles to use (default = 50)}

\item{diagnostics}{Logical indicating whether to print diagnostic information}
}
\value{
List of lists containing normalized methylation data
}
\description{
Normalizes methylation rates using quantile-based normalization to account for
non-linear enzyme efficiency differences.
}
\details{
This function performs a two-step normalization:
1. Within-group normalization to handle technical variation between replicates
2. Between-group normalization to adjust for global efficiency differences
while preserving biological differences. The higher between_alpha preserves
more of the biological differences between groups.
}
