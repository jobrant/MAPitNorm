% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize_normalization.R
\name{visualize_normalization}
\alias{visualize_normalization}
\title{Visualize Normalization Effects}
\usage{
visualize_normalization(
  original_data,
  within_norm,
  full_norm,
  output_file = NULL
)
}
\arguments{
\item{original_data}{Original data list before normalization}

\item{within_norm}{Data after within-group normalization}

\item{full_norm}{Data after both normalizations}

\item{output_file}{Optional file path to save plots}
}
\value{
A grid arrangement of four diagnostic plots
}
\description{
Visualize Normalization Effects
}
