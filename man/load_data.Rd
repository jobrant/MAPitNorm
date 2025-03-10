% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data.R
\name{load_data}
\alias{load_data}
\title{Load and Preprocess Data}
\usage{
load_data(dir_path, type = c("GCH", "HCG"), cores = 1)
}
\arguments{
\item{dir_path}{Path to the directory containing the data files.}

\item{type}{Character string specifying which type of files to load.}

\item{cores}{Number of cores to use for parallel processing. Default is 1
(no parallel processing).}
}
\value{
A list of data.tables containing processed methylation data.
}
\description{
Load methylation data from files and preprocess them.
Files are expected to be named groupName_replicateName_site.tsv
Site is either GCH or HCG.
}
\examples{
\dontrun{
# Load GCH files sequentially
gch_samples <- load_data("path/to/files", type = "GCH")

# Load HCG files using parallel processing
hcg_samples <- load_data("path/to/files", type = "HCG", cores = 4)
}

}
