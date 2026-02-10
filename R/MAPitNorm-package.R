#' @keywords internal
"_PACKAGE"

## Shared imports used across the package
#' @importFrom data.table data.table as.data.table is.data.table fread
#'   fwrite setnames setDT setkey setkeyv haskey copy rbindlist := setDTthreads
#' @importFrom stats median quantile aggregate na.omit
#' @importFrom utils head setTxtProgressBar txtProgressBar combn 
#' @importFrom methods is new
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqinfo<- seqlevels keepSeqlevels
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowRanges colData
#' @importFrom S4Vectors SimpleList
#' @importFrom stringr str_split
#' @importFrom purrr map2 transpose
#' @importFrom parallel detectCores
#' @importFrom progress progress_bar
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.newpage grid.text gpar
#' @importFrom reshape2 dcast
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggplot scale_fill_viridis_c aes geom_point geom_abline
#' geom_density geom_hex geom_bar facet_wrap labs theme theme_minimal element_text
#' @importFrom Rcpp evalCpp
#' 
NULL
