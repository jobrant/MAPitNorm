#' @keywords internal
"_PACKAGE"

## Shared imports used across the package
#' @importFrom data.table data.table as.data.table is.data.table fread
#'   fwrite setnames setDT setkey setkeyv haskey copy rbindlist :=
#' @importFrom stats median quantile
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @importFrom methods is
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqinfo<- seqlevels keepSeqlevels
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowRanges colData
#' @importFrom S4Vectors SimpleList
#' @importFrom stringr str_split
#' @importFrom purrr map2
#' @importFrom parallel detectCores
#' @importFrom progress progress_bar
NULL
