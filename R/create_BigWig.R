#' Export Methylation Data to BigWig Format
#'
#' Exports single-base resolution methylation rate data to BigWig format for
#' genome browser visualization. Can export individual samples or aggregate
#' replicates by group.
#'
#' @param data List object returned by \code{\link{load_data}} containing
#'   methylation data with sample_metadata attribute
#' @param out Directory path for output BigWig files
#' @param genome Genome build. One of "hg38", "hg19", "mm10", "mm39". Default is "hg38"
#' @param sample_name Sample ID to export (for single sample export)
#' @param aggregate_replicates Logical. If TRUE, aggregates replicates within a group.
#'   Default is FALSE
#' @param group_name Group ID to aggregate (required if aggregate_replicates = TRUE)
#' @param min_coverage Minimum coverage threshold for including sites. Default is 1
#' @param type Character string for type of value to include as score in bigwig.
#'        Should be one of "cov" or "meth" (for coverage score or methylation rate score, respectively.)

#'
#' @return Character vector of created file path(s)
#'
#' @details
#' For aggregated replicates, methylation rates are calculated as a coverage-weighted
#' mean: sum(methylated_counts) / sum(total_coverage) across all replicates in the group.
#'
#' @examples
#' \dontrun{
#' # Load data
#' all_samples <- load_data(
#'   dir_path = "path/to/files",
#'   sample_sheet = "samples.csv"
#' )
#'
#' # Export single sample
#' create_bigwig(
#'   data = all_samples,
#'   out = "bigwigs",
#'   sample_name = "allc_M1N2",
#'   type = "meth"
#' )
#'
#' # Export aggregated replicates
#' create_bigwig(
#'   data = all_samples,
#'   out = "bigwigs",
#'   aggregate_replicates = TRUE,
#'   group_name = "M1",
#'   type = "rate"
#' )
#' }
#'
#' @importFrom GenomeInfoDb seqinfo seqinfo<- seqlevels keepSeqlevels
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer export.bw import.bw
#' @importFrom data.table rbindlist copy
#' @export
create_bigwig <- function(
    data,
    type=c("meth","cov"),
    out,
    group_name = NULL,
    sample_name = NULL,
    genome = "hg38",
    aggregate_replicates = FALSE,
    min_coverage = 1)

  {
  # Match and check type argument
  type <- match.arg(type)
  # Only check for the genome-specific package (in Suggests)
  genome_pkg <- switch(genome,
                       "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
                       "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
                       "mm10" = "BSgenome.Mmusculus.UCSC.mm10",
                       "mm39" = "BSgenome.Mmusculus.UCSC.mm39",
                       stop("Unsupported genome. Use 'hg38', 'hg19', 'mm10', or 'mm39'")
  )

  if (!requireNamespace(genome_pkg, quietly = TRUE)) {
    stop(paste0("Package ", genome_pkg, " not installed. Install with: BiocManager::install('", genome_pkg, "')"))
  }

  # Access the BSgenome object from the package namespace
  genome_obj <- getExportedValue(genome_pkg, genome_pkg)
  genome_info <- seqinfo(genome_obj)

    # Get metadata from data object
    metadata <- attr(data, "sample_metadata")
    if (is.null(metadata)) {
      stop("No sample_metadata attribute found in data object")
    }

    # Get the data to export
    if (aggregate_replicates) {
      if (is.null(group_name)) stop("group_name required for aggregation")

      samples_to_agg <- metadata[group_id == group_name, sample_id]
      dt_list <- data[samples_to_agg]
      combined <- rbindlist(dt_list)

      # Aggregate: coverage-weighted mean
      sample_data <- combined[, .(
        rate = sum(mc) / sum(cov),
        cov = sum(cov)
      ), by = .(chr, pos, strand, site)]

      file_base <- paste0(group_name, "_aggregated")
    }
    else {
      if (is.null(sample_name)) stop("sample_name required for single sample export")
      sample_data <- data[[sample_name]]
      file_base <- sample_name
    }

    # Filter by coverage
    sample_data <- sample_data[cov >= min_coverage]

    # Export function
    export <- function(dt, strand_label = "") {
      # Map type to actual column name
      score_col <- if(type == "meth") "rate" else "cov"
      gr <- GRanges(
        seqnames = dt$chr,
        ranges = IRanges(start = dt$pos, width = 1),
        strand = dt$strand,
        score = dt[[score_col]]  # Use the selected type column
      )

      # Only keep chromosomes that exist in genome_info
      gr <- keepSeqlevels(gr, intersect(seqlevels(gr), seqlevels(genome_info)),
                          pruning.mode = "coarse")
      seqinfo(gr) <- genome_info[seqlevels(gr)]

      out_file <- file.path(out,
                            paste0(file_base, strand_label, "_", type, ".bw"))
      export.bw(gr, out_file)
      return(out_file)
    }

    # Export
      return(export(sample_data))

  }

