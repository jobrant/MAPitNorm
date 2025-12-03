#' Export Methylation Data to BigWig Format
#'
#' Exports single-base resolution methylation rate data to BigWig format for
#' genome browser visualization. Can export individual samples or aggregate
#' replicates by group.
#'
#' @param data List object returned by \code{\link{load_data}} containing
#'   methylation data with sample_metadata attribute
#' @param output_dir Directory path for output BigWig files
#' @param genome Genome build. One of "hg38", "hg19", "mm10", "mm39". Default is "hg38"
#' @param sample_name Sample ID to export (for single sample export)
#' @param aggregate_replicates Logical. If TRUE, aggregates replicates within a group.
#'   Default is FALSE
#' @param group_id Group ID to aggregate (required if aggregate_replicates = TRUE)
#' @param min_coverage Minimum coverage threshold for including sites. Default is 1
#' @param split_strand Logical. If TRUE, creates separate BigWig files for plus
#'   and minus strands. Default is FALSE
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
#' export_to_bigwig(
#'   data = all_samples,
#'   output_dir = "bigwigs",
#'   sample_name = "allc_M1N2"
#' )
#'
#' # Export aggregated replicates
#' export_to_bigwig(
#'   data = all_samples,
#'   output_dir = "bigwigs",
#'   aggregate_replicates = TRUE,
#'   group_id = "M1"
#' )
#' }
#'
#' @importFrom rtracklayer export.bw
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqinfo seqlevels keepSeqlevels
#' @importFrom data.table rbindlist copy
#' @export
export_to_bigwig <- function(
    data,
    output_dir,
    genome = "hg38",
    sample_name = NULL,
    aggregate_replicates = FALSE,
    group_id = NULL,
    min_coverage = 1)
  {
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

  genome_obj <- get(genome_pkg)
  genome_info <- seqinfo(genome_obj)

    # Get metadata from data object
    metadata <- attr(data, "sample_metadata")
    if (is.null(metadata)) {
      stop("No sample_metadata attribute found in data object")
    }

    # Get the data to export
    if (aggregate_replicates) {
      if (is.null(group_id)) stop("group_id required for aggregation")

      # Use intermediate variable to avoid column name confusion
      target_group <- group_id
      samples_to_agg <- metadata[group_id == target_group, sample_id]
      dt_list <- data[samples_to_agg]
      combined <- rbindlist(dt_list)

      # Aggregate: coverage-weighted mean
      sample_data <- combined[, .(
        rate = sum(mc) / sum(cov),
        cov = sum(cov)
      ), by = .(chr, pos, strand, site)]

      file_base <- paste0(target_group, "_aggregated")
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
      gr <- GRanges(
        seqnames = dt$chr,
        ranges = IRanges(start = dt$pos, width = 1),
        strand = dt$strand,
        score = dt$rate
      )

      # Only keep chromosomes that exist in genome_info
      gr <- keepSeqlevels(gr, intersect(seqlevels(gr), seqlevels(genome_info)),
                          pruning.mode = "coarse")
      seqinfo(gr) <- genome_info[seqlevels(gr)]

      out_file <- file.path(output_dir,
                            paste0(file_base, strand_label, ".bw"))
      export.bw(gr, out_file)
      return(out_file)
    }

    # Export
      return(export(sample_data))

  }

