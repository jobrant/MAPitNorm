#' Export Methylation Data to BigWig Format
#'
#' Exports single-base resolution methylation data to BigWig format for
#' genome browser visualization. Can export individual samples or aggregate
#' replicates by group. Since BigWig format does not support strand, data
#' is exported as unstranded by default, or optionally split into separate
#' forward and reverse strand files.
#'
#' @param meth_data List object returned by \code{\link{load_data}} containing
#'   methylation data with a \code{sample_metadata} attribute.
#' @param type Character string specifying the score to write.
#'   One of \code{"meth"} (methylation rate) or \code{"cov"} (coverage).
#' @param out Directory path for output BigWig files. Created if it does not exist.
#' @param group_name Group ID to aggregate (required if
#'   \code{aggregate_replicates = TRUE}).
#' @param sample_name Sample ID to export (required if
#'   \code{aggregate_replicates = FALSE}).
#' @param genome Genome build. One of \code{"hg38"}, \code{"hg19"},
#'   \code{"mm10"}, \code{"mm39"}. Default is \code{"hg38"}.
#'   Ignored when \code{chrom_sizes_file} is provided.
#' @param chrom_sizes_file Path to a custom chromosome sizes file
#'   (tab-separated: chr, length). If provided, overrides \code{genome}.
#' @param aggregate_replicates Logical. If \code{TRUE}, aggregates replicates
#'   within a group using coverage-weighted mean rates. Default is \code{FALSE}.
#' @param min_coverage Minimum coverage threshold for including sites.
#'   Default is 1.
#' @param split_strand Logical. If \code{TRUE}, creates separate BigWig files
#'   for forward and reverse strands. If \code{FALSE} (default), all data is
#'   exported as unstranded.
#'
#' @return Character vector of created BigWig file path(s). If
#'   \code{split_strand = TRUE}, returns paths for both strand files.
#'
#' @details
#' BigWig format does not natively support strand information. When
#' \code{split_strand = FALSE} (default), all sites are exported in a single
#' unstranded file. When \code{split_strand = TRUE}, two files are created
#' with \code{_fwd} and \code{_rev} suffixes.
#'
#' For aggregated replicates, methylation rates are calculated as a
#' coverage-weighted mean: \code{sum(mc) / sum(cov)} across all replicates
#' in the group.
#'
#' This function requires the \pkg{rtracklayer} package. Install it from
#' Bioconductor with \code{BiocManager::install("rtracklayer")} if needed.
#'
#' @examples
#' \dontrun{
#' # Load data
#' all_samples <- load_data(
#'   dir_path = "path/to/files",
#'   sample_sheet = "samples.csv"
#' )
#'
#' # Export single sample methylation rates
#' create_bigwig(
#'   meth_data = all_samples,
#'   out = "bigwigs",
#'   sample_name = "allc_M1N2",
#'   type = "meth"
#' )
#'
#' # Export aggregated replicates, split by strand
#' create_bigwig(
#'   meth_data = all_samples,
#'   out = "bigwigs",
#'   aggregate_replicates = TRUE,
#'   group_name = "M1",
#'   type = "meth",
#'   split_strand = TRUE
#' )
#' }
#'
#' @seealso \code{\link{create_bedgraph}} for bedGraph export,
#'   \code{\link{load_data}} for loading methylation data.
#'
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqinfo<- seqlevels keepSeqlevels
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table rbindlist copy fread as.data.table :=
#' @export
create_bigwig <- function(
    meth_data,
    type = c("meth", "cov"),
    out,
    group_name = NULL,
    sample_name = NULL,
    genome = "hg38",
    chrom_sizes_file = NULL,
    aggregate_replicates = FALSE,
    min_coverage = 1,
    split_strand = FALSE)
{
  # Check rtracklayer availability
  
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "The 'rtracklayer' package is required for BigWig export.\n",
      "Install it with: BiocManager::install('rtracklayer')"
    )
  }
  
  type <- match.arg(type)
  
  # Validate inputs
  if (!is.list(meth_data)) {
    stop("'meth_data' must be a list as returned by load_data().")
  }
  
  if (missing(out)) {
    stop("'out' directory must be specified.")
  }
  
  # Create output directory if needed
  out <- normalizePath(out, mustWork = FALSE)
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE)
    message("Created output directory: ", out)
  }
  
  # Load chromosome sizes 
  chrom_sizes <- .load_chrom_sizes(genome, chrom_sizes_file)
  
  # Use NA for genome name when a custom file is provided
  genome_label <- if (!is.null(chrom_sizes_file)) NA_character_ else genome
  
  genome_info <- GenomeInfoDb::Seqinfo(
    seqnames   = chrom_sizes$seqname,
    seqlengths = chrom_sizes$seqlength,
    genome     = genome_label
  )
  
  # Retrieve sample metadata 
  sample_metadata <- attr(meth_data, "sample_metadata")
  if (is.null(sample_metadata)) {
    stop("No 'sample_metadata' attribute found in meth_data. ",
         "Ensure meth_data was created by load_data().")
  }
  # Ensure metadata is a data.table for consistent syntax
  if (!data.table::is.data.table(sample_metadata)) {
    sample_metadata <- data.table::as.data.table(sample_metadata)
  }
  
  # Extract the data to export 
  if (aggregate_replicates) {
    export_dt <- .aggregate_group(meth_data, sample_metadata, group_name)
    file_base <- group_name
  } else {
    export_dt <- .extract_single_sample(meth_data, sample_name)
    file_base <- sample_name
  }
  
  # Chromosome name harmonisation 
  export_dt <- .harmonize_chr_names(export_dt)
  
  # Coverage filter 
  n_before <- nrow(export_dt)
  export_dt <- export_dt[cov >= min_coverage]
  n_after <- nrow(export_dt)
  
  if (n_after == 0) {
    stop("No sites remaining after filtering at min_coverage = ", min_coverage)
  }
  if (n_before > n_after) {
    message(sprintf("Filtered %d sites below min_coverage = %d (%d sites remaining)",
                    n_before - n_after, min_coverage, n_after))
  }
  
  # Write BigWig file(s) 
  if (split_strand) {
    out_files <- .write_bigwig_by_strand(export_dt, type, out, file_base, genome_info)
  } else {
    out_file <- .write_bigwig(export_dt, type, out, file_base,
                              genome_info, strand_label = "")
    out_files <- out_file
  }
  
  message("BigWig export complete.")
  return(invisible(out_files))
}


# Internal helper functions -----------------------------------------------

#' Load chromosome sizes from file or built-in genome
#' @return data.table with seqname and seqlength columns
#' @keywords internal
#' @noRd
.load_chrom_sizes <- function(genome, chrom_sizes_file) {
  if (!is.null(chrom_sizes_file)) {
    if (!file.exists(chrom_sizes_file)) {
      stop("Chromosome sizes file not found: ", chrom_sizes_file)
    }
    return(data.table::fread(chrom_sizes_file,
                             col.names = c("seqname", "seqlength")))
  }
  
  valid_genomes <- c("hg38", "hg19", "mm10", "mm39")
  if (!genome %in% valid_genomes) {
    stop("Unsupported genome: '", genome, "'. ",
         "Use one of: ", paste(valid_genomes, collapse = ", "),
         ", or provide chrom_sizes_file.")
  }
  
  builtin_file <- system.file("extdata",
                              paste0(genome, ".chrom.sizes"),
                              package = "MAPitNorm")
  if (builtin_file == "" || !file.exists(builtin_file)) {
    stop("Built-in chromosome sizes file not found for genome: ", genome)
  }
  
  data.table::fread(builtin_file, col.names = c("seqname", "seqlength"))
}


#' Aggregate replicates within a group
#' @return data.table of aggregated methylation data
#' @keywords internal
#' @noRd
.aggregate_group <- function(meth_data, sample_metadata, group_name) {
  if (is.null(group_name)) {
    stop("'group_name' is required when aggregate_replicates = TRUE.")
  }
  
  samples_in_group <- sample_metadata[group_id == group_name, sample_id]
  
  if (length(samples_in_group) == 0) {
    stop("No samples found for group '", group_name, "'. ",
         "Available groups: ",
         paste(unique(sample_metadata$group_id), collapse = ", "))
  }
  
  # Check all requested samples exist in the data
  missing <- setdiff(samples_in_group, names(meth_data))
  if (length(missing) > 0) {
    stop("Samples listed in metadata but missing from data: ",
         paste(missing, collapse = ", "))
  }
  
  message(sprintf("Aggregating %d replicates for group '%s': %s",
                  length(samples_in_group), group_name,
                  paste(samples_in_group, collapse = ", ")))
  
  combined <- data.table::rbindlist(meth_data[samples_in_group])
  
  combined[, .(
    rate = sum(mc) / sum(cov),
    cov  = sum(cov),
    mc   = sum(mc)
  ), by = .(chr, pos, strand, site)]
}


#' Extract a single sample with validation
#' @return data.table copy of the sample
#' @keywords internal
#' @noRd
.extract_single_sample <- function(meth_data, sample_name) {
  if (is.null(sample_name)) {
    stop("'sample_name' is required when aggregate_replicates = FALSE.")
  }
  
  if (!sample_name %in% names(meth_data)) {
    stop("Sample '", sample_name, "' not found in data. ",
         "Available samples: ",
         paste(names(meth_data), collapse = ", "))
  }
  
  message("Exporting sample: ", sample_name)
  data.table::copy(meth_data[[sample_name]])
}


#' Harmonize chromosome names to UCSC style
#' @return data.table with chr column updated in place
#' @keywords internal
#' @noRd
.harmonize_chr_names <- function(dt) {
  sample_chr <- unique(dt$chr)
  has_chr_prefix <- any(grepl("^chr", sample_chr))
  
  if (!has_chr_prefix) {
    message("Converting Ensembl-style chromosome names to UCSC style (adding 'chr' prefix)...")
    dt[, chr := paste0("chr", chr)]
    dt[chr == "chrMT", chr := "chrM"]
  }
  
  return(dt)
}


#' Write a single BigWig file
#' @return File path of created BigWig
#' @keywords internal
#' @noRd
.write_bigwig <- function(dt, type, out, file_base, genome_info, strand_label = "") {
  score_col <- if (type == "meth") "rate" else "cov"
  
  # BigWig does not support strand — set to unstranded
  gr <- GenomicRanges::GRanges(
    seqnames = dt$chr,
    ranges   = IRanges::IRanges(start = dt$pos, width = 1),
    strand   = "*",
    score    = dt[[score_col]]
  )
  
  # Restrict to chromosomes present in genome info
  shared_levels <- intersect(
    GenomeInfoDb::seqlevels(gr),
    GenomeInfoDb::seqlevels(genome_info)
  )
  
  if (length(shared_levels) == 0) {
    stop("No overlapping chromosome names between data and genome info. ",
         "Check chromosome naming convention.")
  }
  
  dropped <- setdiff(GenomeInfoDb::seqlevels(gr), shared_levels)
  if (length(dropped) > 0) {
    message("Dropping chromosomes not in genome info: ",
            paste(dropped, collapse = ", "))
  }
  
  gr <- GenomeInfoDb::keepSeqlevels(gr, shared_levels, pruning.mode = "coarse")
  GenomeInfoDb::seqinfo(gr) <- genome_info[GenomeInfoDb::seqlevels(gr)]
  
  out_file <- file.path(out, paste0(file_base, strand_label, "_", type, ".bw"))
  
  rtracklayer::export.bw(gr, out_file)
  message("Wrote: ", out_file)
  
  return(out_file)
}


#' Write separate BigWig files for forward and reverse strands
#' @return Character vector of file paths
#' @keywords internal
#' @noRd
.write_bigwig_by_strand <- function(dt, type, out, file_base, genome_info) {
  out_files <- character(0)
  
  fwd_dt <- dt[strand == "+"]
  rev_dt <- dt[strand == "-"]
  
  if (nrow(fwd_dt) > 0) {
    fwd_file <- .write_bigwig(fwd_dt, type, out, file_base,
                              genome_info, strand_label = "_fwd")
    out_files <- c(out_files, fwd_file)
  } else {
    message("No forward-strand sites to export.")
  }
  
  if (nrow(rev_dt) > 0) {
    rev_file <- .write_bigwig(rev_dt, type, out, file_base,
                              genome_info, strand_label = "_rev")
    out_files <- c(out_files, rev_file)
  } else {
    message("No reverse-strand sites to export.")
  }
  
  # Handle unstranded sites
  unstranded_dt <- dt[!strand %in% c("+", "-")]
  if (nrow(unstranded_dt) > 0) {
    message(sprintf("Note: %d sites with unrecognized strand values exported separately.",
                    nrow(unstranded_dt)))
    unstrand_file <- .write_bigwig(unstranded_dt, type, out, file_base,
                                   genome_info, strand_label = "_unstranded")
    out_files <- c(out_files, unstrand_file)
  }
  
  if (length(out_files) == 0) {
    stop("No sites to export on any strand.")
  }
  
  return(out_files)
}