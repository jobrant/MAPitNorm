# Load all samples using C++ implementation (default)
# On HiPerGator
#all_samples <- load_data(
#  dir_path = "/blue/cancercenter-dept/shared/MAPitNorm/allc_files/",
#  sample_sheet = "/blue/cancercenter-dept/shared/MAPitNorm/M-Series_batches.csv"
#)

all_samples <- load_data(
  dir_path = "tests/testthat/test-data/",
  sample_sheet = "tests/testthat/test-data/M-Series_batches.csv"
)

# Load data
all_samples <- load_data(
  dir_path = "tests/testthat/test-data/",
  sample_sheet = "tests/testthat/test-data/M-Series_batches.csv"
)

# Prepare test data with chr prefix
test_samples <- lapply(all_samples, function(dt) {
  dt_copy <- copy(dt[chr == "1" & pos < 1000000])
  dt_copy[, chr := paste0("chr", chr)]
  return(dt_copy)
})
attr(test_samples, "sample_metadata") <- attr(all_samples, "sample_metadata")

# Run all 4 test combinations and store results
test_results <- list(
  single_meth = create_bigwig(
    data = test_samples,
    out = "tests/testthat/test-results/",
    genome = "hg38",
    sample_name = "allc_M1N2",
    type = "meth"
  ),
  single_cov = create_bigwig(
    data = test_samples,
    out = "tests/testthat/test-results/",
    genome = "hg38",
    sample_name = "allc_M1N2",
    type = "cov"
  ),
  agg_meth = create_bigwig(
    data = test_samples,
    out = "tests/testthat/test-results/",
    genome = "hg38",
    aggregate_replicates = TRUE,
    group_name = "M1",
    type = "meth"
  ),
  agg_cov = create_bigwig(
    data = test_samples,
    out = "tests/testthat/test-results/",
    genome = "hg38",
    aggregate_replicates = TRUE,
    group_name = "M1",
    type = "cov"
  )
)

# Define validation tests
validate_bigwig <- function(bw_file, test_name) {
  cat("\n=== Testing:", test_name, "===\n")

  # Check file exists
  if (!file.exists(bw_file)) {
    cat("ERROR: File does not exist\n")
    return(FALSE)
  }
  cat("File exists:", bw_file, "\n")

  # Load BigWig
  bw <- import.bw(bw_file)
  cat("Loaded", length(bw), "positions\n")

  # Test specific position: chr1:67231
  test_pos <- 67231
  bw_at_pos <- bw[seqnames(bw) == "chr1" & start(bw) == test_pos]

  if (length(bw_at_pos) > 0) {
    cat("Position chr1:", test_pos, "= ", bw_at_pos$score, "\n")

    # Get expected value from test data
    is_agg <- grepl("^agg_", test_name)
    type_type <- if(grepl("_meth$", test_name)) "meth" else "cov"

    if (is_agg) {
      # Aggregate validation
      metadata <- attr(test_samples, "sample_metadata")
      m1_samples <- metadata[group_id == "M1", sample_id]
      m1_at_pos <- lapply(m1_samples, function(s) {
        test_samples[[s]][chr == "chr1" & pos == test_pos]
      })
      m1_combined <- rbindlist(m1_at_pos)

      if (type_type == "meth") {
        expected <- sum(m1_combined$mc) / sum(m1_combined$cov)
      } else {
        expected <- sum(m1_combined$cov)
      }
    } else {
      # Single sample validation
      original <- test_samples$allc_M1N2[chr == "chr1" & pos == test_pos]
      # Map type to actual column name
      score_col <- if(type_type == "meth") "rate" else "cov"
      expected <- original[[score_col]]
    }

    match <- all.equal(expected, bw_at_pos$score)
    cat("Expected:", expected, "| Match:", match, "\n")
    return(isTRUE(match))
  } else {
    cat("Position not found in BigWig\n")
    return(FALSE)
  }
}

# Run validation on all results
cat("\n########## VALIDATION RESULTS ##########\n")
validation_results <- lapply(names(test_results), function(name) {
  validate_bigwig(test_results[[name]], name)
})
names(validation_results) <- names(test_results)

# Summary
cat("\n########## SUMMARY ##########\n")
print(validation_results)
all_passed <- all(unlist(validation_results))
cat("\nAll tests passed:", all_passed, "\n")

# Test 5: Export all groups at once (just test it runs)
cat("\n########## TESTING ALL GROUPS EXPORT ##########\n")
metadata <- attr(test_samples, "sample_metadata")
tryCatch({
  for (group in unique(metadata$group_id)) {
    create_bigwig(
      data = test_samples,
      out = "tests/testthat/test-results/",
      genome = "hg38",
      aggregate_replicates = TRUE,
      group_name = group
    )
  }
  cat("All groups exported successfully\n")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})
