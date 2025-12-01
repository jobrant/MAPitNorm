# Load all samples using C++ implementation (default)
all_samples <- load_data(
  dir_path = "/blue/cancercenter-dept/shared/MAPitNorm/allc_files/",
  sample_sheet = "/blue/cancercenter-dept/shared/MAPitNorm/M-Series_batches.csv"
)

# Test the Function
# Prepare test data with chr prefix and maintain structure
# Subsample chr1, first 1M bases from all samples
test_samples <- lapply(all_samples, function(dt) {
  dt_copy <- copy(dt[chr == "1" & pos < 1000000])
  dt_copy[, chr := paste0("chr", chr)]
  return(dt_copy)
})

# Attach metadata
attr(test_samples, "sample_metadata") <- attr(all_samples, "sample_metadata")

# Test 1: Single sample export
export_to_bigwig(
  data = test_samples,
  output_dir = ".",
  genome = "hg38",
  sample_name = "allc_M1N2"
)

# Test 2: Aggregate replicates (all M1 samples)
export_to_bigwig(
  data = test_samples,
  output_dir = ".",
  genome = "hg38",
  aggregate_replicates = TRUE,
  group_id = "M1"
)

# Test 3: Export all groups at once
metadata <- attr(test_samples, "sample_metadata")
for (group in unique(metadata$group_id)) {
  export_to_bigwig(
    data = test_samples,
    output_dir = ".",
    genome = "hg38",
    aggregate_replicates = TRUE,
    group_id = group
  )
}

# Checks
# Load it back into R
bw <- import.bw("allc_M1N2.bw")
head(bw)

# View first few positions from original data
test_data[12:20]

# Example: position chr1:67231 has rate=0.3333333
# Check the BigWig has it
bw[seqnames(bw) == "chr1" & start(bw) == 67231]  # positions from the data

# Also check in IGV genome browser

# Test aggregation math
# Load aggregated BigWig
bw_agg <- import.bw("M1_aggregated.bw")

# Pick a position to check - chr1:69426
test_pos <- 69426

# Get the value from the aggregated BigWig
bw_value <- bw_agg[seqnames(bw_agg) == "chr1" & start(bw_agg) == test_pos]
bw_value$score

# Manually calculate what it should be from the original M1 samples
metadata <- attr(test_samples, "sample_metadata")
m1_samples <- metadata[group_id == "M1", sample_id]

# Get this position from all M1 replicates
m1_at_pos <- lapply(m1_samples, function(s) {
  test_samples[[s]][chr == "chr1" & pos == test_pos]
})

# Combine and show
m1_combined <- rbindlist(m1_at_pos, idcol = "sample")
print(m1_combined)

# Calculate aggregated rate manually: sum(mc) / sum(cov)
manual_rate <- sum(m1_combined$mc) / sum(m1_combined$cov)
cat("\nManual calculation:", manual_rate, "\n")
cat("BigWig value:", bw_value$score, "\n")
cat("Match:", all.equal(manual_rate, bw_value$score), "\n")


