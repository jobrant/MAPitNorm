# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

createBedgraphCpp <- function(sample_df, type, out_path, group_name = NULL, sample_name = NULL) {
    .Call(`_MAPitNorm_createBedgraphCpp`, sample_df, type, out_path, group_name, sample_name)
}

readMethylationFile <- function(filename) {
    .Call(`_MAPitNorm_readMethylationFile`, filename)
}

readMethylationFiles <- function(filenames) {
    .Call(`_MAPitNorm_readMethylationFiles`, filenames)
}

testCpp <- function() {
    .Call(`_MAPitNorm_testCpp`)
}

