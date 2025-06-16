#include <Rcpp.h>
#include <fstream>
#include <string>
#include <iomanip>
using namespace Rcpp;

// [[Rcpp::export]]
String createBedgraphCpp(DataFrame sample_df, 
                         std::string type,
                         std::string out_path,
                         Nullable<String> group_name = R_NilValue,
                         Nullable<String> sample_name = R_NilValue) {

    // Validate type
    if (type != "meth" && type != "cov") {
      stop("Invalid type: must be 'meth' or 'cov'");
    }

  // Determine the value column based on type
  std::string value_col = (type == "meth") ? "rate" : "cov";

  // Validate data frame columns
  CharacterVector col_names = sample_df.names();
  std::vector<std::string> required_cols = {"chr", "pos", value_col};
  for (const auto& col : required_cols) {
    if (!std::any_of(col_names.begin(), col_names.end(),
                     [&col](const String& name) { return name == col; })) {
      stop("Missing required column: " + col);
    }
  }
  
  // Extract vectors from data frame
  CharacterVector chr_vec = sample_df["chr"];
  IntegerVector pos_vec = sample_df["pos"];
  NumericVector value_vec = sample_df[value_col];

  // Validate vector lengths
  int n = chr_vec.size();
  if (pos_vec.size() != n || value_vec.size() != n) {
    stop("Column lengths do not match");
  }
  
  // Build file name
  std::string file_prefix;
  if (sample_name.isNotNull()) {
    file_prefix = as<std::string>(sample_name);
    if (group_name.isNotNull()) {
      file_prefix = as<std::string>(group_name) + "_" + file_prefix;
    }
  } else {
    file_prefix = "sample_" + std::to_string(Rcpp::runif(1, 0, 1000000)[0]);
  }
  
  std::string output_file = out_path + "/" + file_prefix + "_" + type + ".bedGraph";
  
  // Open output file
  std::ofstream outfile(output_file);
  if (!outfile.is_open()) {
    stop("Could not open file for writing: " + output_file);
  }
  
  // Write bedGraph data
  outfile << std::fixed << std::setprecision(6); // Match R's precision
  for (int i = 0; i < n; ++i) {
    outfile << "chr" << as<std::string>(chr_vec[i]) << "\t"
            << pos_vec[i] - 1 << "\t" // BED format is 0-based
            << pos_vec[i] << "\t"
            << value_vec[i] << "\n";
  }
  
  outfile.close();
  
  // Verify file was written
  if (!outfile.good()) {
    stop("Error writing to file: " + output_file);
  }
  
  return output_file;
}