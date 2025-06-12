#include <Rcpp.h>
#include <fstream>
#include <string>
using namespace Rcpp;

// [[Rcpp::export]]
String createBedgraphCpp(DataFrame sample_df, 
                         std::string type,
                         std::string out_path,
                         Nullable<String> group_name = R_NilValue,
                         Nullable<String> sample_name = R_NilValue) {
  
  // Extract vectors from data frame
  CharacterVector chr_vec = sample_df["chr"];
  NumericVector pos_vec = sample_df["pos"];
  
  // Determine values based on the type
  NumericVector value_vec;
  if (type == "meth") {
    value_vec = sample_df["rate"];
  } else {
    value_vec = sample_df["cov"];
  }
  
  // Build file name
  std::string file_prefix;
  if (sample_name.isNotNull()) {
    if (group_name.isNotNull()) {
      file_prefix = as<std::string>(group_name) + "_" + as<std::string>(sample_name);
    } else {
      file_prefix = as<std::string>(sample_name);
    }
  } else {
    file_prefix = "sample_" + std::to_string(std::time(nullptr));
  }
  
  std::string output_file = out_path + "/" + 
                           file_prefix + "_" + 
                           type + ".bedGraph";
  
  // Open output file
  std::ofstream outfile(output_file);
  if (!outfile.is_open()) {
    stop("Could not open file for writing: " + output_file);
  }
  
  // Write bedGraph data line by line
  int n = chr_vec.size();
  for (int i = 0; i < n; i++) {
    outfile << "chr" << as<std::string>(chr_vec[i]) << "\t"
            << pos_vec[i] - 1 << "\t"  // BED format is 0-based
            << pos_vec[i] << "\t"
            << value_vec[i] << "\n";
  }
  
  outfile.close();
  return output_file;
}