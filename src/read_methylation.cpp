#include <Rcpp.h>
#include <fstream>
#include <string>
#include <zlib.h>
using namespace Rcpp;

// Helper function to read gzipped files
std::vector<std::string> readGzippedLines(const char* filename) {
  std::vector<std::string> lines;
  gzFile file = gzopen(filename, "rb");
  
  if (!file) {
    Rcpp::stop("Could not open file: %s", filename);
  }
  
  char buffer[1024];
  while (gzgets(file, buffer, sizeof(buffer)) != Z_NULL) {
    lines.push_back(std::string(buffer));
  }
  
  gzclose(file);
  return lines;
}

extern "C" {
  // [[Rcpp::export]]
  DataFrame readMethylationFile(std::string filename) {
    Rcout << "Reading file: " << filename << "\n";
    // Read all lines from the gzipped file
    std::vector<std::string> lines = readGzippedLines(filename.c_str());
    
    // Preallocate vectors for data
    int n_lines = lines.size();
    CharacterVector chr(n_lines);
    IntegerVector pos(n_lines);
    CharacterVector strand(n_lines);
    CharacterVector site(n_lines);
    IntegerVector mc(n_lines);
    IntegerVector cov(n_lines);
    
    // Parse each line
    for (int i = 0; i < n_lines; i++) {
      std::string line = lines[i];
      std::istringstream iss(line);
      std::string chr_val, strand_val, site_val;
      int pos_val, mc_val, cov_val;
      
      // Read tab-separated values
      if (std::getline(iss, chr_val, '\t') && 
          (iss >> pos_val) && iss.ignore() &&
          std::getline(iss, strand_val, '\t') &&
          std::getline(iss, site_val, '\t') &&
          (iss >> mc_val) && iss.ignore() &&
          (iss >> cov_val)) {
          
        // Store in vectors
        chr[i] = chr_val;
        pos[i] = pos_val;
        strand[i] = strand_val;
        site[i] = site_val;
        mc[i] = mc_val;
        cov[i] = cov_val;
      }
    }
    
    // Create and return data frame
    DataFrame df = DataFrame::create(
      _["chr"] = chr,
      _["pos"] = pos,
      _["strand"] = strand,
      _["site"] = site,
      _["mc"] = mc,
      _["cov"] = cov
    );
    
    return df;
  }


  // [[Rcpp::export]]
  List readMethylationFiles(CharacterVector filenames) {
    Rcout << "Entering readMethylationFiles with " << filenames.size() << " files\n";
    int n_files = filenames.size();
    List results(n_files);
    
    for (int i = 0; i < n_files; i++) {
      std::string filename = as<std::string>(filenames[i]);
      Rcout << "Processing file: " << filename << "\n";
      Rcpp::checkUserInterrupt(); // Allow user to cancel
      results[i] = readMethylationFile(filename);
    }
    
    Rcout << "Exiting readMethylationFiles\n";
    return results;
  } 

  // [[Rcpp::export]]
  void testCpp() {
    Rcout << "Test C++ function called successfully\n";
  }

}

