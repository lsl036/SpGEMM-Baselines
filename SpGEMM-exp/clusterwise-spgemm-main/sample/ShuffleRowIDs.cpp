#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <iterator>
#include <ctime>
#include <cmath>
#include <string>
#include <sstream>
#include <random>
#include <chrono>

#include <filesystem>
namespace fs = std::filesystem;

#include "../utility.h"
#include "../CSR.h"
#include "../multiply.h"

#include "../hash_mult.h"
#include "sample_common.hpp"

using namespace std;

#define VALUETYPE double
#define INDEXTYPE int64_t

int main(int argc, char *argv[]) {
  bool save_to_file = false;
  INDEXTYPE cluster_size = 8;
  CSR<INDEXTYPE, VALUETYPE> A_csr, B_csr;
  string output_path_str = "";

  // Parse command line arguments
  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-s" || arg == "--save") {
      save_to_file = true;
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      output_path_str = argv[++i];
    }
  }

  if (argc < 3) {
    cout << "Usage: ./shuffle {gen|binary|text} <matrix-file> [OPTIONS]" << endl;
    cout << "Options:" << endl;
    cout << "  -s, --save              Save shuffled row IDs to file" << endl;
    cout << "  -o, --output <path>     Output directory path (overrides SHUFFLED_DATA_PATH)" << endl;
    cout << "" << endl;
    cout << "Examples:" << endl;
    cout << "  ./shuffle text matrix.mtx -s" << endl;
    cout << "  ./shuffle text matrix.mtx -s -o /path/to/output" << endl;
    return -1;
  }

  // Determine output path: command line argument > environment variable
  if (output_path_str.empty()) {
    char* output_path_env = std::getenv("SHUFFLED_DATA_PATH");
    if (!output_path_env) {
      cerr << "Error: Output path not specified. Use -o/--output option or set SHUFFLED_DATA_PATH environment variable." << endl;
      return -1;
    }
    output_path_str = output_path_env;
  }

  // Create output directory if it doesn't exist
  if (save_to_file) {
    fs::path output_dir(output_path_str);
    if (!fs::exists(output_dir)) {
      fs::create_directories(output_dir);
      cout << "Created output directory: " << output_path_str << endl;
    }
  }

  /* Generating input matrices based on argument */
  SetInputMatricesAsCSR(A_csr, argv, cluster_size);
  A_csr.Sorted();

  vector<INDEXTYPE> row_id(A_csr.rows);
  for (int i = 0; i < A_csr.rows; i += 1) row_id[i] = i;

  double start, end, sec;
  start = omp_get_wtime();
  mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
  shuffle(row_id.begin(), row_id.end(), rng);
  end = omp_get_wtime();
  sec = (end - start);

  string filename = parseFileNameFromPath(argv[2]);
  cout << "Dataset " << filename << " takes: " << sec << " seconds." << endl;

  if(save_to_file) {
    string filename_suffix = create_filename_suffix(argv);
    string outfile = output_path_str + '/' + filename_suffix;

    std::fstream file(outfile, std::ios::out);
    for (INDEXTYPE i = 0; i < A_csr.rows; ++i) {
      file << row_id[i] << endl;
    }
    file.close();
    cout << "[DONE] Writing to file!" << endl;
  }

  A_csr.make_empty();

  return 0;
}
