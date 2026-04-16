#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <climits>
#include <algorithm>
#include <map>
#include <cassert>
using namespace std;

int main(int argc, char *argv[]) {
  string input, output;
  map<int, vector<int>> partition_to_reorder_map;

  input = argv[1];
  output = argv[2];

  std::ifstream infile(input);
  if (!infile.is_open()) {
    std::cout << "Couldn't open file " << input << std::endl;
    std::exit(-2);
  }

  int row_id = 0, part_id;
  while (infile >> part_id) {
    partition_to_reorder_map[part_id].push_back(row_id);
    row_id += 1;
  }
  infile.close();

  ofstream outfile(output);
  if(!outfile) {
    cout << "Cannot open the " << output << " output file!" << endl;
    exit(-3);
  }

  int last_part = -1;
  for(auto const &item : partition_to_reorder_map) {
    assert(item.first > last_part && "Partition ids should be accessed in sorted order.");
    for(auto const &r : item.second) {
      outfile << r << "\n";
    }
    last_part = item.first;
  }
  outfile.close();

  cout << "MAX Partition-id: " << last_part << endl;

  return 0;
}
