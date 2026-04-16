#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <ctime>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <random>
#include <cstring>
#include <typeinfo>
//#include "mmio.h"
using namespace std;

///g++ mtx_to_el.cc -o mtx_to_el -std=c++11

typedef int64_t IT;
typedef double NT;

int main(int argc, char *argv[]) {
  string input, output;

  if (argc < 3) {
    cout
        << "Normal usage: ./mtx_to_el <input-mtx> <output-el>"
        << endl;
    return -1;
  }

  input = argv[1];          // input file
  output = argv[2];  // output file

//  std::ifstream infile(input);
//  if (!infile.is_open()) {
//    std::cout << "Couldn't open file " << input << std::endl;
//    std::exit(-2);
//  }

  ofstream out_file(output);
  if(!out_file) {
    cout << "Cannot open the unique-single-output file!" << endl;
    exit(-3);
  }

  bool isSymmetric = false;
  ifstream infile(input.c_str());
  char line[1024];
  char c = infile.get();
  while(c == '%')
  {
    infile.getline(line,1024);
    if (strstr(line, "symmetric")) {
      isSymmetric = true;
      cout << "Matrix is symmetric" << endl;
    }
    c = infile.get();
  }
  infile.unget();

  infile.getline(line,1024);
//  cout << line << endl;
  IT m,n,nnz;
  if (typeid(IT) == typeid(int)) {
    sscanf(line, "%d %d %d", &m, &n, &nnz);
  }
  else if (typeid(IT) == typeid(long long int)) {
    sscanf(line, "%lld %lld %lld", &m, &n, &nnz);
  }
  else {
    sscanf(line, "%lld %lld %lld", &m, &n, &nnz);
  }

  if (isSymmetric) {
    nnz *= 2;
  }

//  cout << "come so far..." << endl;
//  cout << m << " " << n << " " << nnz << endl;
  int64_t row, col;
  double val;
//    cout << "mem allocation dome..." << endl;
  if (infile.is_open())
  {
    IT cnz = 0;	// current number of nonzeros
    while (! infile.eof() && cnz < nnz)
    {
      infile.getline(line,1024);
      char *ch = line;
      row = (IT)(atoi(ch));
      ch = strchr(ch, ' ');
      ch++;
      col = (IT)(atoi(ch));
      ch = strchr(ch, ' ');
      if (ch != NULL) {
        ch++;
        /* Read third word (value data)*/
        val = (NT)(atof(ch) + 1.0);
        ch = strchr(ch, ' ');
      }
      else {
        val = 1.0;
      }
      row--;
      col--;
//      cout << row << " " << col << endl;
//      exit(-1);
      out_file << row << " " << col << endl;
      if (isSymmetric) {
        if (col != row) {
          cnz++;
          out_file << col << " " << row << endl;
        }
        else {
          nnz--;
        }
      }
      ++cnz;
    }
    assert(cnz == nnz);
  }

  out_file.flush();
  out_file.close();
  infile.close();

  return 0;
}

