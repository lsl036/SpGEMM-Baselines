#include <iostream>
#include <set>

#include "sparsebase/format/csr.h"
#include "sparsebase/format/format_order_one.h"
#include "sparsebase/format/format_order_two.h"
#include "sparsebase/object/object.h"
#include "sparsebase/reorder/amd_reorder.h"
#include "sparsebase/bases/iobase.h"
#include "sparsebase/bases/reorder_base.h"

using namespace std;
using namespace sparsebase;
using namespace sparsebase::reorder;
using namespace sparsebase::bases;
using namespace io;
using namespace bases;
using namespace format;
using namespace reorder;

using vertex_type = int;
using edge_type = int;
using value_type = float;

int main(int argc, char *argv[]) {
  if (argc < 2) {
//    cout << "Usage: ./amd_order <uedgelist_file>\n";
    cout << "Usage: ./amd_order <mtx_file> <reorder_output_file>\n";
    cout
        << "Hint: You can use the edgelist: examples/data/com-dblp.uedgelist\n";
    return 1;
  }
  cout << "F t re  s sp r e!" << endl;
  string file_name = argv[1];
  string out_file_name = argv[2];
  context::CPUContext cpu_context;

  cout << "********************************" << endl;

//  cout << "Reading graph from " << file_name << "..." << endl;
//  // The name of the edge list file in disk
//  std::string matrix_filename(argv[1]);
//  // Read the edge list file into a CSR object
//  CSR<vertex_type, edge_type, value_type>* csr = IOBase::ReadMTXToCSR<vertex_type, edge_type, value_type>(matrix_filename);
//
//  // get_dimensions() returns a vector with the dimension of
//  // each order of the format object
//  vertex_type num_rows = csr->get_dimensions()[0];
//  vertex_type num_columns = csr->get_dimensions()[1];
//  edge_type num_non_zeros = csr->get_num_nnz();
//
//  std::cout << "Matrix has "
//            << num_rows << " rows, "
//            << num_columns << " columns, and "
//            << num_non_zeros << " non-zeros" << std::endl;

  cout << "Reading graph from " << file_name << "..." << endl;
//  object::Graph<vertex_type, edge_type, value_type> g;
//  g.ReadConnectivityFromEdgelistToCSR(file_name);
//
//  cout << "Number of vertices: " << g.n_ << endl;
//  cout << "Number of edges: " << g.m_ << endl;

  auto csr = sparsebase::bases::IOBase::ReadMTXToCSR<vertex_type, edge_type, value_type>(file_name, true);
  cout << "Number of vertices: " << csr->get_dimensions()[0] << endl;
  cout << "Number of edges: " << csr->get_num_nnz() << endl;

  cout << "********************************" << endl;

  cout << "Generating AMD ordering..." << endl;

  // Context representing the GPU with ID 0 in the system
  context::CPUContext cpu;

  // Create a parameters object to store special parameters specific
  AMDReorderParams p{};

  // AMD reordering
  // Create an inverse permutation array  of the matrix `csr`
//  auto *con =
//      g.get_connectivity()
//          ->AsAbsolute<format::CSR<vertex_type, edge_type, value_type>>();
//  vertex_type num_rows = con->get_dimensions()[0];
//  vertex_type num_columns = con->get_dimensions()[0];
//  edge_type num_non_zeros = con->get_num_nnz();
  vertex_type num_rows = csr->get_dimensions()[0];
  vertex_type num_columns = csr->get_dimensions()[1];
  edge_type num_non_zeros = csr->get_num_nnz();
  std::cout << "Matrix has "
            << num_rows << " rows, "
            << num_columns << " columns, and "
            << num_non_zeros << " non-zeros" << std::endl;
//  vertex_type * amd_reorder = bases::ReorderBase::Reorder<AMDReorder>(p, con, {&cpu}, true);
  vertex_type * amd_reorder = bases::ReorderBase::Reorder<AMDReorder>(p, csr, {&cpu}, true);

  cout << "********************************" << endl;

  cout << "Checking the correctness of the ordering..." << endl;
  bool order_is_correct = true;
  set<vertex_type> ids;
  for (vertex_type i = 0; i < num_rows && order_is_correct; i++) {
    vertex_type i_order = amd_reorder[i];
    if (i_order < num_rows && ids.find(i_order) == ids.end()) {
      ids.insert(i_order);
    } else {
      cout << "AMD ordering is incorrect!";
      order_is_correct = false;
      return 1;
    }
  }
  if (ids.size() > num_rows) {
    cout << "AMD ordering is incorrect!";
    order_is_correct = false;
  }
  if (order_is_correct) {
    cout << "Order is correct!" << endl;
  } else {
    cout << "AMD ordering is incorrect!";
    order_is_correct = false;
    return 1;
  }

  // Save the result in file
  std::ofstream amd_order_file(out_file_name.c_str());
  assert(amd_order_file.is_open() && "AMD order output file is not open!");
  std::copy(&amd_reorder[0], &amd_reorder[num_rows], std::ostream_iterator<vertex_type>(amd_order_file, "\n"));
  amd_order_file.close();

  return 0;
}
