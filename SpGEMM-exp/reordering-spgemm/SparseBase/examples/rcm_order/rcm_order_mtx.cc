#include <iostream>
#include <set>

#include "sparsebase/bases/iobase.h"
#include "sparsebase/bases/reorder_base.h"
#include "sparsebase/format/csr.h"
#include "sparsebase/format/format_order_one.h"
#include "sparsebase/format/format_order_two.h"
#include "sparsebase/object/object.h"
#include "sparsebase/reorder/rcm_reorder.h"

using namespace std;
using namespace sparsebase;
using namespace io;
using namespace bases;
using namespace format;
using namespace reorder;

//using vertex_type = unsigned int;
//using edge_type = unsigned int;
//using value_type = unsigned int;
using vertex_type = int64_t;
using edge_type = int64_t;
using value_type = double;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: ./rcm_order <uedgelist_file>\n";
    cout
        << "Hint: You can use the edgelist: examples/data/com-dblp.uedgelist\n";
    return 1;
  }
  cout << "F t re  s sp r e!" << endl;
  string file_name = argv[1];
//  context::CPUContext cpu_context;

  cout << "********************************" << endl;

  cout << "Reading graph from " << file_name << "..." << endl;

//  object::Graph<vertex_type, edge_type, value_type> g;
//  g.ReadConnectivityFromEdgelistToCSR(file_name);
//  cout << "Number of vertices: " << g.n_ << endl;
//  cout << "Number of edges: " << g.m_ << endl;

  auto csr =
      sparsebase::bases::IOBase::ReadMTXToCSR<vertex_type, edge_type,
                                              value_type>(file_name, true);
  cout << "Number of vertices: " << csr->get_dimensions()[0] << endl;
  cout << "Number of edges: " << csr->get_num_nnz() << endl;

  vertex_type num_rows = csr->get_dimensions()[0];
  edge_type num_non_zeros = csr->get_num_nnz();

  cout << "********************************" << endl;

  cout << "Generating RCM ordering..." << endl;

  // A context representing the host system
  context::CPUContext cpu_context;

  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

  reorder::RCMReorder<vertex_type, edge_type, value_type> orderer;
//  auto *con =
//      g.get_connectivity()
//          ->AsAbsolute<format::CSR<vertex_type, edge_type, value_type>>();
//  vertex_type *order = orderer.GetReorder(con, {&cpu_context}, false);
//  auto xadj = con->get_row_ptr();
//  auto adj = con->get_col();
//  vertex_type n = con->get_dimensions()[0];
  vertex_type *order = orderer.GetReorder(csr, {&cpu_context}, false);

  chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
  chrono::duration<float> elapsed_time = chrono::duration_cast<chrono::duration<float>>(end - start);
  std::cout << "Reordering takes: " << elapsed_time.count() << " seconds" << std::endl;

//  auto xadj = csr->get_row_ptr();
//  auto adj = csr->get_col();
  vertex_type n = csr->get_dimensions()[0];

  cout << "********************************" << endl;

  cout << "Checking the correctness of the ordering..." << endl;
  bool order_is_correct = true;
  set<vertex_type> ids;
  for (vertex_type i = 0; i < n && order_is_correct; i++) {
    vertex_type i_order = order[i];
    if (i_order < n && ids.find(i_order) == ids.end()) {
      ids.insert(i_order);
    } else {
      cout << "RCM ordering is incorrect!";
      order_is_correct = false;
      return 1;
    }
  }
  if (ids.size() > n) {
    cout << "RCM ordering is incorrect!";
    order_is_correct = false;
  }
  if (order_is_correct) {
    cout << "Order is correct!" << endl;
  } else {
    cout << "RCM ordering is incorrect!";
    order_is_correct = false;
    return 1;
  }
  return 0;
}
