#include <iostream>
#include <set>

#include "sparsebase/bases/iobase.h"
#include "sparsebase/bases/reorder_base.h"
#include "sparsebase/format/csr.h"
#include "sparsebase/format/format_order_one.h"
#include "sparsebase/format/format_order_two.h"
#include "sparsebase/object/object.h"
#include "sparsebase/reorder/slashburn_reorder.h"

using namespace std;
using namespace sparsebase;
using namespace io;
using namespace bases;
using namespace format;
using namespace reorder;

using vertex_type = int64_t;
using edge_type = int64_t;
using value_type = double;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: ./slashburn_order <matrix_market_file> "
            "<reorder_output_file> <save_to_file>\n";
    cout << "Hint: You can use the edgelist: examples/data/com-dblp.mtx\n";
    return 1;
  }
  cout << "F t re  s sp r e!" << endl;
  string file_name = argv[1];
  string out_file_name = argv[2];
  bool save_to_file = false;
  if(argc == 4) save_to_file = (atoi(argv[3]) == 0) ? false : true;

  cout << "********************************" << endl;

  cout << "Reading graph from " << file_name << "..." << endl;
  auto csr =
      sparsebase::bases::IOBase::ReadMTXToCSR<vertex_type, edge_type,
          value_type>(file_name, true);
  cout << "Number of vertices: " << csr->get_dimensions()[0] << endl;
  cout << "Number of edges: " << csr->get_num_nnz() << endl;

  vertex_type num_rows = csr->get_dimensions()[0];
  edge_type num_non_zeros = csr->get_num_nnz();

  cout << "********************************" << endl;
  cout << "Generating SlashBurn ordering..." << endl;

  // A context representing the host system
  context::CPUContext cpu_context;

  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  // Create a parameters object to store special parameters specific
  int hubset_k_size =
      (int)(0.02 * (double)num_rows);  // 0.02 is used in RabbitOrder paper;
                                       // 0.001 is used in BEAR paper;
  if (hubset_k_size < 1) hubset_k_size = 1;
  SlashburnReorderParams sb_params(hubset_k_size, false, true);

  vertex_type *sb_reorder = bases::ReorderBase::Reorder<SlashburnReorder>(
      sb_params, csr, {&cpu_context}, true);

  chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
  chrono::duration<float> elapsed_time = chrono::duration_cast<chrono::duration<float>>(end - start);
  std::cout << "Reordering takes: " << elapsed_time.count() << " seconds" << std::endl;
  cout << "********************************" << endl;

  cout << "Checking the correctness of the ordering..." << endl;
  bool order_is_correct = true;
  set<vertex_type> ids;
  for (vertex_type i = 0; i < num_rows && order_is_correct; i++) {
    vertex_type i_order = sb_reorder[i];
    if (i_order < num_rows && ids.find(i_order) == ids.end()) {
      ids.insert(i_order);
    } else {
      cout << "SlashBurn ordering is incorrect!";
      order_is_correct = false;
      return 1;
    }
  }
  if (ids.size() > num_rows) {
    cout << "SlashBurn ordering is incorrect!";
    order_is_correct = false;
  }
  if (order_is_correct) {
    cout << "Order is correct!" << endl;
  } else {
    cout << "SlashBurn ordering is incorrect!";
    order_is_correct = false;
    return 1;
  }

  // Save the result in file
  if (save_to_file) {
    std::ofstream sb_order_file(out_file_name.c_str());
    assert(sb_order_file.is_open() && "SlashBurn order output file is not open!");
    std::copy(&sb_reorder[0], &sb_reorder[num_rows],
              std::ostream_iterator<vertex_type>(sb_order_file, "\n"));
    sb_order_file.close();
  }

  delete sb_reorder;

  return 0;
}