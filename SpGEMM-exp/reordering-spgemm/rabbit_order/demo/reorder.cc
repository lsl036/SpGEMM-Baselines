//
// A demo program of reordering using Rabbit Order.
//
// Author: ARAI Junya <arai.junya@lab.ntt.co.jp> <araijn@gmail.com>
//

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/count.hpp>
#include "../rabbit_order.hpp"
#include "edge_list.hpp"
//#include <iostream>

using rabbit_order::vint;
typedef std::vector<std::vector<std::pair<vint, float> > > adjacency_list;

vint count_unused_id(const vint n, const std::vector<edge_list::edge>& edges) {
  std::vector<char> appears(n);
  for (size_t i = 0; i < edges.size(); ++i) {
    appears[std::get<0>(edges[i])] = true;
    appears[std::get<1>(edges[i])] = true;
  }
  return static_cast<vint>(boost::count(appears, false));
}

template<typename RandomAccessRange>
adjacency_list make_adj_list(const vint n, const RandomAccessRange& es) {
  using std::get;

  // Symmetrize the edge list and remove self-loops simultaneously
  std::vector<edge_list::edge> ss(boost::size(es) * 2);
  #pragma omp parallel for
  for (size_t i = 0; i < boost::size(es); ++i) {
    auto& e = es[i];
    if (get<0>(e) != get<1>(e)) {
      ss[i * 2    ] = std::make_tuple(get<0>(e), get<1>(e), get<2>(e));
      ss[i * 2 + 1] = std::make_tuple(get<1>(e), get<0>(e), get<2>(e));
    } else {
      // Insert zero-weight edges instead of loops; they are ignored in making
      // an adjacency list
      ss[i * 2    ] = std::make_tuple(0, 0, 0.0f);
      ss[i * 2 + 1] = std::make_tuple(0, 0, 0.0f);
    }
  }

  // Sort the edges
  __gnu_parallel::sort(ss.begin(), ss.end());

  // Convert to an adjacency list
  adjacency_list adj(n);
  #pragma omp parallel
  {
    // Advance iterators to a boundary of a source vertex
    const auto adv = [](auto it, const auto first, const auto last) {
      while (first != it && it != last && get<0>(*(it - 1)) == get<0>(*it))
        ++it;
      return it;
    };

    // Compute an iterator range assigned to this thread
    const int    p      = omp_get_max_threads();
    const size_t t      = static_cast<size_t>(omp_get_thread_num());
    const size_t ifirst = ss.size() / p * (t)   + std::min(t,   ss.size() % p);
    const size_t ilast  = ss.size() / p * (t+1) + std::min(t+1, ss.size() % p);
    auto         it     = adv(ss.begin() + ifirst, ss.begin(), ss.end());
    const auto   last   = adv(ss.begin() + ilast,  ss.begin(), ss.end());

    // Reduce edges and store them in std::vector
    while (it != last) {
      const vint s = get<0>(*it);

      // Obtain an upper bound of degree and reserve memory
      const auto maxdeg = 
          std::find_if(it, last, [s](auto& x) {return get<0>(x) != s;}) - it;
      adj[s].reserve(maxdeg);

      while (it != last && get<0>(*it) == s) {
        const vint t = get<1>(*it);
        float      w = 0.0;
        while (it != last && get<0>(*it) == s && get<1>(*it) == t)
          w += get<2>(*it++);
        if (w > 0.0)
          adj[s].push_back({t, w});
      }

      // The actual degree can be smaller than the upper bound
      adj[s].shrink_to_fit();
    }
  }

  return adj;
}

adjacency_list read_graph(const std::string& graphpath) {
  const auto edges = edge_list::read(graphpath);

  // The number of vertices = max vertex ID + 1 (assuming IDs start from zero)
  const auto n =
      boost::accumulate(edges, static_cast<vint>(0), [](vint s, auto& e) {
          return std::max(s, std::max(std::get<0>(e), std::get<1>(e)) + 1);});

  if (const size_t c = count_unused_id(n, edges)) {
    std::cerr << "WARNING: " << c << "/" << n << " vertex IDs are unused"
              << " (zero-degree vertices or noncontiguous IDs?)\n";
  }

  return make_adj_list(n, edges);
}

template<typename InputIt>
typename std::iterator_traits<InputIt>::difference_type
count_uniq(const InputIt f, const InputIt l) {
  std::vector<typename std::iterator_traits<InputIt>::value_type> ys(f, l);
  return boost::size(boost::unique(boost::sort(ys)));
}

double compute_modularity(const adjacency_list& adj, const vint* const coms) {
  const vint  n    = static_cast<vint>(adj.size());
  const auto  ncom = count_uniq(coms, coms + n);
  double      m2   = 0.0;  // total weight of the (bidirectional) edges

  std::unordered_map<vint, double[2]> degs(ncom);  // ID -> {all, loop}
  degs.reserve(ncom);

  #pragma omp parallel reduction(+:m2)
  {
    std::unordered_map<vint, double[2]> mydegs(ncom);
    mydegs.reserve(ncom);

    #pragma omp for
    for (vint v = 0; v < n; ++v) {
      const vint  c = coms[v];
      auto* const d = &mydegs[c];
      for (const auto e : adj[v]) {
        m2      += e.second;
        (*d)[0] += e.second;
        if (coms[e.first] == c) (*d)[1] += e.second;
      }
    }

    #pragma omp critical
    {
      for (auto& kv : mydegs) {
        auto* const d = &degs[kv.first];
        (*d)[0] += kv.second[0];
        (*d)[1] += kv.second[1];
      }
    }
  }
  assert(static_cast<intmax_t>(degs.size()) == ncom);

  double q = 0.0;
  for (auto& kv : degs) {
    const double all  = kv.second[0];
    const double loop = kv.second[1];
    q += loop / m2 - (all / m2) * (all / m2);
  }

  return q;
}

void detect_community(adjacency_list adj, const std::string &commu_path, bool save_to_file=false) {
  auto _adj = adj;  // copy `adj` because it is used for computing modularity

  std::cerr << "Detecting communities...\n";
  const double tstart = rabbit_order::now_sec();
  //--------------------------------------------
  auto       g = rabbit_order::aggregate(std::move(_adj));
  const auto c = std::make_unique<vint[]>(g.n());
  #pragma omp parallel for
  for (vint v = 0; v < g.n(); ++v)
    c[v] = rabbit_order::trace_com(v, &g);
  //--------------------------------------------
  std::cout << "Runtime for community detection [sec]: "
            << rabbit_order::now_sec() - tstart << std::endl;

  // Save the result in file
  if(save_to_file) {
    std::ofstream commu_file(commu_path.c_str());
    assert(commu_file.is_open() && "Community output file is not open!");
    std::copy(&c[0], &c[g.n()], std::ostream_iterator<vint>(commu_file, "\n"));
    commu_file.close();
  }
  else { // Print the result to stdout
    std::copy(&c[0], &c[g.n()], std::ostream_iterator<vint>(std::cout, "\n"));
  }

  std::cerr << "Computing modularity of the result...\n";
  const double q = compute_modularity(adj, c.get());
  std::cerr << "Modularity: " << q << std::endl;
}

void reorder(adjacency_list adj, const std::string &reorder_path, const std::string &commu_path, bool save_to_file = false) {
  std::cerr << "Generating a permutation...\n";
  const double tstart = rabbit_order::now_sec();
  //--------------------------------------------
  const auto g = rabbit_order::aggregate(std::move(adj));
  const auto p = rabbit_order::compute_perm(g, reorder_path, commu_path, save_to_file);
  //--------------------------------------------
  std::cout << "Runtime for permutation generation [sec]: "
            << rabbit_order::now_sec() - tstart << std::endl;

  // Print the result
//  if(!save_to_file) {
//    std::copy(&p[0], &p[g.n()], std::ostream_iterator<vint>(std::cout, "\n"));
//  }
}

void PrintHelp(const char *program_name)
{
    std::cerr << "Rabbit Order Graph Reordering Tool\n";
    std::cerr << "==================================\n\n";
    std::cerr << "Usage: " << program_name << " <graph-file> <reorder-file> <community-file> [OPTIONS]\n\n";
    std::cerr << "Required Arguments:\n";
    std::cerr << "  <graph-file>         Path to the edge-list format file (.el)\n";
    std::cerr << "  <reorder-file>       Path to save the reordering result (.rabbitorder)\n";
    std::cerr << "  <community-file>     Path to save the community detection result (.off)\n\n";
    std::cerr << "Optional Arguments:\n";
    std::cerr << "  -s, --save           Save the result to files (default: print to stdout)\n";
    std::cerr << "  -c, --community      Print community IDs instead of a new ordering\n";
    std::cerr << "  -h, --help           Show this help message\n\n";
    std::cerr << "Examples:\n";
    std::cerr << "  " << program_name << " graph.el reorder.rabbitorder community.off -s\n";
    std::cerr << "  " << program_name << " graph.el reorder.rabbitorder community.off -c -s\n";
    std::cerr << "  " << program_name << " graph.el reorder.rabbitorder community.off\n\n";
}

// ./reorder GRAPH_FILE REORDER_FILE COMMUNITY_FILE -s -c
int main(int argc, char* argv[]) {
  using boost::adaptors::transformed;

  // Check for help option first
  if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
    PrintHelp(argv[0]);
    return EXIT_SUCCESS;
  }

  // Parse command-line arguments
  if (argc < 4) {
    std::cerr << "Error: Missing required arguments.\n\n";
    PrintHelp(argv[0]);
    exit(EXIT_FAILURE);
  }

  const std::string graphpath = argv[1];
  const std::string reorder_path = argv[2];
  const std::string commu_path = argv[3];
  
  bool save_to_file = false;
  bool commode = false;
  
  // Parse optional arguments
  for (int i = 4; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "-s" || arg == "--save") {
      save_to_file = true;
    } else if (arg == "-c" || arg == "--community") {
      commode = true;
    } else if (arg == "-h" || arg == "--help") {
      PrintHelp(argv[0]);
      return EXIT_SUCCESS;
    } else {
      std::cerr << "Error: Unknown option: " << arg << "\n\n";
      PrintHelp(argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  std::cerr << "Number of threads: " << omp_get_max_threads() << std::endl;

  std::cerr << "Reading an edge-list file: " << graphpath << std::endl;
  auto       adj = read_graph(graphpath);
  const auto m   =
      boost::accumulate(adj | transformed([](auto& es) {return es.size();}),
                        static_cast<size_t>(0));
  std::cerr << "Number of vertices: " << adj.size() << std::endl;
  std::cerr << "Number of edges: "    << m          << std::endl;

  if (commode)
    detect_community(std::move(adj), commu_path, save_to_file);
  else
    reorder(std::move(adj), reorder_path, commu_path, save_to_file);

  return EXIT_SUCCESS;
}

