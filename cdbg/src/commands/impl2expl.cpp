// std
#include <fstream>  // ofstream
#include <iostream>  // cerr, endl
#include <string>
#include <tuple>  // tie
#include <vector>
// local
#include "../handle_graph.hpp"  // print_graph
#include "cdbg/cdbg.hpp"  // CDBG
#include "cdbg/io/implicit_stream.hpp"  // load_implicit


using std::cerr;
using std::endl;
using std::ofstream;
using std::string;
using std::tie;
using std::vector;
using cdbg::CDBG;
using cdbg::io::load_implicit;


namespace cdbg {
namespace commands {


bool impl2expl(const string& filename_graph, const string& filename_output)
{
  CDBG g = load_implicit(filename_graph);
  vector<node> graph;
  vector<uint64_t> start_nodes;
  tie(graph, start_nodes) = g.get_explicit_representation();
  ofstream output(filename_output+".dot");
  ofstream output_start_nodes(filename_output+".start_nodes.txt");
  if (!output.is_open() || !output_start_nodes.is_open()) {
    cerr << "Could not open '" << filename_output << ".dot' or '";
    cerr << filename_output << ".start_nodes.txt' for writing." << endl;
    return false;
  } else {
    print_graph(graph, start_nodes, output, output_start_nodes);
  }
  return true;
}


}  // commands
}  // cdbg
