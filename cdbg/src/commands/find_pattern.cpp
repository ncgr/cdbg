// std
#include <chrono>  // duration_cast, high_resolution_clock, milliseconds
#include <iostream>  // cout, endl
// local
#include "cdbg/cdbg.hpp"  // CDBG
#include "cdbg/io/implicit_stream.hpp"  // load_implicit

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::cout;
using std::endl;
using cdbg::io::load_implicit;


namespace cdbg {
namespace commands {


void find_pattern(const string& filename_graph, const string& filename_pattern)
{
  // Load Graph and input
  CDBG g = load_implicit(filename_graph);
  ifstream patternfile(filename_pattern);
  string p;
  vector<string> pattern;
  auto t1 = high_resolution_clock::now();
  auto t2 = high_resolution_clock::now();
  auto time_search_pattern = t2-t2;
  auto time_documents_of_first_node = t2-t2;
  auto time_documents_of_all_nodes = t2-t2;
  uint64_t number_patterns = 0;
  uint64_t number_found = 0;
  while (patternfile >> p) {
    ++number_patterns;
    vector<uint64_t> node_sequences;
    uint64_t tmp;
    t1 = high_resolution_clock::now();
    tie(node_sequences, tmp) = g.find_nodes(p);
    t2 = high_resolution_clock::now();
    time_search_pattern += t2-t1;
    if (node_sequences.size()) {
      cout << "Pattern '" << p << "' occurs in the following nodes: ";
      cout << node_sequences[0];
      for (uint64_t i = 1; i < node_sequences.size(); ++i) {
        cout << ", " << node_sequences[i];
      }
      cout << endl;
      ++number_found;
      t1 = high_resolution_clock::now();
      vector<uint64_t> seq = g.sequences_in_node(node_sequences.front());
      t2 = high_resolution_clock::now();
      time_documents_of_first_node += t2-t1;
      t1 = high_resolution_clock::now();
      for (const auto& nodeid : node_sequences) {
        vector<uint64_t> seq = g.sequences_in_node(nodeid);
        cout << "Node " << nodeid << " corresponds to a substring that occurs in the following sequences: ";
        cout << seq[0];
        for (uint64_t i = 1; i < seq.size(); ++i) {
          cout << ", " << seq[i];
        }
        cout << endl;
      }
      t2 = high_resolution_clock::now();
      time_documents_of_all_nodes += t2-t1;
    } else {
      cout << "Pattern '" << p << "' does not occur." << endl;
    }
    cout << endl;
  }
  cout << "Found " << number_found << " of " << number_patterns << " Pattern" << endl;
  cout << setw(10) << duration_cast<milliseconds>(time_search_pattern).count() << "ms to search pattern." << endl;
  cout << setw(10) << duration_cast<milliseconds>(time_documents_of_first_node).count() << "ms to list documents of the first node." << endl;
  cout << setw(10) << duration_cast<milliseconds>(time_documents_of_all_nodes).count() << "ms to list documents of all nodes." << endl;
}


}  // commands
}  // cdbg
