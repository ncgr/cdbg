#ifndef HANDLE_GRAPH
#define HANDLE_GRAPH

// std
#include <iostream>  // cerr, cout, endl, ostream
#include <vector>


using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::vector;


namespace cdbg {


template<class t_node>
void print_graph(
  const vector<t_node>& graph,
  const vector<uint64_t>& start_nodes,
  ostream& out=cout,
  ostream& out2=cout)
{
  // Print start nodes
  for (const auto& node : start_nodes) {
    out2 << node << endl;
  }
  uint64_t labels = 0;
  uint64_t edges = 0;
  out << "digraph G {\n";
  uint64_t node_number = 0;
  for (const auto& node : graph) {
    labels += graph[node_number].pos_list.size();
    edges += graph[node_number].adj_list.size();
    // Print node_number
    out << "  " << node_number << " ";
    // Print label
    out << "[label=\"";
    // if(node.pos_list.size()) // ToDo: Should not be necessary
    {
      out << node.pos_list[node.pos_list.size()-1];
      for (uint64_t j = node.pos_list.size()-2; j < node.pos_list.size()-1; --j) {
        out << "," << node.pos_list[j];
      }
    }
    out << ":" << node.len << "\"]\n";
    // Print edges
    for (uint64_t j = node.adj_list.size()-1; j < node.adj_list.size(); --j) {
      out << "  " << node_number << " -> " << node.adj_list[j] << "\n";
    }
    ++node_number;
  }
  out << "}" << endl;
}


template<class t_node>
bool verify_graph(vector<t_node> graph, uint64_t k, vector<uint64_t>& sequenzes) // call-by-value
{
  uint64_t start_pos = 1;
  bool done = false;
  while (!done) {
    // Search node start with position start_pos
    uint64_t start_node=graph.size();
    for (uint64_t i = 0; i < graph.size(); ++i) {
      if (graph[i].pos_list.size() > 0 && graph[i].pos_list.back() == start_pos) {
        start_node = i;
        break;
      }
    }
    if (start_node == graph.size()) {  // Next sequenze could not be found
      // Check if graph is now empty
      for (uint64_t i = 0; i < graph.size(); ++i) {
        if (graph[i].adj_list.size() || graph[i].pos_list.size()) {
          cerr << "Graph is not empty after covering " << sequenzes.size() << " sequenzes" << endl;
          cerr << "node " << i << " has " << graph[i].adj_list.size() << " outgoing edges and " << graph[i].pos_list.size() << " entry points" << endl;
          return false;
        }
      }
      done = true;
    } else {
      // Check current sequenze
      uint64_t cur_node = start_node;
      uint64_t cur_pos = start_pos;
      while (graph[cur_node].adj_list.size()) {
        uint64_t next_node = graph[cur_node].adj_list.back();
        graph[cur_node].adj_list.pop_back();
        graph[cur_node].pos_list.pop_back();
        uint64_t next_pos = cur_pos + graph[cur_node].len - k + 1;
        if (next_pos != graph[next_node].pos_list.back()) {
          cerr << "I was at node " << cur_node << " with len=" << graph[cur_node].len << " at position " << cur_pos << " and walked to node " << next_node << " - but could not find the position " << next_pos << " but found " << graph[next_node].pos_list.back() << " instead." << endl;
          return false;
        }
        cur_node = next_node;
        cur_pos = next_pos;
      }
      // Enter last node and setup for next sequenze
      graph[cur_node].pos_list.pop_back();
      cur_pos += graph[cur_node].len;
      sequenzes.emplace_back(cur_pos - start_pos);
      start_pos = cur_pos;
    }
  }
  return true;
}


template<class t_node>
vector<uint8_t> restore_text(
  vector<t_node> graph,
  const vector<uint8_t>& text,
  uint64_t k) // call-by-value
{
  vector<uint8_t> res;
  vector<uint64_t> sequenzes;
  res.reserve(text.size());
  uint64_t start_pos = 1;
  bool done = false;
  while (!done) {
    // Search node start with position start_pos
    uint64_t start_node=graph.size();
    for (uint64_t i = 0; i < graph.size(); ++i) {
      if (graph[i].pos_list.size() > 0 && graph[i].pos_list.back() == start_pos) {
        start_node = i;
        break;
      }
    }
    if (start_node == graph.size()) {  // Next sequenze could not be found
      // Check if graph is now empty
      for (uint64_t i = 0; i < graph.size(); ++i) {
        if (graph[i].adj_list.size() || graph[i].pos_list.size()) {
          cerr << "Graph is not empty after covering " << sequenzes.size() << " sequenzes" << endl;
          cerr << "node " << i << " has " << graph[i].adj_list.size() << " outgoing edges and " << graph[i].pos_list.size() << " entry points" << endl;
          res.resize(0);
          return res;
        }
      }
      res.emplace_back(0);
      done = true;
    } else {
      if (res.size()) {
        res.emplace_back(1);
      }
      // Check current sequenze
      uint64_t cur_node = start_node;
      uint64_t cur_pos = start_pos;
      while (graph[cur_node].adj_list.size()) {
        uint64_t next_node = graph[cur_node].adj_list.back();
        graph[cur_node].adj_list.pop_back();
        graph[cur_node].pos_list.pop_back();
        uint64_t next_pos = cur_pos + graph[cur_node].len - k + 1;
        if (next_pos != graph[next_node].pos_list.back()) {
          cerr << "I was at node " << cur_node << " with len=" << graph[cur_node].len << " at position " << cur_pos << " and walked to node " << next_node << " - but could not find the position " << next_pos << " but found " << graph[next_node].pos_list.back() << " instead." << endl;
          res.resize(1);
          return res;
        }
        // Restore Text
        uint64_t text_pos = graph[cur_node].pos_list.front();
        for (uint64_t i = 0; i < graph[cur_node].len-k+1; ++i) {
          res.emplace_back(text[text_pos-1+i]);
        }
        cur_node = next_node;
        cur_pos = next_pos;
      }
      // Restore Text
      uint64_t text_pos = graph[cur_node].pos_list.front();
      for (uint64_t i = 0; i < graph[cur_node].len-1; ++i) {
        res.emplace_back(text[text_pos-1+i]);
      }
      // Enter last node and setup for next sequenze
      graph[cur_node].pos_list.pop_back();
      cur_pos += graph[cur_node].len;
      sequenzes.emplace_back(cur_pos - start_pos);
      start_pos = cur_pos;
    }
  }
  return res;
}


}  // cdbg


#endif
