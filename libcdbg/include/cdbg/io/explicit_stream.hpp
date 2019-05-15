#ifndef EXPLICIT_STREAM_HPP
#define EXPLICIT_STREAM_HPP


// std
#include <fstream>  // ifstream, ofstream
#include <iostream>  // cerr, endl
#include <string>
#include <vector>


using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;


namespace cdbg {
namespace io {


template<class t_node>
void store_graph(const vector<t_node>& graph, const string& filename)
{
  ofstream out(filename);
  uint64_t graph_size = graph.size();
  out.write((char*)&graph_size, sizeof(graph_size));
  for (const auto& node : graph) {
    uint64_t len = node.len;
    uint64_t size_adj_list = node.adj_list.size();
    uint64_t size_pos_list = node.pos_list.size();
    out.write((char*)&len, sizeof(len));
    out.write((char*)&size_adj_list, sizeof(size_adj_list));
    out.write((char*)node.adj_list.data(), size_adj_list*sizeof(node.adj_list[0]));
    out.write((char*)&size_pos_list, sizeof(size_pos_list));
    out.write((char*)node.pos_list.data(), size_pos_list*sizeof(node.pos_list[0]));
  }
  if (!out) {
    cerr << "Something went wrong - storage didn't work as expected" << endl;
  }
}


template<class t_node>
void load_graph(vector<t_node>& graph, const string& filename)
{
  ifstream inp(filename);
  uint64_t graph_size=0;
  inp.read((char*)&graph_size, sizeof(graph_size));
  graph.resize(graph_size);
  for (uint64_t i = 0; i < graph.size(); ++i) {
    uint64_t len=0;
    uint64_t size_adj_list=0;
    uint64_t size_pos_list=0;
    inp.read((char*)&len, sizeof(len));
    graph[i].len = len;
    inp.read((char*)&size_adj_list, sizeof(size_adj_list));
    graph[i].adj_list.resize(size_adj_list);
    inp.read((char*)graph[i].adj_list.data(), size_adj_list*sizeof(graph[i].adj_list[0]));
    inp.read((char*)&size_pos_list, sizeof(size_pos_list));
    graph[i].pos_list.resize(size_pos_list);
    inp.read((char*)graph[i].pos_list.data(), size_pos_list*sizeof(graph[i].pos_list[0]));
  }
  if (!inp) {
    cerr << "Something went wrong - reading didn't work as expected" << endl;
  }
}


}  // io
}  // cdbg


#endif
