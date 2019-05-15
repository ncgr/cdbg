// std
#include <iostream>  // cerr
#include <string>
// local
#include "cdbg/cdbg.hpp"  // CDBG
#include "cdbg/io/implicit_stream.hpp"  //load_implicit


using std::cerr;
using std::string;
using cdbg::io::load_implicit;


namespace cdbg {
namespace commands {


void print_graph_details(const string& graphfile)
{
  cerr << endl << graphfile << ":" << endl;
  CDBG g = load_implicit(graphfile);
  g.print_graph_statistics(cerr);
}


}  // commands
}  // cdbg
