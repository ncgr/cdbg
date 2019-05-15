// std
#include <string>
// sdsl
#include <sdsl/io.hpp>  // load_from_file
// local
#include "cdbg/cdbg.hpp"  // CDBG


using std::string;
using cdbg::CDBG;
using sdsl::load_from_file;


namespace cdbg {
namespace io {


CDBG load_implicit(const string& filename) {
  CDBG g;
  load_from_file(g, filename);
  return g;
}


}  // io
}  // cdbg
