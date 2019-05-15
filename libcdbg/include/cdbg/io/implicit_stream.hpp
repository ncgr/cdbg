#ifndef IMPLICIT_STREAM_HPP
#define IMPLICIT_STREAM_HPP


// std
#include <string>
// local
#include "cdbg/cdbg.hpp"  // CDBG


using std::string;
using cdbg::CDBG;


namespace cdbg {
namespace io {


CDBG load_implicit(const string& filename);


}  // io
}  // cdbg


#endif
