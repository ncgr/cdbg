// std
#include <algorithm>  // min_element
#include <fstream>  // ifstream, ofstream
#include <string>
#include <tuple>  // tie
#include <vector>  // begin, end
// sdsl
#include <sdsl/config.hpp>  // cache_config
#include <sdsl/util.hpp>  // sdsl::util
// local
#include "cdbg/cdbg.hpp"  // CDBG
#include "../create_datastructures.hpp"  // create_bwt, create_da, create_sa,
                                         // create_text


using std::begin;
using std::end;
using std::ifstream;
using std::min_element;
using std::ofstream;
using std::string;
using std::tie;
using std::vector;
using sdsl::cache_config;


namespace cdbg {
namespace commands {


void construct(
  const string& inputfile,
  const string& outputfile,
  const string& kfilename,
  bool with_document_array)
{
  uint64_t min_length = 0;
  // Create datastructures
  cache_config config(true, ".", "tmp");
  {
    // Get input
    vector<uint64_t> sequences = create_text(config, inputfile, true);
    auto min = min_element(begin(sequences), end(sequences));
    if (min != end(sequences)) {
      min_length = *min;
    }
    // Get sa
    bool fast = false; // space-efficient!
    create_sa(config, fast);
    // Get bwt
    create_bwt(config);
    // Get document array
    if (with_document_array) {
      create_da(config, sequences);
    }
  }
  // Read k-values
  ifstream kfile(kfilename);
  uint64_t k;
  while (kfile >> k) {
    if (min_length < k) {
      cerr << "k=" << k << " must smaller than sequence length";
      cerr << ", but in input file '" << inputfile << "' there is a sequence with length ";
      cerr << min_length << " - this k-values will be skipped." << endl;
    } else {
      // Create graph
      CDBG g(config, k, with_document_array);
      // Store graph
      ofstream out(outputfile+".k"+to_string(k)+".bin");
      g.serialize(out);
    }
  }
  // Delete files
  if (config.delete_files) {
    sdsl::util::delete_all_files(config.file_map);
  }
}


}  // commands
}  // cdbg
