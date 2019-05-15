#ifndef CREATE_DATASTRUCTURES_HPP
#define CREATE_DATASTRUCTURES_HPP

// std
#include <string>
#include <vector>
// sdsl
#include <sdsl/config.hpp>  // cache_config


using std::string;
using std::vector;
using sdsl::cache_config;


namespace cdbg {


vector<uint64_t> create_text(cache_config&, const string&, bool=true);
bool k_smaller_than_sequences(
  const vector<uint64_t>&,
  const string&,
  const string&);
void create_sa(cache_config&, const bool);
void create_bwt(cache_config&);
void create_da(cache_config&, const vector<uint64_t>&);
uint64_t create_datastructures(
  cache_config&,
  const string&,
  const string&,
  const bool);


}  // cdbg


#endif
