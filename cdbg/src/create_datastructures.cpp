// std
#include <algorithm>  // lower_bound
#include <iostream>  // cerr, endl
#include <fstream>  // ifstream
#include <string>
#include <vector>  // begin, end, vector
// sdsl
#include <sdsl/config.hpp>  // sdsl::conf, SE_SAIS, cache_config
#include <sdsl/construct.hpp>  // contains_no_zero_symbol, load_vector_from_file
#include <sdsl/construct_bwt.hpp>  // construct_bwt
#include <sdsl/construct_config.hpp>  // construct_config
#include <sdsl/construct_sa.hpp>  // construct_sa
#include <sdsl/int_vector.hpp>  // int_vector
#include <sdsl/int_vector_buffer.hpp>  // int_vector_buffer
#include <sdsl/io.hpp>  // cache_file_exists, cache_file_name, load_from_cache,
                        // register_cache_file, store_to_cache


using std::begin;
using std::cerr;
using std::end;
using std::endl;
using std::ifstream;
using std::lower_bound;
using std::string;
using std::vector;
using sdsl::SE_SAIS;
using sdsl::cache_config;
using sdsl::cache_file_exists;
using sdsl::cache_file_name;
using sdsl::construct_bwt;
using sdsl::construct_config;
using sdsl::construct_sa;
using sdsl::contains_no_zero_symbol;
using sdsl::int_vector;
using sdsl::int_vector_buffer;
using sdsl::load_from_cache;
using sdsl::load_vector_from_file;
using sdsl::register_cache_file;
using sdsl::store_to_cache;


namespace cdbg {


vector<uint64_t>
create_text(
  cache_config& config,
  const string& inputfile,
  bool calc_sequences=true)
{
  vector<uint64_t> sequences;
  // (1) Check, if the text is cached
  if (!cache_file_exists(sdsl::conf::KEY_TEXT, config)) {
    int_vector<8> text;
    uint64_t num_bytes = 1;
    load_vector_from_file(text, inputfile, num_bytes);
    bool skip = false;
    uint64_t target = 0;
    uint64_t sequence_length = 0;
    for (uint64_t i = 0; i < text.size(); ++i) {
      if (text[i] == '>') {
        skip = true;
        if (i > 0) {
          text[target] = 1;  // delimiter character $
          ++target;
          sequences.emplace_back(sequence_length);
          sequence_length = 0;
        }
      } else if (text[i] == '\n') {
        skip = false;
      } else if (!skip) {
        text[target] = text[i];  // shift text to cover previous description
        ++target;
        ++sequence_length;
      }
    }
    sequences.emplace_back(sequence_length);
    text.resize(target);
    if (contains_no_zero_symbol(text, inputfile)) {
      append_zero_symbol(text);
      store_to_cache(text, sdsl::conf::KEY_TEXT, config);
    }
  } else {
    config.delete_files = false;
    if (calc_sequences) {
      int_vector<8> text;
      load_from_cache(text, sdsl::conf::KEY_TEXT, config);
      for (uint64_t i = 0, len = 0; i < text.size(); ++i) {
        if (text[i] <= 1) {
          sequences.emplace_back(len);
          len = 0;
        } else {
          ++len;
        }
      }
    }
  }
  register_cache_file(sdsl::conf::KEY_TEXT, config);
  return sequences;
}


bool k_smaller_than_sequences(
  const vector<uint64_t>& sequences,
  const string& inputfile,
  const string& kfilename)
{
  if (!sequences.size()) {
    cerr << "Could not find a sequence." << endl;
    return false;
  }
  ifstream kfile(kfilename);
  uint64_t k;
  while (kfile >> k) {
    for (const auto& length : sequences) {
      if (length < k) {
        cerr << "k=" << k << " must smaller than sequence length, but in input file '" << inputfile << "' there is a sequence with length " << length << "." << endl;
        return false;
      }
    }
  }
  return true;
}


void create_sa(cache_config& config, const bool fast)
{
  // (2) Check, if the suffix array is cached
  if (!cache_file_exists(sdsl::conf::KEY_SA, config)) {
    if (!fast) {
      construct_config::byte_algo_sa = SE_SAIS;
    }
    construct_sa<8>(config);
  } else {
    config.delete_files = false;
  }
  register_cache_file(sdsl::conf::KEY_SA, config);
}


void create_bwt(cache_config& config)
{
  // (3) Check, if bwt is cached
  if (!cache_file_exists(sdsl::conf::KEY_BWT, config)) {
    construct_bwt<8>(config);
  } else {
    config.delete_files = false;
  }
  register_cache_file(sdsl::conf::KEY_BWT, config);
}


void create_da(cache_config& config, const vector<uint64_t>& sequences)
{
  // Load SA
  if (!cache_file_exists(sdsl::conf::KEY_SA, config)) {
    create_sa(config, true);
  }
  int_vector_buffer<> sa(cache_file_name(sdsl::conf::KEY_SA, config));
  // Convert length to endpositions
  vector<uint64_t> endpos = sequences;
  for (uint64_t i = 1; i < endpos.size(); ++i) {
    endpos[i] += (1+endpos[i-1]);
  }
  // Create Document Array
  uint8_t bit_width = 1;
  for (uint64_t j = 2; j < sequences.size(); j *= 2) {
    ++bit_width;
  }
  int_vector<> da(sa.size(), 0, bit_width);
  // Fill Document Array
  for (uint64_t i = 0; i < sa.size(); ++i) {
    uint64_t sa_value = sa[i];
    da[i] = lower_bound(begin(endpos), end(endpos), sa_value) - begin(endpos);
  }
  // Store Document Array
  store_to_cache(da, "DA", config);
  return;
}


uint64_t create_datastructures(
  cache_config& config,
  const string& inputfile,
  const string& kfilename,
  const bool fast)
{
  // Get input
  vector<uint64_t> sequences = create_text(config, inputfile, true);
  // Check input
  if (!k_smaller_than_sequences(sequences, inputfile, kfilename)) {
    return 1;
  }
  // Get sa
  create_sa(config, fast);
  // Get bwt
  create_bwt(config);
  return 0;
}


}  // cdbg
