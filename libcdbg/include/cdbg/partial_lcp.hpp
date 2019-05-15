#ifndef PARTIAL_LCP_HPP
#define PARTIAL_LCP_HPP

// std
#include <queue>
#include <utility>  // swap
#include <vector>
// sdsl
#include <sdsl/int_vector.hpp>  // bit_vector, int_vector
#include <sdsl/rank_support_v.hpp>  // rank_support_v
#include <sdsl/select_support_mcl.hpp>  // select_support_mcl
#include <sdsl/select_support_scan.hpp>  // select_support_scan
#include <sdsl/structure_tree.hpp>  // structure_tree_node
#include <sdsl/util.hpp>  // sdsl::util
#include <sdsl/wt_helper.hpp>  // int_tree
#include <sdsl/wt_huff.hpp>  // wt_huff


using std::queue;
using std::swap;
using std::vector;
using sdsl::bit_vector;
using sdsl::int_tree;
using sdsl::int_vector;
using sdsl::rank_support_v;
using sdsl::select_support_mcl;
using sdsl::select_support_scan;
using sdsl::structure_tree_node;
using sdsl::wt_huff;


namespace cdbg {


enum lcp_value_enum {gt_k=0, lt_k=1, eq_k=2};


template<class t_wt>
int_vector<2> construct_partial_lcp(
  const t_wt& wt_bwt,
  vector<uint64_t>& C,
  uint64_t k)
{
  typedef int_vector<>::size_type size_type;
  lcp_value_enum marker = gt_k;
  int_vector<2> lcp(wt_bwt.size(), marker);
  marker = lt_k;
  
  uint64_t n = wt_bwt.size();  // Input length
  size_type lcp_value = 0;  // Current LCP value
  
  // Declare needed variables
  size_type intervals = 0;  // Number of intervals which are currently stored
  size_type intervals_new = 0;  // Number of new intervals
  
  queue<size_type> q;  // Queue for storing the intervals
  vector<bit_vector> dict(2);  // Bit_vector for storing the intervals
  size_type source = 0, target = 1;  // Defines which bit_vector is source and which is target
  bool queue_used = true;  // Defines whether a queue (true) or the bit_vectors (false) was used to store intervals
  size_type use_queue_and_wt = n/2048;  // If intervals < use_queue_and_wt, then use queue and wavelet tree
                                        // else use dictionary and wavelet tree
  
  size_type quantity;  // Quantity of characters in interval
  vector<unsigned char> cs(wt_bwt.sigma);  // List of characters in the interval
  vector<size_type> rank_c_i(wt_bwt.sigma);  // Number of occurrence of character in [0 .. i-1]
  vector<size_type> rank_c_j(wt_bwt.sigma);  // Number of occurrence of character in [0 .. j-1]
  
  // Save position of first LCP-value
  lcp[0] = marker;
  
  // Calculate first intervals
  interval_symbols(wt_bwt, 0, n, quantity, cs, rank_c_i, rank_c_j);
  for (size_type i = 0; i < quantity; ++i) {
    unsigned char c = cs[i];
    if (c == 1) {
      continue;
    }
    size_type a_new = C[c] + rank_c_i[i];
    size_type b_new = C[c] + rank_c_j[i];
    // Save LCP value and corresponding interval if not seen before
    if (!lcp[b_new]) {
      lcp[b_new] = marker;
      // Save interval
      q.push(a_new);
      q.push(b_new);
      ++intervals;
    }
  }
  for (size_type i = C[1]; i < C[2]; ++i) {
    size_type a_new = i;
    size_type b_new = i+1;
    // Save LCP value and corresponding interval if not seen before
    if (!lcp[b_new]) {
      lcp[b_new] = marker;
      // Save interval
      q.push(a_new);
      q.push(b_new);
      ++intervals;
    }
  }
  ++lcp_value;
  // Calculate LCP positions
  while (intervals && lcp_value <= k) {
    if (lcp_value == k) {
      marker = eq_k;
    }
    if (intervals < use_queue_and_wt && !queue_used) {
      sdsl::util::clear(dict[target]);
      // Copy from bitvector to queue
      size_type a2 = sdsl::util::next_bit(dict[source], 0);
      size_type b2 = sdsl::util::next_bit(dict[source], a2+1);
      while (b2 < dict[source].size()) {
        q.push((a2-1)>>1);
        q.push(b2>>1);
        // Get next interval
        a2 = sdsl::util::next_bit(dict[source], b2+1);
        b2 = sdsl::util::next_bit(dict[source], a2+1);
      }
      sdsl::util::clear(dict[source]);
    }
    if (intervals >= use_queue_and_wt && queue_used) {
      dict[source].resize(2*(n+1));
      sdsl::util::set_to_value(dict[source], 0);
      // Copy from queue to bitvector
      while (!q.empty()) {
        dict[source][(q.front()<<1)+1 ] = 1; q.pop();
        dict[source][(q.front()<<1)   ] = 1; q.pop();
      }
      dict[target].resize(2*(n+1));
      sdsl::util::set_to_value(dict[target], 0);
    }
  
    if (intervals < use_queue_and_wt) {
      queue_used = true;
      intervals_new = 0;
      while (intervals) {
        // Get next interval
        size_type a = q.front(); q.pop();
        size_type b = q.front(); q.pop();
        --intervals;
        interval_symbols(wt_bwt, a, b, quantity, cs, rank_c_i, rank_c_j);
  	    for (size_type i = 0; i < quantity; ++i) {
  		  unsigned char c = cs[i];
  		  size_type a_new = C[c] + rank_c_i[i];
  		  size_type b_new = C[c] + rank_c_j[i];
  		  // Save LCP value and corresponding interval if not seen before
  		  if (!lcp[b_new]) {
  		    lcp[b_new] = marker;
  		    // Save interval
  		    q.push(a_new);
  		    q.push(b_new);
  		    ++intervals_new;
  		  }
  	    }
      }
      intervals = intervals_new;
    } else {
      queue_used = false;
      intervals = 0;
      // Get next interval
      size_type a2 = sdsl::util::next_bit(dict[source], 0);
      size_type b2 = sdsl::util::next_bit(dict[source], a2+1);
      while (b2 < dict[source].size()) {
        interval_symbols(wt_bwt, ((a2-1)>>1), (b2>>1), quantity, cs, rank_c_i,
                           rank_c_j);
  	    for (size_type i = 0; i < quantity; ++i) {
  	      unsigned char c = cs[i];
  	      size_type a_new = C[c] + rank_c_i[i];
  	      size_type b_new = C[c] + rank_c_j[i];
  	      // Save LCP value if not seen before
  	      if (!lcp[b_new]) {
  	        lcp[b_new] = marker;
  	        // Save interval
  	        dict[target][(a_new<<1)+1] = 1;
  	        dict[target][(b_new<<1)  ] = 1;
  	        ++intervals;
  	      }
  	    }
        // Get next interval
        a2 = sdsl::util::next_bit(dict[source], b2+1);
        b2 = sdsl::util::next_bit(dict[source], a2+1);
      }
      swap(source, target);
      sdsl::util::set_to_value(dict[target], 0);
    }
    ++lcp_value;
  }
  return lcp;
}


}  // cdbg


#endif
