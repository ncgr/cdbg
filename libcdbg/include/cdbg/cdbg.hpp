#ifndef CDBG_HPP
#define CDBG_HPP

// std
#include <fstream>  // ifstream
#include <iomanip>  // setw
#include <iostream>  // cerr, endl, istream, ostream
#include <limits>  // numeric_limits
#include <stack>
#include <string>  // string, to_string
#include <tuple>
#include <utility>  // move
// sdsl
#include <sdsl/bit_vector_il.hpp>  // bit_vector_il
#include <sdsl/config.hpp>  // sdsl::conf, cache_config
#include <sdsl/construct.hpp>  // construct
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>  // read_member, write_member
#include <sdsl/structure_tree.hpp>  // structure_tree
#include <sdsl/util.hpp>  // sdsl::util
// local
#include "partial_lcp.hpp"


using std::cerr;
using std::endl;
using std::ifstream;
using std::istream;
using std::move;
using std::numeric_limits;
using std::ostream;
using std::setw;
using std::stack;
using std::string;
using std::to_string;
using std::tuple;
using std::vector;
using sdsl::bit_vector_il;
using sdsl::cache_config;
using sdsl::construct;
using sdsl::int_vector_buffer;
using sdsl::read_member;
using sdsl::structure_tree;
using sdsl::write_member;


namespace cdbg {


struct node
{
  uint64_t lb;
  uint64_t rb;
  uint64_t len;
  bool exit_node;
  vector<uint64_t> adj_list;
  vector<uint64_t> pos_list;
  node(uint64_t _lb=0, uint64_t _rb=0, uint64_t _len=0, bool _exit_node=false) :
    lb(_lb), rb(_rb), len(_len), exit_node(_exit_node) { }
};


template<
  class t_wt=wt_huff<bit_vector, rank_support_v<>, select_support_mcl<1>,
    select_support_mcl<0>>,
  class t_bv1=bit_vector,
  class t_bv3=bit_vector,
  class t_wt_doc=wt_huff<bit_vector, rank_support_v<>, select_support_scan<1>,
    select_support_scan<0>, int_tree<>>>
struct compressed_debruijn_graph
{
  typedef uint64_t size_type;

  private:
    struct node_c
    {
      uint64_t lb;  // lb of prefix kmer
      uint64_t len;
      uint64_t size;
      uint64_t first_lb;  // lb of suffix kmer
      node_c(
        uint64_t _lb=0,
        uint64_t _len=0,
        uint64_t _size=0,
        uint64_t _first_lb=0) :
        lb(_lb), len(_len), size(_size), first_lb(_first_lb) { }
    };

    uint64_t m_k;
    t_wt m_wt_bwt;
    vector<uint64_t> m_carray;
    vector<node_c> m_nodes;
    uint64_t m_right_max;
    vector<uint64_t> m_stop_nodes;
    t_bv1 m_bv1;
    t_bv3 m_bv3;
    typename t_bv1::rank_1_type m_bv1_rank;
    typename t_bv3::rank_1_type m_bv3_rank;
    t_wt_doc m_wt_doc;

    void create_carray()
    {
      m_carray = vector<uint64_t>(256, 0);
      for (uint64_t i = 0, sum = 0; i < 256; ++i) {
      	m_carray[i] = sum;
      	sum += m_wt_bwt.rank(m_wt_bwt.size(), i);
      }
    }

    void detect_nodes(cache_config& config)
    {
      // Create int_vector<2> that indicates if the lcp value is smaller, eqal
      // or greater than k
      int_vector<2> lcp_k = construct_partial_lcp<t_wt>(m_wt_bwt, m_carray, m_k);
      // It should be possible to stream lcp_k from file ... but memory peak is
      // during construct_partial_lcp
      stack<uint64_t> indexes;
      bit_vector bv1(m_wt_bwt.size(), 0);
      bit_vector bv3(m_wt_bwt.size(), 0);
      // Add right_maximal nodes
      bool open=false;
      uint64_t kvalue=0;
      uint64_t lb=0;
      uint64_t last_change=0;
      vector<uint64_t> lf = m_carray;
      int_vector_buffer<8> bwt(cache_file_name(sdsl::conf::KEY_BWT, config));
      for (uint64_t i = 1; i < lcp_k.size(); ++i) {
        ++lf[bwt[i-1]];
        if (lcp_k[i] == gt_k || lcp_k[i] == eq_k) {
          open = true;
          if (lcp_k[i] == eq_k) {
          	kvalue = i;
          }
        } else {
          if (open) {
            if (kvalue > lb) {
              bv1[lb] = true;
              bv1[i-1] = true;
              m_nodes.emplace_back(node_c(lb, m_k, i-lb, lb));
            }
            if (last_change > lb) {
              for (uint64_t j = lb; j <= i-1; ++j) {
                uint8_t c = bwt[j];
                bv3[lf[c]-1] = true;
              }
            }
            open = false;
          }
          lb = i;
        }
        if (bwt[i] != bwt[i-1] || bwt[i] <= 1) {
          last_change = i;
        }
      }
      if (open) {
        ++lf[bwt[lcp_k.size()-1]];
        if (kvalue > lb) {
          bv1[lb] = true;
          bv1[lcp_k.size()-1] = true;
          m_nodes.emplace_back(node_c(lb, m_k, lcp_k.size()-lb, lb));
        }
        if (last_change > lb) {
          for (uint64_t j = lb; j <= lcp_k.size()-1; ++j) {
            uint8_t c = bwt[j];
            bv3[lf[c]-1] = true;
          }
        }
      }
      // Add Endnodes
      for (uint64_t i = 0; i < m_carray[2]; ++i) {
        m_stop_nodes.emplace_back(m_nodes.size());
        m_nodes.emplace_back(node_c(i, 1, 1, i));
        bv3[i] = 0;
      }
      open = false;
      for (uint64_t i = 0; i < bv1.size(); ++i) {
        if (open) {
          bv3[i] = 0;
          if (bv1[i]) {
            open = false;
          }
        } else if (bv1[i]) {
          bv3[i] = 0;
          open = true;
        }
      }
      // Create rank support
      m_bv1 = t_bv1(move(bv1));
      m_bv3 = t_bv3(move(bv3));
      sdsl::util::init_support(m_bv1_rank, &m_bv1);
      sdsl::util::init_support(m_bv3_rank, &m_bv3);
    }

    void complete_nodes()
    {
      uint64_t quantity;
      vector<uint8_t> cs(m_wt_bwt.sigma);  // List of characters in the interval
      vector<uint64_t> rank_c_i(m_wt_bwt.sigma);  // Number of occurrence of character in [0 .. i-1]
      vector<uint64_t> rank_c_j(m_wt_bwt.sigma);  // Number of occurrence of character in [0 .. j-1]
      stack<uint64_t> order;
      for (uint64_t i = 0, undef = numeric_limits<uint64_t>::max();
      i < m_nodes.size(); ++i) {
        // Determine which node should be completed next
        uint64_t nodeid = i;
        if (i >= m_right_max) {
          nodeid = order.top();
          order.pop();
        }
        uint64_t cur_lb = m_nodes[nodeid].lb;
        uint64_t cur_rb = m_nodes[nodeid].lb+m_nodes[nodeid].size-1;
        uint64_t cur_len = m_nodes[nodeid].len;
        bool extend = true;
        while (extend) {
          extend = false;
          m_wt_bwt.interval_symbols(
            cur_lb,
            cur_rb+1,
            quantity,
            cs,
            rank_c_i,
            rank_c_j);
          for (uint64_t j = 0; j < quantity; ++j) {
            uint8_t c = cs[j];
            uint64_t lb = m_carray[c] + rank_c_i[j];
            uint64_t rb = m_carray[c] + rank_c_j[j] - 1;
            uint64_t ones = m_bv1_rank(lb+1);
            uint64_t node_number = undef;
            if (ones % 2 == 0 && m_bv1[lb] == 0) {
              // no-op
            } else {
              node_number = (ones-1)/2;
            }
            if (node_number != undef) {
              m_nodes[nodeid].lb = cur_lb;
              m_nodes[nodeid].len = cur_len;
            } else if (c <= 1) {  // c == sentinal
              m_nodes[nodeid].lb = cur_lb;
              m_nodes[nodeid].len = cur_len;
            } else {
              if (quantity == 1) {
                extend = true;
                cur_len++;
                cur_lb = lb;
                cur_rb = rb;
              } else {
                uint64_t next_node_id = m_right_max + m_bv3_rank(lb);
                order.push(next_node_id);
                m_nodes[next_node_id] = node_c(lb, m_k, rb-lb+1, lb);
                m_nodes[nodeid].lb = cur_lb;
                m_nodes[nodeid].len = cur_len;
              }
            }
          }
        }
      }
    }

  public:

    compressed_debruijn_graph() {}

    compressed_debruijn_graph(
      cache_config& config,
      uint64_t k,
      bool with_document_array) : m_k(k)
    {
      // Create WT of the BWT
      construct(m_wt_bwt, cache_file_name(sdsl::conf::KEY_BWT, config));
      // Create C-array (needed for interval_symbols)
      create_carray();
      // Detect and create nodes incl. bit vectors for calculation node numbers
      detect_nodes(config);
      // Add space for nodes not ending with an right maximal kmer
      m_right_max = m_nodes.size();
      uint64_t lmax = m_bv3_rank(m_bv3.size());
      m_nodes.resize(m_right_max + lmax);
      // Complete nodes
      complete_nodes();
      // Load Document Array
      if (with_document_array) {
        construct(m_wt_doc, cache_file_name("DA", config));
      }
    }

    tuple<vector<node>, vector<uint64_t>> get_explicit_representation() const
    {
      vector<node> graph(m_nodes.size());
      uint64_t d = m_carray[2];
      vector<uint64_t> start_nodes(d);
      uint64_t pos = m_wt_bwt.size()+1; // 1 indexed
      for (uint64_t s = 0, i = 0; s < d; ++s) {
        uint64_t prev_node_number = m_right_max-m_carray[2]+i;
        uint64_t idx = m_nodes[prev_node_number].lb;
        pos -= m_nodes[prev_node_number].len;
        graph[prev_node_number].len = m_nodes[prev_node_number].len;
        graph[prev_node_number].pos_list.emplace_back(pos);
        graph[prev_node_number].exit_node = true;
        auto res = m_wt_bwt.inverse_select(idx);
        i = m_carray[res.second] + res.first;
        while (res.second > 1) {  // c != sentinal
          uint64_t ones = m_bv1_rank(i+1);
          uint64_t node_number = 0;
          if (ones % 2 == 0 && m_bv1[i] == 0) {
            node_number = m_right_max + m_bv3_rank(i);
          } else {
            node_number = (ones-1)/2;
          }
          idx = m_nodes[node_number].lb + (i-m_nodes[node_number].first_lb);
          pos -= (m_nodes[node_number].len-m_k+1);
          // Store all information in new graph
          graph[node_number].adj_list.emplace_back(prev_node_number);
          graph[node_number].pos_list.emplace_back(pos);
          graph[node_number].len = m_nodes[node_number].len;
          // Go Node back
          prev_node_number = node_number;
          res = m_wt_bwt.inverse_select(idx);
          i = m_carray[res.second] + res.first;
        }
        start_nodes[d-1-s] = prev_node_number;
      }
      return make_tuple(move(graph), move(start_nodes));
    }

    void print_graph_statistics(ostream& out) const
    {
      // Print graph overview
      {
        uint64_t labels = 0;
        uint64_t edges = 0;
        uint64_t repeat_nodes = 0;
        uint64_t unique_nodes = 0;
        for (const auto& node : m_nodes) {
          labels += node.size;
          edges += node.size;
          if (node.size == 1) {
            ++unique_nodes;
          } else {
            ++repeat_nodes;
          }
        }
        edges -= m_stop_nodes.size();  // Stop nodes have no outgoing edge
        out << "       nodes=" << setw(10) << m_nodes.size() << endl;
        out << "unique nodes=" << setw(10) << unique_nodes << endl;
        out << "repeat nodes=" << setw(10) << repeat_nodes << endl;
        out << "      labels=" << setw(10) << labels << endl;
        out << "       edges=" << setw(10) << edges << endl;
        out << endl;
      }
      // Print lengths - statistics
      {
        uint64_t sum_length = 0;
        uint64_t sum_length_unique_nodes = 0;
        uint64_t sum_length_repeat_nodes = 0;
        uint64_t unique_nodes = 0;
        uint64_t repeat_nodes = 0;
        uint64_t max_length = 0;
        for (const auto& node: m_nodes) {
          sum_length += node.len;
          if (node.size == 1) {
            ++unique_nodes;
            sum_length_unique_nodes += node.len;
          } else {
            ++repeat_nodes;
            sum_length_repeat_nodes += node.len;
          }
          if (max_length < node.len) {
            max_length = node.len;
          }
        }
        out << "max_length=" << max_length << endl;
        out << "avg_length=" << sum_length/(double)m_nodes.size() << endl;
        out << "avg_length(unique_nodes)=" << sum_length_unique_nodes/(double)unique_nodes << endl;
        out << "avg_length(repeat_nodes)=" << sum_length_repeat_nodes/(double)repeat_nodes << endl;
        if (false) {
          vector<uint64_t> lengths(101, 0);
          uint64_t diff = (max_length+1)/(lengths.size()-1);
          if (diff == 0) {
            diff = 1;
          }
          for (const auto& node: m_nodes) {
            ++lengths[node.len/diff];
          }
          for (uint64_t i = 0; i < lengths.size(); ++i) {
            out << setw(10) << (i*diff) << " - " << setw(10) << ((i+1)*diff) << ": " << lengths[i] << " nodes" << endl;
          }
        }
        out << endl;
        out << "==========================================================================" << endl;
        out << endl;
      }
      // Print node-to-sequence relation details
      {
        uint64_t seq_number = m_wt_doc.sigma;
        // k = node_pair[i][j] means:
        // There are k nodes that contain (at least) sequence i and j
        vector<vector<uint64_t>> node_pair(seq_number, vector<uint64_t>(seq_number, 0));
        /// k = count_node_sequences[j] means:
        // There are k nodes that contain exactly j sequences
        vector<uint64_t> count_node_sequnces(seq_number+1, 0);
        uint64_t quantity;
        vector<uint64_t> cs(m_wt_doc.sigma);  // List of sequences in the interval
        vector<uint64_t> rank_c_i(m_wt_doc.sigma);  // Number of occurrence of character in [0 .. i-1]
        vector<uint64_t> rank_c_j(m_wt_doc.sigma);  // Number of occurrence of character in [0 .. j-1]
        for (const auto& node : m_nodes) {
          m_wt_doc.interval_symbols(node.lb, node.lb+node.size, quantity, cs, rank_c_i, rank_c_j);
          ++count_node_sequnces[quantity];
          for (uint64_t i = 0; i < quantity; ++i) {
            uint64_t seq_a = cs[i];
            for (uint64_t j = i; j < quantity; ++j) {
              uint64_t seq_b = cs[j];
              ++node_pair[seq_a][seq_b];
            }
          }
        }
        // "Sort"
        for (uint64_t i = 0; i < node_pair.size(); ++i) {
          for(uint64_t j = i+1; j < node_pair[i].size(); ++j) {
            node_pair[i][j] += node_pair[j][i];
            node_pair[j][i] = 0;
          }
        }
        // Print Results
        for (uint64_t i = 0; i < count_node_sequnces.size(); ++i) {
          out << setw(10) << count_node_sequnces[i] << " nodes covers exactly " << setw(3) << i << " sequences." << endl;
        }
        out << endl;
        for (uint64_t i = 0; i < node_pair.size(); ++i) {
          for (uint64_t j = 0; j < node_pair[i].size(); ++j) {
            out << setw(11) << node_pair[i][j];
          }
          out << endl;
        }
        out << endl;
      }
    }

    uint64_t get_k() const
    {
      return m_k;
    }


    vector<uint64_t> get_stop_nodes() const
    {
      return m_stop_nodes;
    }

    // Find all nodes that contains the pattern s
    tuple<vector<uint64_t>, uint64_t> find_nodes(const string& s) const
    {
      vector<uint64_t> result;
      assert(s.size() >= m_k);
      if (s.size() < m_k) {
        cerr << "Pattern length must at least " << m_k << endl;
        return make_tuple(result, -1);
      }
      // Find sa-interval
      uint64_t i = 0;
      uint64_t j = m_wt_bwt.size()-1;
      uint64_t pos = s.size();
      while (i < j && pos > s.size()-m_k) {
        pos--;
        uint8_t c = s[pos];
        i = m_carray[c] + m_wt_bwt.rank(i  , c);
        j = m_carray[c] + m_wt_bwt.rank(j+1, c)-1;
      }
      while (i <= j && pos > s.size()-m_k) {
        pos--;
        uint8_t c = s[pos];
        auto res = m_wt_bwt.inverse_select(i);
        if (res.second == c) {
          i = m_carray[c] + res.first;
          j = i;
        } else {
          i = j+1; // Not Found!
        }
      }
      if (i > j) {
        return make_tuple(result, -1);
      }
      uint64_t lb = i;
      uint64_t rb = j;
      // Find nodeid of end node (node that contains the suffix of length k)
      uint64_t nodeid = m_wt_bwt.size()+1;
      uint64_t l = 0;
      while (nodeid > m_wt_bwt.size()) {
        uint64_t ones = m_bv1_rank(i+1);
        if (ones % 2 == 1 || m_bv1[i] == 1) {
          nodeid = (ones-1)/2;
        } else {
          uint64_t ones_i = m_bv3_rank(i);
          uint64_t ones_j = m_bv3_rank(j+1);
          if (ones_i != ones_j) {
            nodeid = m_right_max + ones_i;
          } else {
            uint8_t c = 0;
            if (l < m_k) {
              c = s[s.size()-m_k+l];
            } else {
              while (i >= m_carray[c]) {
                ++c;
              }
              --c;
            }
            if (i == j) {
              i = m_wt_bwt.select(i-m_carray[c]+1, c); // ilf
              j = i;
            } else {
              i = m_wt_bwt.select(i-m_carray[c]+1, c); // ilf
              j = m_wt_bwt.select(j-m_carray[c]+1, c); // ilf
            }
            ++l;
            if (i < m_carray[2]) {  // Found sentinal
              nodeid = m_right_max-m_carray[2]+i;
              l -= (m_k-1);
            }
          }
        }
      }
      l = m_nodes[nodeid].len - l - m_k;  // Start position of suffix in current node
      // Add nodeid to path
      result.emplace_back(nodeid);
      // Find preceeding nodes
      i = lb;
      j = rb;
      while (i<=j && pos > 0) {
        // Add character to the front of the suffix
        pos--;
        uint8_t c = s[pos];
        i = m_carray[c] + m_wt_bwt.rank(i  , c);
        j = m_carray[c] + m_wt_bwt.rank(j+1, c)-1; // if i == j this can be done better with inverse_select!
        // auto res = m_wt_bwt.inverse_select(i);
        // Check if I'm in a new node
        if (l == 0) {
          uint64_t ones = m_bv1_rank(i+1);
          if (ones % 2 == 1 || m_bv1[i] == 1) {
            nodeid = (ones-1)/2;
          } else {
            nodeid = m_right_max + m_bv3_rank(i);
          }
          l = m_nodes[nodeid].len - m_k;
        } else {
          --l;
        }
        // Add nodeid to path
        result.emplace_back(nodeid);
      }
      // If pattern was not found
      if (i>j) {
        result.resize(0);
        l = 0;
      }
      return make_tuple(move(result), l);
    }

	// Find all sequences that occur in a node
	// Precondition: nodeid is valid
    vector<uint64_t> sequences_in_node(const uint64_t nodeid) const
    {
      assert(nodeid < m_nodes.size());
      auto node = m_nodes[nodeid];
      vector<uint64_t> result;
      uint64_t quantity;
      vector<uint64_t> cs(m_wt_doc.sigma);       // List of sequences in the interval
      vector<uint64_t> rank_c_i(m_wt_doc.sigma); // Number of occurrence of character in [0 .. i-1]
      vector<uint64_t> rank_c_j(m_wt_doc.sigma); // Number of occurrence of character in [0 .. j-1]
      m_wt_doc.interval_symbols(node.lb, node.lb+node.size, quantity, cs, rank_c_i, rank_c_j);
      for (uint64_t i = 0; i < quantity; ++i) {
        result.emplace_back(cs[i]);
      }
      return result;
    }


    template<class t_pod>
    uint64_t serialize_vector_pod(
      const vector<t_pod>& vpod,
      ostream& out,
      structure_tree_node* v,
      string name="") const
    {
      structure_tree_node* child = structure_tree::add_child(v, name, "vector<POD>");
      uint64_t written_bytes = 0;
      uint64_t size=vpod.size();
      written_bytes += write_member(size, out);
      out.write((char*)vpod.data(), vpod.size()*sizeof(vpod[0]));
      written_bytes += vpod.size()*sizeof(vpod[0]);
      structure_tree::add_size(child, written_bytes);
      return written_bytes;
    }

    template<class t_pod> void load_vpod(vector<t_pod>& vpod, istream& in)
    {
      uint64_t size;
      read_member(size, in);
      vpod.resize(size);
      in.read((char*)vpod.data(), size*sizeof(vpod[0]));
    }

    //! Serialize into stream
    size_type serialize(
      ostream& out,
      structure_tree_node* v=nullptr,
      string name="") const
    {
      structure_tree_node* child = structure_tree::add_child(v, name, sdsl::util::class_name(*this));
      size_type written_bytes = 0;
      written_bytes += write_member(m_k, out, child, "k");
      written_bytes += m_wt_bwt.serialize(out, child, "wt_bwt");
      written_bytes += serialize_vector_pod(m_carray, out, child, "c_array");
      written_bytes += serialize_vector_pod(m_nodes, out, child, "nodes");
      written_bytes += write_member(m_right_max, out, child, "right_max");
      written_bytes += serialize_vector_pod(m_stop_nodes, out, child, "stop_nodes");
      written_bytes += m_bv1.serialize(out, child, "bv1");
      written_bytes += m_bv3.serialize(out, child, "bv3");
      written_bytes += m_bv1_rank.serialize(out, child, "bv1_rank");
      written_bytes += m_bv3_rank.serialize(out, child, "bv3_rank");
      written_bytes += m_wt_doc.serialize(out, child, "wt_doc");
      structure_tree::add_size(child, written_bytes);
      return written_bytes;
    }

    //! Load sampling from disk
    void load(istream& in)
    {
      read_member(m_k, in);
      m_wt_bwt.load(in);
      load_vpod(m_carray, in);
      load_vpod(m_nodes, in);
      read_member(m_right_max, in);
      load_vpod(m_stop_nodes, in);
      m_bv1.load(in);
      m_bv3.load(in);
      m_bv1_rank.load(in, &m_bv1);
      m_bv3_rank.load(in, &m_bv3);
      m_wt_doc.load(in);
    }
};


// WTTYPE
#ifdef WTBV
  // bv
  typedef wt_huff<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>> WT_TYPE;
#elif WTBV5
  // bv5
  typedef wt_huff<bit_vector, rank_support_v5<>, select_support_mcl<1>, select_support_mcl<0>> WT_TYPE;
#elif WTIL
  // bvil
  typedef wt_huff<bit_vector_il<>, rank_support_il<>, select_support_il<1>, select_support_il<0>> WT_TYPE;
#elif WTRRR15
  // rrr15
  typedef wt_huff<rrr_vector<15>, rrr_vector<15>::rank_1_type, rrr_vector<15>::select_1_type, rrr_vector<15>::select_0_type> WT_TYPE;
#elif WTRRR63
  // rrr63
  typedef wt_huff<rrr_vector<63>, rrr_vector<63>::rank_1_type, rrr_vector<63>::select_1_type, rrr_vector<63>::select_0_type> WT_TYPE;
#elif WTSD
  // sd
  typedef wt_huff<sd_vector<>, sd_vector<>::rank_1_type, sd_vector<>::select_1_type> WT_TYPE;
#else
  // bv
  typedef wt_huff<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>> WT_TYPE;
#endif


// BV1
#ifdef BV1BV
  // bv
  typedef bit_vector BV1_TYPE;
#elif BV1IL
  // bvil
  typedef bit_vector_il<> BV1_TYPE;
#elif BV1RRR15
  // rrr15
  typedef rrr_vector<15> BV1_TYPE;
#elif BV1RRR63
  // rrr63
  typedef rrr_vector<63> BV1_TYPE;
#elif BV1SD
  // sd
  typedef sd_vector<> BV1_TYPE;
#else
  // bv
  typedef bit_vector BV1_TYPE;
#endif


// BV3
#ifdef BV3BV
  // bv
  typedef bit_vector BV3_TYPE;
#elif BV3IL
  // bvil
  typedef bit_vector_il<> BV3_TYPE;
#elif BV3RRR15
  // rrr15
  typedef rrr_vector<15> BV3_TYPE;
#elif BV3RRR63
  // rrr63
  typedef rrr_vector<63> BV3_TYPE;
#elif BV3SD
  // sd
  typedef sd_vector<> BV3_TYPE;
#else
  // bv
  typedef bit_vector BV3_TYPE;
#endif


typedef compressed_debruijn_graph<WT_TYPE, BV1_TYPE, BV3_TYPE> CDBG;


}  // cdbg


#endif
