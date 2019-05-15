// std
#include <iomanip>  // setw
#include <iostream>  // endl, left;
#include <string>
// GNU
#include <getopt.h>  // getopt_long, no_argument, option, required_argument
// local
#include "commands/construct.hpp"
#include "commands/find_pattern.hpp"
#include "commands/impl2expl.hpp"
#include "commands/print_graph_details.hpp"

using std::cerr;
using std::endl;
using std::left;
using std::setw;
using std::string;


struct options_t {
  string inputfile;
  string outputfile;
  string kfile;
  string graphfile;
  string patternfile;
};


void print_command(const string& command, const string& description)
{
  cerr << setw(25) << left << "   " << command << " - " << description << endl;
}


void print_option(const string& arg, const string& description)
{
  cerr << setw(10) << arg << " " << description << endl;
}


void usage(const string& program, const string& command="")
{
  cerr << "usage: " << program << " ";
  if (command == "") {
    cerr << "<command> <args>" << endl;
    cerr << endl;
    cerr << "Where command must be one of the following list:" << endl;
    print_command("construct", " - Construct the compressed de bruijn graph");
    print_command("print_graph_details", " - Print graph details");
    print_command("find_pattern", " - Finding pattern in the pan-genome");
    print_command("impl2expl", " - Convert to explicit representation");
  } else {
    cerr << command << " options" << endl;
    cerr << endl;
    if (command == "construct") {
      print_option("-i, --inputfile=INFILE", "the input file");
      print_option("-o, --outputfile=OUTFILE", "the output file");
      print_option("-k, --kfile=KFILE", "text file containing k values");
    } else if(command == "print_graph_details") {
      print_option("-g, --graphfile=GRAPHFILE", "graph file, created via construct command");
    } else if(command == "find_pattern") {
      print_option("-g, --graphfile=GRAPHFILE", " graph file, created via construct command");
      print_option("-p, --patternfile=PATTERNFILE", " pattern file, containing pattern");
    } else if(command == "impl2expl") {
      cerr << "Program will create OUTFILE.dot and OUTFILE.start_nodes.txt" << endl;
      cerr << endl;
      print_option("-g, --graphfile=GRAPHFILE", " graph file, created via construct command");
      print_option("-o, --outputfile=OUTFILE", " the output file");
    }
  }
  cerr << endl;
}


void check_argument_given(
  const string& program,
  const string& command,
  const string& value,
  const string name)
{
  if (value == "") {
    usage(program, command);
    cerr << "ERROR: No " << name << " given." << endl;
    exit(EXIT_FAILURE);
  }
}


void call_construct(const string& program, const options_t& opts)
{
  check_argument_given(program, "construct", opts.inputfile, "inputfile");
  check_argument_given(program, "construct", opts.outputfile, "outputfile");
  check_argument_given(program, "construct", opts.kfile, "kfile");
  cdbg::commands::construct(opts.inputfile, opts.outputfile, opts.kfile, true);
}


void call_print_graph_details(const string& program, const options_t& opts)
{
  check_argument_given(program, "print_graph_details", opts.graphfile, "graphfile");
  cdbg::commands::print_graph_details(opts.graphfile);
}


void call_find_pattern(const string& program, const options_t& opts)
{
  check_argument_given(program, "find_pattern", opts.graphfile, "graphfile");
  check_argument_given(program, "find_pattern", opts.patternfile, "patternfile");
  cdbg::commands::find_pattern(opts.graphfile, opts.patternfile);
}


void call_impl2expl(const string& program, const options_t& opts)
{
  check_argument_given(program, "impl2expl", opts.graphfile, "graphfile");
  check_argument_given(program, "impl2expl", opts.outputfile, "outputfile");
  if (!cdbg::commands::impl2expl(opts.graphfile, opts.outputfile)) {
    exit(1);
  }
}


options_t parse_args(int argc, char* argv[])
{
  options_t opts;
  const char* const short_opts = "i:o:k:g:p:h";
  static struct option long_opts[] =
  {
    {"inputfile", required_argument, nullptr, 'i'},
    {"outputfile", required_argument, nullptr, 'o'},
    {"kfile", required_argument, nullptr, 'k'},
    {"graphfile", required_argument, nullptr, 'g'},
    {"patternfile", required_argument, nullptr, 'p'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, no_argument, nullptr, 0}
  };
  int64_t c;
  while ((c = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
    switch (c) {
      case 'i':
        opts.inputfile = string(optarg);
        break;
      case 'o':
        opts.outputfile = string(optarg);
        break;
      case 'k':
        opts.kfile = string(optarg);
        break;
      case 'g':
        opts.graphfile = string(optarg);
        break;
      case 'p':
        opts.patternfile = string(optarg);
        break;
      default:
        usage(argv[0], argv[1]);
        break;
    }
  }
  return opts;
}


int main(int argc, char* argv[]) {
  if (argc <= 1) {
    usage(argv[0]);
    return 1;
  }
  string command = string(argv[1]);
  options_t opts = parse_args(argc, argv);
  if (command == "construct") {
    call_construct(argv[0], opts);
  } else if (command == "print_graph_details") {
    call_print_graph_details(argv[0], opts);
  } else if (command == "find_pattern") {
    call_find_pattern(argv[0], opts);
  } else if(command == "impl2expl") {
    call_impl2expl(argv[0], opts);
  } else {
    usage(argv[0], command);
    return 1;
  }
  return 0;
}
