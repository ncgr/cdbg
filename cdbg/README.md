cdbg is a refactoring of the "a4" algorithm from _A representation of a
compressed de Bruijn graph for pan-genome analysis that enables search_.
The tool builds an implicit representation of a  compressed de Bruijn graph
(a Burrows-Wheeler transform with auxiliary data structures) for a given FASTA
file and _k_-mer size and saves it as a ` .bin` file.
The implicit representation is more space efficient than an explicit de Bruijn
graph and it enables efficient sequence searching.

Though this is a refactoring of an existing tool, the existing functionality
will be revised and new functionality will be added.
So the tool may contain bugs and the instructions on how to use it are subject
to change.

## Compiling

cdbg is implemented in C++11.
The `Makefile` is configured to use g++ and assumes it is in your path.

cdbg's only dependencies are the
[Succinct Data Structure Library](https://github.com/simongog/sdsl-lite)
(sdsl-lite) and our own libcdbg, both of which can be easily installed on any
Unix style system.

Once sdsl-lite and libcdbg are installed, you can compile cdbg with the `make`
command.
This will create an executable `cdbg` file in the `bin/` directory.


## Usage

First, construct the implicit representation of the compressed de Bruijn graph:
```
./cdbg construct --inputfile=input.fa --outputfile=example --kfile=kfile.txt
```
where `kfile.txt` is a file containing a single positive integer per line for
each _k_-mer size you want to generate a graph for.
For instance, if `kfile.txt` contains a single line with the value `100`, then
the previous command will create a `example.k100.bin` file.

To see graph statistics use:
```
./cdbg print_graph_details --graphfile=example.k100.bin
```

Search the graph for a sequence as follows:
```
./cdbg find_pattern --graphfile=example.k100.bin --patternfile=pattern.txt
```
where `pattern.txt` is a file containing a single sequence per line with length
greater than or equal to the _k_-mer size of the graph you're searching.

Generate an explicit representation (`.dot` file) from the implicit
representation as follows:
```
./cdbg impl2expl --graphfile=example.k100.bin --outputfile=example.k100
```
This will create a `example.k100.dot` file and a `example.k100.start_nodes.txt`
file.
