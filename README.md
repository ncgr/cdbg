**cdbg** (compressed de Bruijn graph) is a refactoring of the code presented in
[1] (original code available at
<https://www.uni-ulm.de/in/theo/research/seqana.html>).
Specifically, it breaks the code into a dynamic library and a command line
program.
The dynamic library provides developers the graph data structure and an API for
reading/writing the structure from disk, and the command line program provides
users an interface for constructing and interacting with compressed de Bruijn
graphs.

The goal of this work is to enhance the adoptability of the compressed de Bruijn
graph data model from [1].
By providing a program for efficiently constructing compressed de Bruijn graphs
and a library for utilizing them, developers are relieved of implementing their
own de Bruijn graph data structure and construction algorithms while enabling
users to perform multiple analyses on a single graph.

## Compiling
TODO


## Usage
TODO


## References
[1]: Beller, Timo, and Enno Ohlebusch. "A representation of a compressed de Bruijn graph for pan-genome analysis that enables search." _Algorithms for Molecular Biology_ 11.1 (2016): 20.
