In [1], the authors present an efficient algorithm for constructing a
representation of a compressed de Bruijn graph for pangenome analysis that
enables search.
**cdbg** is a refactoring of [1]'s code (original code available at
<https://www.uni-ulm.de/in/theo/research/seqana.html>) that breaks it into a
dynamic library and a command line program.
The dynamic library provides developers the graph data structure and an API for
reading/writing the structure from disk, and the command line program provides
users an interface for constructing and interacting with compressed de Bruijn
graphs.

The goal of this work is to enhance the adoptability of the compressed de Bruijn
graph data model from [1].
By providing a program for efficiently constructing compressed de Bruijn graphs
and a library for utilizing them, developers are relieved of implementing their
own de Bruijn graph data structure and construction algorithms while allowing
users to perform multiple analyses on a single graph.

## Algorithm and Representation
The graph construction algorithm uses an FM-index to build a compressed de
Bruijn graph directly (rather than generating all possible k-mers) from a set of
sequences in O(n log &sigma;) time, where n is the total length of the sequences
and &sigma; is the size of the alphabet.

The program outputs an _implicit_ representation of the graph that is linear in
the number of graph nodes and is not effected by the number of graph edges or
k-mer occurrences.
This is advantageous because the representation is lossless and is likely much
smaller than the total length of the sequences (sub-linear).
Additionally, the implicit representation supports searching the graph for
sequences of length geq k in worst-case O((m+l) log &sigma;) time, where m is
the length of the search sequence and l is the length of the longest compressed
k-mer chain in the graph.

## Compiling
TODO


## Usage
TODO


## References
[1]: Beller, Timo, and Enno Ohlebusch. "A representation of a compressed de Bruijn graph for pan-genome analysis that enables search." _Algorithms for Molecular Biology_ 11.1 (2016): 20.
