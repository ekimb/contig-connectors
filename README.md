
# Details about current implementation code and organization

### Search

The current implementation takes a sketching-based approach for containment query.  

`indexing` : Briefly, ordered minimizer sketches are computed for each input sequence with the specified k-mer length and density.
The sketches are serialized to file.  This indexer code is written in Rust and exists in the `rust-cc`
subdirectory.  To build this code, change into the `rust-cc` directory and run `cargo build --release`.

`query` : The query code reads the indexed sketches and bins sketches by the hashes they contain.  The query is sketched and checked against 
the bins that may have a hit for the query.  Within the bin, the minimum Hamming distance is computed against the bin elements (brute force for now but this can easily be improved). 
Queries are finially filtered by (gapless) ANI and output.  The query code is written in C++ and lives in the `src` subdirectory.  To build this code create a `build` directory,
change into it, and execute `cmake .. && cmake --build . --config release`.

### Utilities 

The `utils` directory contains a few utilities to help with formatting the output of different tools as well as comparing predictions to ground truth containments.  Currently, 
there are 3 scripts in this directory:

`filter_contain.py` : This takes input from either MashMap format or PAF format (e.g. from minimap2) and filters the hits for containment as determined
by our definition (at least 95% of the query is covered by an alignment of at least 95% sequence identity).

`filter_blast.py` : This should be merged with the above for functionality, but this currently filters results in BLAST6 format to retain just hits that
match the adopted definition of containment.

`compare_hits.py` : From a ground truth and predicted set of hits in BLAST6 format (with at least the first 3 columns populated), this script computes precision and 
and recall of the predicted containments compared to the truth.  This comparison ignores order (i.e. p contained in q is treated identically to q contained in p) and 
explicitly filters out any reported self containments.
=======

# Benchmarking Metagenomic Contig Matches
Last modified on: September 30, 2021

This project is part of the Petabyte-Scale Sequence Search: Metagenomics Benchmarking Codeathon, hosted virtually from Monday, September 27, 2021 to Friday, October 1, 2021. 

## Problem Statement

Given a query contig, can we find contigs in other samples that are completely contained within it, or that completely contain it?

This is biologically relevant for instances when you have interesting metagenome contigs and want to determine if the same contigs have been observed in any other metagenome datasets. 

Within the codeathon, we identified a way to benchmark contig containments (ie. if you use different tools to identify containments, how close to the "truth" are the results). Detailed codeathon project organization is found in the [wiki](https://github.com/NCBI-Codeathons/psss-team2/wiki).


## Team

Mihai Pop, PhD- 
University of Maryland

Rob Patro, PhD- 
University of Maryland

Jackie Michaelis, PhD- 
University of Maryland

Nicholas Cooley- 
University of Pittsburgh

Barış Ekim- 
Massachusetts Institute of Technology (MIT)

Priyanka Ghosh- 
National Institutes of Health (NIH)

Harihara Subrahmaniam Muralidharan- 
University of Maryland

Amatur Rahman- 
Pennsylvania State University

Vinicius Salazar-
University of Melbourne, Australia

Michael Shaffer- 
Colorado State University

Andrew Tritt- 
Lawrence Berkeley National Laboratory
