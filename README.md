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
