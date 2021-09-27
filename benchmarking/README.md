## Code for the benchmarking experiments

---

### Notes

**Day 1:**
Write code that takes something similary formatted to the gold standard and see how it goes. What is the input and output? Clearly define them and let other teams know. Containment is just one of many possible relationships that can happen? Consider workflow system for running benchmarking experiments, *e.g.* Snakemake. Have an output and compare it to the gold standard. How to include potential novel implementation, *e.g.* Docker container.

**Snakemake workflow specification**
- Inputs:
    - File with ground truth data (BLAST tabular output, may omit some columns) - [docs](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6)
    - Query contig (FASTA)
    - Reference data (FASTA)
- Steps:
    - Running novel implementation
    - Post process output of novel implementation
    - Generate report
        - FP, FN
        - Precision & recall
        - Performance metrics (Runtime & Memory usage)
        - Score ?
    - Store results in database
    - MongoDB / TinyDB ?
- Output:
    - Report with accuracy and performance metrics

