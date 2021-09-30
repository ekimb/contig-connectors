#!/usr/bin/env python3

"""
filter a BLAST6 output file to contain only hits
that meet our definition for containment.
"""

import argparse
import sys

def main(args):
    pi = args.input
    of = args.out

    ld = {}
    with open(args.lens, 'r') as ifile:
        for l in ifile:
            toks = l.rstrip().split()
            ld[toks[0]] = int(toks[1])
    
    def frac_covered(r):
        return float(r.qend - r.qstart) / r.qlen

    ofile = open(of, 'w')
    ofile.write('qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n')
    with open(pi) as ifile:
        for record in ifile:
            toks = record.rstrip().split()
            if toks[0] != toks[1]:
                if ld[toks[0]] >= 500 and ld[toks[1]] >= 500:
                    aln_len = int(toks[3])
                    fc = aln_len / float(min(ld[toks[0]], ld[toks[1]]))
                    if (fc >= 0.95) and (float(toks[2]) >= 95.0):
                        ofile.write(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter minimap2 output to find high-quality containment.')
    parser.add_argument('--input', help='input blast 6')
    parser.add_argument('--out', help='output tab-delimited BLAST file')
    parser.add_argument('--lens', help='file with contig lengths')
    args = parser.parse_args()
    main(args)
