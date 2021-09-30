from pafpy import PafFile
import argparse
#import pandas as pd
import sys

"""
Filters a PAF file for containment according to our
criteria.

NOTE: this should ideally be merged with the script to 
filter a BLAST6 file for containment.
"""


#FAI_HEADER = ['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH', 'QUALOFFSET']

def write_filtered_paf(pi, ofile):
    def frac_covered(r):
        return max(r.query_coverage, r.target_coverage)

    with PafFile(pi) as paf:
        for record in paf:
            if record.qname != record.tname:
                if record.qlen >= args.minlen:
                    if (frac_covered(record) >= args.fract) and (record.blast_identity() >= args.ident):
                        ofile.write(f"{record.qname}\t{record.tname}\t{record.blast_identity()}\t*\t*\t*\t*\t*\t*\t*\t*\t*\n")

def main(argsj):
    pi = args.paf
    of = args.tblast

    ofile = open(of, 'w')
    ofile.write('qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n')
    write_filtered_paf(pi, ofile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter minimap2 output to find high-quality containment.')
    parser.add_argument('--paf', help='input PAF file')
    parser.add_argument('--tblast', help='output tab-delimited BLAST file')
    parser.add_argument('--minlen', default=500, type=int, help='minimum contig length to consider')
    parser.add_argument('--ident', default=0.95, type=float, help='required percent identity')
    parser.add_argument('--fract', default=0.95, type=float, help='fraction of the contained contig that must be covered by the match')
    args = parser.parse_args()
    main(args)
