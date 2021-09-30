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
mash_header_list = ['qId', 'qlen', 'qstart', 'qend', 'strand', 'tId', 'tlen', 'tstart', 'tend', 'ani']
mash_header = dict(zip(mash_header_list, range(len(mash_header_list))))

def write_filtered_paf(pi, ofile):
    def frac_covered(r):
        return max(r.query_coverage, r.target_coverage)

    with PafFile(pi) as paf:
        for record in paf:
            if record.qname != record.tname:
                if record.qlen >= args.minlen:
                    if (frac_covered(record) >= args.fract) and (record.blast_identity() >= args.ident):
                        ofile.write(f"{record.qname}\t{record.tname}\t{record.blast_identity()}\t*\t*\t*\t*\t*\t*\t*\t*\t*\n")

def write_filtered_mash(pi, ofile):
    def frac_covered(r):
        qlen = int(r[mash_header['qlen']])
        tlen = int(r[mash_header['tlen']])

        qspan = int(r[mash_header['qend']]) - int(r[mash_header['qstart']]) 
        tspan = int(r[mash_header['tend']]) - int(r[mash_header['tstart']]) 
        
        query_cov = float(qspan) / qlen
        target_cov = float(tspan) / tlen
        return max(query_cov, target_cov)

    with open(pi, 'r') as ifile:
        for l in ifile:
            toks = l.rstrip().split()
            if toks[mash_header['qId']] != toks[mash_header['tId']]:
                if int(toks[mash_header['qlen']]) >= args.minlen:
                    ani = (float(toks[mash_header['ani']]) / 100.0)
                    if (frac_covered(toks) >= args.fract) and (ani >= args.ident):
                        ofile.write(f"{toks[mash_header['qId']]}\t{toks[mash_header['tId']]}\t{ani}\t*\t*\t*\t*\t*\t*\t*\t*\t*\n")

def main(args):
    if args.fmt == "paf":
        pi = args.input
        of = args.output
        ofile = open(of, 'w')
        ofile.write('qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n')
        write_filtered_paf(pi, ofile)
    elif args.fmt == "mash":
        pi = args.input
        of = args.output
        ofile = open(of, 'w')
        ofile.write('qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n')
        write_filtered_mash(pi, ofile)
    else:
        print(f"ERROR: don't know input format {args.fmt}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter minimap2 output to find high-quality containment.')
    parser.add_argument('--input', help='input file')
    parser.add_argument('--output', help='output tab-delimited BLAST file')
    parser.add_argument('--fmt', help="input format {paf, mash}")
    parser.add_argument('--minlen', default=500, type=int, help='minimum contig length to consider')
    parser.add_argument('--ident', default=0.95, type=float, help='required percent identity')
    parser.add_argument('--fract', default=0.95, type=float, help='fraction of the contained contig that must be covered by the match')
    args = parser.parse_args()
    main(args)
