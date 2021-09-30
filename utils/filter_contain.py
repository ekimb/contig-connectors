from pafpy import PafFile
import argparse
import sys


def main(argsj):
    pi = args.paf
    of = args.tblast

    def frac_covered(r):
        return max(r.query_coverage, r.target_coverage)

    ofile = open(of, 'w')
    ofile.write('qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n')
    with PafFile(pi) as paf:
        for record in paf:
            if record.qname != record.tname:
                if record.qlen >= args.minlen:
                    if (frac_covered(record) >= args.fract) and (record.blast_identity() >= args.ident):
                        print(record)
                        ofile.write(f"{record.qname}\t{record.tname}\t{record.blast_identity()}\t*\t*\t*\t*\t*\t*\t*\t*\t*\n")
                    else:
                        print(record.qname, record.tname, frac_covered(record), record.blast_identity())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter minimap2 output to find high-quality containment.')
    parser.add_argument('--paf', help='input PAF file')
    parser.add_argument('--tblast', help='output tab-delimited BLAST file')
    parser.add_argument('--minlen', default=500, type=int, help='minimum contig length to consider')
    parser.add_argument('--ident', default=0.95, type=float, help='required percent identity')
    parser.add_argument('--fract', default=0.95, type=float, help='fraction of the contained contig that must be covered by the match')
    args = parser.parse_args()
    main(args)
