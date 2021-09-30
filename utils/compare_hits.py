#!/usr/bin/env python3

import argparse
import sys

"""
Parse two BLAST6 format files and compute the precision 
and recall of a predicted set of containments against 
the truth.
"""

#FAI_HEADER = ['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH', 'QUALOFFSET']

def parse_tab6(tf):
	res = []
	first = True
	with open(tf, 'r') as ifile:
		for l in ifile:
			if first:
				first = False
				continue
			toks = l.rstrip().split('\t')
			if toks[0] != toks[1]:
				res.append(tuple(sorted((toks[0], toks[1]))))
	return res

def main(args):
	tf = args.truth
	pf = args.pred
	truth = set(parse_tab6(tf))
	pred = set(parse_tab6(pf))

	fn = len(truth - pred)
	fp = len(pred - truth)
	tp = len(truth & pred)
	print(f"tp: {tp}, fn: {fn}, fp: {fp}")
	print("precision : {}".format(tp / (tp+fp)))
	print("recall : {}".format(tp / (tp+fn)))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='compare to ground truth containments.')
	parser.add_argument('--truth', help='truth TAB6 BLAST file')
	parser.add_argument('--pred', help='predicted TAB6 BLAST file')
	#parser.add_argument('--fai', help='faidx file (to get lengths)')
	args = parser.parse_args()
	main(args)
