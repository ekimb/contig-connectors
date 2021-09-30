from pafpy import PafFile
import argparse
import sys

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

	missed = len(truth - pred)
	fp = len(pred - truth)
	tp = len(truth & pred)
	print(f"tp: {tp}, fn: {missed}, fp: {fp}")
	print("precision : {}".format(tp / (tp+fp)))
	print("recall : {}".format(tp / (tp+missed)))

    """
	d = {}
	with open(args.lens) as ifile:
		for l in ifile:
			toks = l.rstrip().split()
			d[toks[0]] = int(toks[1])

	fn =  list(truth - pred)
	num = 0
	for f in fp:
		if num > 10:
			break
		l1 = d[f[0]]
		l2 = d[f[1]]
		if l1 >= 500 and l2 >= 500:
			print(f"{f[0]} : {l1} <=> {f[1]} {l2}")	
			num += 1
	"""

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='compare to ground truth containments.')
	parser.add_argument('--truth', help='truth TAB6 BLAST file')
	parser.add_argument('--pred', help='predicted TAB6 BLAST file')
	parser.add_argument('--lens', help='lengths file')
	args = parser.parse_args()
	main(args)
