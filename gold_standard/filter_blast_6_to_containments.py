"""Filter a blast out format 6/8 file to only containments"""
import argparse
import pandas as pd

FAI_HEADER = ['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH', 'QUALOFFSET']
QUERY_INDEX_PATH = 'query/query.ml500.fna.fai'
REFERENCE_INDEX_PATH = 'reference/reference.ml500.fna.fai'

BOUTFMT6_COLUMNS = ['qId', 'tId', 'seqIdentity', 'alnLen', 'mismatchCnt', 'gapOpenCnt', 'qStart', 'qEnd', 'tStart',
                    'tEnd', 'eVal', 'bitScore']
RESULTS_PATH = 'result.b6'


PERCENT_IDENTITY = 95
COVERED_LENGTH = .95

OUTPUT_PATH = 'results.containments.b6'


def check_containment(row, query_index, reference_index, percent_identity=PERCENT_IDENTITY, covered_length=COVERED_LENGTH):
    """Checks if a row from a blast out format 6 file is a containment
    Takes in a row from a blast out format 6 table, a DataFrames with query sequence and reference sequence data.
    """
    if (row['qId'] != row['tId']) and (row['seqIdentity'] >= percent_identity):
        query_covered = row['alnLen']/float(query_index.loc[row['qId'], 'LENGTH'])
        reference_covered = row['alnLen']/float(reference_index.loc[row['tId'], 'LENGTH'])
        if query_covered >= covered_length or reference_covered >= covered_length:
            return True
        else:
            return False
    else:
        return False


def find_containments(results_path, query_index_path, reference_index_path, output_path,
                      percent_identity=PERCENT_IDENTITY, covered_length=COVERED_LENGTH,
                      as_decimal=False, has_header=False):
    if as_decimal:
        percent_identity = percent_identity/100

    query_index = pd.read_csv(query_index_path, sep='\t', header=None, names=FAI_HEADER)
    query_index = query_index.set_index('NAME')
    reference_index = pd.read_csv(reference_index_path, sep='\t', header=None, names=FAI_HEADER)
    reference_index = reference_index.set_index('NAME')


    if has_header:
        search_results = pd.read_csv(results_path, sep='\t', header=True)
        search_results.columns = BOUTFMT6_COLUMNS
    else:
        search_results = pd.read_csv(results_path, sep='\t', header=None, names=BOUTFMT6_COLUMNS)

    is_containment = search_results.apply(check_containment, axis=1,
                                          query_index=query_index, reference_index=reference_index,
                                          percent_identity=percent_identity, covered_length=covered_length)

    search_results.loc[is_containment].to_csv(output_path, sep='\t', index=None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--results_path', type=str, help='BLAST format 6 table of alignments to check for containments.')
    parser.add_argument('-q', '--query_index_path', type=str, help='Index of query sequences in samtools faidx format.')
    parser.add_argument('-r', '--reference_index_path', type=str, help='Index of reference sequences in samtools faidx format.')
    parser.add_argument('-o', '--output_path', type=str, help='Path to write containments in BLAST format 6.')
    parser.add_argument('--percent_identity', type=float, default=PERCENT_IDENTITY,
                        help='Minimum percent identity to consider an alignment to be a containment. [0-100]')
    parser.add_argument('--covered_length', type=float, default=COVERED_LENGTH,
                        help='Minimum percent covered of shorter contig to consider an alignment to be a containment. [0-1]')
    parser.add_argument('--as_decimal', action='store_true', default=False,
                        help="Treat percent_identity as decimal. Divides percent_identity by 100.")
    parser.add_argument('--has_header', action='store_true', default=False,
                        help="Input BLAST format 6 table has headers. Header will be changed to standard column names.")

    args = parser.parse_args()
    find_containments(args.results_path, args.query_index_path, args.reference_index_path, args.output_path, args.percent_identity,
                      args.covered_length, args.as_decimal, args.has_header)
