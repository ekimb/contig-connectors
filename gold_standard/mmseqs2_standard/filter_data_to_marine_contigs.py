import pandas as pd

QUERY_SAMPLE = 'nmdc:mga04781'

def make_S3_path(sample):
    return 's3://psss-metagenomics-codeathon-data/marine/%s/assembly/%s_contigs.fna' % (sample, sample.replace(':', '_'))

open('query_paths.txt', 'w').write(make_S3_path(QUERY_SAMPLE) + '\n')

metadata = pd.read_csv('data.tsv', sep='\t')
metadata = metadata.query("TAXA == 'marine metagenome'")
metadata = metadata.query("MGA_ID != @QUERY_SAMPLE")

full_paths = [make_S3_path(i) for i in metadata['MGA_ID'].values]
open('reference_paths.txt', 'w').write('%s\n' % '\n'.join(full_paths))
