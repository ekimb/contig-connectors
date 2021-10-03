import pandas as pd

QUERY_SAMPLE = 'nmdc:mga04781'
QUERY_AWS_PATH = 's3://psss-metagenomics-codeathon-data/marine/nmdc:mga04781/assembly/nmdc_mga04781_contigs.fna'
metadata = pd.read_csv('data.tsv', sep='\t')
metadata = metadata.query("MGA_ID != @QUERY_SAMPLE")
metadata = metadata.set_index('MGA_ID')


open('query_paths.txt', 'w').write(QUERY_AWS_PATH + '\n')

full_paths = ['%s/assembly/%s_contigs.fna' % (aws_path, sample.replace(':', '_'))
              for sample, aws_path in metadata['AWS_PATH'].iteritems()]
open('reference_paths.txt', 'w').write('%s\n' % '\n'.join(full_paths))
