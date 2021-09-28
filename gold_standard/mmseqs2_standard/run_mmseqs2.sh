# set up conda
wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
export MINICONDA_PREFIX="$HOME/miniconda"
bash miniconda.sh -b -p $MINICONDA_PREFIX
export PATH="$MINICONDA_PREFIX/bin:$PATH"
conda config --set always_yes yes
conda update -q conda
conda config --add channels bioconda
conda config --add channels conda-forge
conda info -a

conda install wget pandas mmseqs2 bbmap

# get data sheet
wget https://raw.githubusercontent.com/NCBI-Codeathons/psss-datasets/master/data.tsv

# script filters data down to files with S3 paths for all references and query
# TODO: update path to main repo or switch to cloning repo instead of pulling specific file
wget https://raw.githubusercontent.com/shafferm/psss-team2/main/gold_standard/mmseqs2_standard/filter_data_to_marine_contigs.py
python filter_data_to_marine_contigs.py

# download all data
mkdir query
while read query_path; do
    aws s3 cp $query_path query
done < query_paths.txt
cat query/*.fna > query/query.fna
reformat.sh in=query/query.fna out=query/query.ml500.fna ml=500

mkdir reference
while read reference_path; do
    aws s3 cp $reference_path reference
done < reference_paths.txt
cat reference/*.fna > reference/reference.fna
reformat.sh in=reference/reference.fna out=reference/reference.ml500.fna ml=500

# run mmseqs2
