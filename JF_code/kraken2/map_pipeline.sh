#!/bin/bash

# read filtering and mapping for lab2

# make sure there's no column (:) in the reference sequence names as downstream tools don't like this
ref=ref.fna
cpus=8
database_output_path=/customized_db_path
######################
### PIPELINE START ###
######################

source /opt/conda/etc/profile.d/conda.sh
conda activate kraken2_test

cd /seq/folder/
out_dir1=/output/folder/
mkdir -p $out_dir1
for fastq in *.fastq;
do
        fn=${fastq%.fastq}
        echo "Mapping $fastq..."
        kraken2 --db $database_output_path --threads 6 $fastq > $out_dir1/$fn.txt
done



conda deactivate
