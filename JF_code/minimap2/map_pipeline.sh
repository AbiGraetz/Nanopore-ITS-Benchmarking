#!/bin/bash

# read filtering and mapping for lab2

# make sure there's no column (:) in the reference sequence names as downstream tools don't like this
ref=ref.fna
cpus=8

######################
### PIPELINE START ###
######################

source /opt/conda/etc/profile.d/conda.sh
conda activate env

cd /seq/folder/
out_dir1=/output/folder
# mkdir -p $out_dir1
for fastq in *.fastq;
do
        fn=${fastq%.fastq}
        echo "Mapping $fastq..."
        minimap2 -x map-ont --secondary=no $ref $fastq -o $out_dir1/${fn}.paf
done


conda deactivate

