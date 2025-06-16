#!/bin/bash

# read filtering and mapping for lab2

# make sure there's no column (:) in the reference sequence names as downstream tools don't like this

cpus=2

ref=ref.fasta


######################
### PIPELINE START ###
######################

source /opt/conda/etc/profile.d/conda.sh
conda activate aodp_support


cd /seq/folder
out_path=/output/folder/
mkdir -p $out_path
for fasta in *.fasta;
do
        fn=${fasta%.fasta}
        echo "Mapping $fasta..."
        aodp --match=${fn}.fasta --match-output=$out_path/${fn}.match --oligo-size=20 --reverse-complement=yes $ref
done

conda deactivate
