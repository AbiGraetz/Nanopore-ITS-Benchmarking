#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate mycoai


for i in {0..9}; do

    cd /path/to/sequences/files
    for fasta in *.fasta;
    do
        fn=${fasta%.fasta}
        mycoai-classify $fasta --model /path/to/model_${i}.pt --out /output/path/${fn}_${i}.csv
    done


done


conda deactivate
