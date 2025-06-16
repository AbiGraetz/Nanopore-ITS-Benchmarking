#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate mycoai


cd /path/to/seq/folder
for fasta in *.fasta;
do
        fn=${fasta%.fasta}
    mkdir -p /output/folder

    mycoai-classify $fasta --model model.pt --out /output/folder/${fn}.csv
    
done


conda deactivate
