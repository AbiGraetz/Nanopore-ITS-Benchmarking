#!/bin/bash



source /opt/conda/etc/profile.d/conda.sh
conda activate vuthuyduong3

out_dir1=/output/folder
mkdir -p $out_dir1
python3 classifyCNN_test_overfitting.py -i test_sequences.fasta

    
conda deactivate


