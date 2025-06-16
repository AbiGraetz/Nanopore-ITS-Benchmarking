#!/bin/bash



source /opt/conda/etc/profile.d/conda.sh
conda activate vuthuyduong_cnn_env

out_dir1=/output/folder
mkdir -p $out_dir1
python3 classifyCNN.py -i Mock1_test_sequences.fasta -o /path/to/Mock1.classified


conda deactivate


