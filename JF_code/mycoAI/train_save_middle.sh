#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate mycoai



python3 train_model_with_mid.py /path/to/train_sequences.fasta --out /path/to/model.pt


conda deactivate

