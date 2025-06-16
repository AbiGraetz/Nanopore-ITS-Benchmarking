#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate mycoai

mycoai-train /path/to/train_sequences.fasta --out /path/to/output/model.pt




conda deactivate
