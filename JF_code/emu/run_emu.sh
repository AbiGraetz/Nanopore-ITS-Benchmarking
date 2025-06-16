#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate emu_test

cd /seq/folder
for file in *.fasta; do
emu abundance "$file" --db path/to/my_db --output-dir /output/folder --keep-files --keep-counts --output-unclassified 
done


conda deactivate
