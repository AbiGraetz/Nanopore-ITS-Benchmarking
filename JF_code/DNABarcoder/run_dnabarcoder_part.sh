#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate dnabarcoder_test

# parameter -ml was modified when we doing the experiment on different ref database 
python dnabarcoder.py predict -i ref.fasta -c ref.classification -st 0.7 -et 1 -s 0.001 -rank species -higherrank genus,family,order,class,phylum -ml 800 -sim dnabarcoder/CBSITS.sim
python dnabarcoder.py best -i dnabarcoder/consensus.cutoffs.json -c ref.classification
cd dnabarcoder


python dnabarcoder.py search -i path_to_combined_sequences.fasta -r ref.fasta -ml 800

python dnabarcoder.py classify -i dnabarcoder/path_to_combined_sequences.consensus_BLAST.bestmatch -c ref.classification -cutoff 0.994 -rank species -confidence 0.8334

mkdir -p /output/path/
mv dnabarcoder/path_to_combined_sequences* /output/path/


conda deactivate
