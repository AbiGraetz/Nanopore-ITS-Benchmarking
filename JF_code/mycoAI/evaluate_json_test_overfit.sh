#!/bin/bash



source /opt/conda/etc/profile.d/conda.sh
conda activate jupyter_visulize

for i in {0..9}; do
    
    python3 evaluate_json_test_overfit.py /classified/files/folder/ /output_json_file_${i}.json ${i}

    
done

conda deactivate


