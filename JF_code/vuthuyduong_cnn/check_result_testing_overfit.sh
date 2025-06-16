#!/bin/bash



source /opt/conda/etc/profile.d/conda.sh
conda activate jupyter_visulize

for i in {0..19}; do
    python3 check_result_part_test_overfit.py Mock_${i}.classified /output_json_${i}.json 

done



conda deactivate


