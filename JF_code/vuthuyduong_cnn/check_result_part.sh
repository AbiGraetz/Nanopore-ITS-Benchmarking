#!/bin/bash



source /opt/conda/etc/profile.d/conda.sh
conda activate jupyter_visulize
Mock_ind=1
python3 check_result_part.py /path/to/Mock1.classified output_json_file.json ${Mock_ind}

conda deactivate


