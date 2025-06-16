#!/bin/bash



source /opt/conda/etc/profile.d/conda.sh
conda activate jupyter_visulize


python3 analyze_result_part.py /out/folder output_json_file.json


conda deactivate


