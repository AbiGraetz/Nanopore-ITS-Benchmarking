#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate env


python3 check_result_part.py 'output_json_file.json' '/minimap2_output_path'
