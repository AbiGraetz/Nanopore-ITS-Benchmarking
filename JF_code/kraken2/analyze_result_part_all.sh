#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate jupyter_gui



python3 analyze_result_part.py /kraken2/output/folder/ output_json_file.json


conda deactivate


