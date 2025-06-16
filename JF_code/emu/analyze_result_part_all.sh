#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate jupyter_gui



python3 analyze_result_part.py /path/to/emu/output/folder/ output_json_file.json

conda deactivate


