
source /opt/conda/etc/profile.d/conda.sh
conda activate jupyter_visulize
Mock_ind=1
python3 check_result_part.py /path_to_BLAST.species.classified /Mock1_json_output.json ${Mock_ind}


conda deactivate
