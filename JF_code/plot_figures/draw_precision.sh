#!/bin/bash

Qmin_scores=("min15" "min17")
Mocks=("Mock1" "Mock2" "Mock3")
for Mock in "${Mocks[@]}"; do
    for Qmin_score in "${Qmin_scores[@]}"; do
        python3 draw_precision_gsd.py common_species.json species_ncbi_dict.json \
            ../taxonomy.json gsd_path.txt $Mock $Qmin_score species \
            /output/folder/ golden 
    done
done
