#!/bin/bash

rank_levels=("species" "genus" "phylum")
Qmin_scores=("min15" "min17")
Mocks=("Mock1" "Mock2" "Mock3")
for rank_level in "${rank_levels[@]}"; do
    for Mock in "${Mocks[@]}"; do
        for Qmin_score in "${Qmin_scores[@]}"; do
            python3 draw_recall_gsd.py common_species.json species_ncbi_dict.json \
                ../taxonomy.json  gsd_path.txt $Mock $Qmin_score $rank_level \
                /output/folder/path/ golden 
        done
    done
done
