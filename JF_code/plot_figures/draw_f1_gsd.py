import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json
import pandas as pd
import os
import sys
import csv

out_path = sys.argv[1]

files_prefix = []

for root, dirs, files in os.walk(out_path):
    for file in files:
        if file.endswith("_golden_precision.png"):
            files_prefix.append(file[:-len("_precision.png")])

files_prefix.sort()



print(files_prefix)

colors = [
    "#7570b3",  
    "#1b9e77",  
    "#d95f02", 
    "#e7298a",  
    "#66a61e",  
    "#e6ab02",  
    "#a6761d",  
    "#666666",  
    "#2c7fb8",  
    "#e31a1c"   
]

for file_prefix in files_prefix:
    recall_rate_file = file_prefix + "_recall_rate.csv"
    precision_rate_file = file_prefix + "_precision_rate.csv"
    f_out_data = open(out_path + file_prefix + "_f1_score.csv", 'w')
    dataset_index = [file_prefix.split('_')[2], file_prefix.split('_')[3]]
    ref_set = "golden"


    plt.figure(figsize=(18, 6))

    methods = []
    i = 0
    with open(out_path + recall_rate_file, 'r', encoding='utf-8-sig') as fr, open(out_path + precision_rate_file, 'r', encoding='utf-8-sig') as fp:
        reader_r = csv.reader(fr)
        reader_p = csv.reader(fp)
        for row1, row2 in zip(reader_r, reader_p):
            assert row1[0] == row2[0]
            assert len(row1) == len(row2)
            methods.append(row1[0])
            f_out_data.write(row1[0])
            f_out_data.write(",")
            recall_rates = np.array(row1[1:], dtype=float)
            precision_rates = np.array(row2[1:], dtype=float)
            # f1_scores = 2 * precision_rates * recall_rates / (precision_rates + recall_rates)
            denominator = precision_rates + recall_rates
            f1_scores = np.where(
                denominator == 0, 
                0, 
                2 * precision_rates * recall_rates / denominator
            )

            # print(denominator)
            # print(f1_scores)
            
            
            f_out_data.write(",".join(map(str, f1_scores)))
            f_out_data.write("\n")

            plt.scatter(range(i * (len(f1_scores) + 11) + 1, i * (len(f1_scores) + 11) + len(f1_scores) + 1), f1_scores, \
                color=colors[i], s = 15)
            plt.scatter(np.mean(np.array(range(i * (len(f1_scores) + 11) + 1, i * (len(f1_scores) + 11) + len(row1[1:]) + 1))), np.mean(f1_scores), marker='s', \
                s=30, edgecolors='black', linewidths=1, color=colors[-1])
            
            i += 1

        plt.title('Species level on ' + dataset_index[0] + ' dataset with Q' + dataset_index[1] + ' mapped to ' + ref_set, fontsize=18)
        plt.ylabel('F1 score', fontsize=16)
        plt.ylim(0, 1)
        plt.xlim(0, 520)

    plt.xticks(np.arange(1, i+1) * (len(row1[1:]) + 12) - len(row1[1:])//2 - 20, methods, fontsize=16)
    plt.legend(loc='best')

    save_filepath = out_path + "species_level_" + dataset_index[0] + "_" + dataset_index[1] + "_" + ref_set + "_f1.png"
    plt.savefig(save_filepath)


    f_out_data.close()