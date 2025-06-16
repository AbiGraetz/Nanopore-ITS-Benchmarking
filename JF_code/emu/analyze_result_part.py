import os
import pandas as pd
import json
import sys

output_files = []
output_file_path = sys.argv[1]
output_file_name = sys.argv[2]
species_list = []
for filename in os.listdir(output_file_path):
    file_path = os.path.join(output_file_path, filename)
    species_list.append(filename.split('.')[2])
    if os.path.isfile(file_path) and filename[-4:] == ".tsv":
        output_files.append(filename)
output_files.sort()
# print(output_files)


ref_dict = {}
f = open("../taxonomy3.txt", 'r')

lines = f.readlines()
for line in lines:
    splited = line[:-1].split(',')
    ref_dict[splited[1]] = ','.join(splited[2:3] + splited[4:])

    
f.close()

# for i in range(len(output_files)):
results_all = {}

for i in range(len(output_files)):
    output_file = output_files[i]
    correct_species = output_file.split('.')[2]
    result_dict = {}
    
    with open(output_file_path + output_file, 'r') as file:
        for line in file:
            # Split the line into components
            components = line.strip().split('\t')
#             print(components)
            
            if components[2] in ref_dict.keys():
                result_dict[ref_dict[components[2]]] = float(components[-1])
    results_all[correct_species] = result_dict
#             if components[2] in species_dict.keys():
#                 if species_dict[components[2]] == correct_species:
#                     recall_count[correct_species] += 1
#                 precision_count[species_dict[components[2]]] += 1


with open(output_file_name, 'w') as json_file:
    json.dump(results_all, json_file)
print("done")

