import csv
import os
import pandas as pd
import json
import sys


csv_file_path = sys.argv[1]
output_file_name = sys.argv[2]
ind_csv = sys.argv[3]

csv_files = set()
for filename in os.listdir(csv_file_path):
    file_path = os.path.join(csv_file_path, filename)
    if os.path.isfile(file_path) and filename.endswith(str(ind_csv) + ".csv"):
        csv_files.add(filename)
csv_files = list(csv_files)
csv_files.sort()
# print(csv_files)
print(len(csv_files))


ref_dict = {}
with open("../taxonomy.json", 'r') as json_file:
    ref_dict_base = json.load(json_file)

for key in ref_dict_base.keys():
    ref_dict[key] = ','.join(["Fungi"] + ref_dict_base[key].split('\t') + [key])


cou1 = {}
cou2 = {}
all_species = set()
for i in range(len(csv_files)):
    file = csv_files[i]
    cou1[file.split('.')[2]] = 0
    cou2[file.split('.')[2]] = 0
    all_species.add(file.split('.')[2])



results_all = {}
for file in csv_files:
    result_dict = {}
    key_id_set = set()
    species_name = file.split('.')[2]
    lines = open(csv_file_path + file, "r")
    correct_number = 0
    sequencs_ids = set()
    for line in lines:
        sequence_id = line[:-1].split(',')[0]
        predicted_species = line[:-1].split(',')[-1]
        
        if predicted_species in all_species:
            if predicted_species in key_id_set:
                result_dict[ref_dict[predicted_species]] += 1
            else:
                result_dict[ref_dict[predicted_species]] = 1
                key_id_set.add(predicted_species)
            
            sequencs_ids.add(sequence_id)
    results_all[species_name] = result_dict


with open(output_file_name, 'w') as json_file:
    json.dump(results_all, json_file)


