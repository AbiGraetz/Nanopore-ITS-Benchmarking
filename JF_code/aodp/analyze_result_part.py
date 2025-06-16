import sys 
import pandas as pd
import os
import json

ref_dict = {}
with open("../taxonomy.json", 'r') as json_file:
    ref_dict_base = json.load(json_file)

for key in ref_dict_base.keys():
    ref_dict[key] = ','.join(["Fungi"] + ref_dict_base[key].split('\t') + [key])

output_files = []
output_file_path = sys.argv[1]
output_file_name = sys.argv[2]
for filename in os.listdir(output_file_path):
    file_path = os.path.join(output_file_path, filename)
    if os.path.isfile(file_path):
        output_files.append(filename)
output_files.sort()
print(output_files)


results_all = {}
for file in output_files:
    species_name = file.split('.')[2]
    result_dict = {}
    key_id_set = set()
    lines = open(output_file_path + file, "r")
    correct_number = 0
    sequencs_ids = set()
    for line in lines:
        sequence_id = line.split('\t')[0]
        predicted_species = line.split('\t')[1]
        if sequence_id not in sequencs_ids and not predicted_species == '-':
            
            if predicted_species in key_id_set:
                result_dict[ref_dict[predicted_species]] += 1
            else:
                result_dict[ref_dict[predicted_species]] = 1
                key_id_set.add(predicted_species)
            sequencs_ids.add(sequence_id)
    results_all[species_name] = result_dict

with open(output_file_name, 'w') as json_file:
    json.dump(results_all, json_file)




