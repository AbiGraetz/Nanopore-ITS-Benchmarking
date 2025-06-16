import sys
import os
import json
import pandas as pd


output_files = []
output_file_path = sys.argv[1]
output_file_name = sys.argv[2]
for filename in os.listdir(output_file_path):
    file_path = os.path.join(output_file_path, filename)
    if os.path.isfile(file_path):
        output_files.append(filename)
output_files.sort()
print(output_files)

ref_dict = {}
f = open("../taxonomy3.txt", 'r')

lines = f.readlines()
for line in lines:
    splited = line[:-1].split(',')
    ref_dict[splited[1]] = ','.join(splited[2:3] + splited[4:])

    
f.close()

ref_dict_high_level = {}
f = open("taxonomy3_high_level.txt", 'r')

lines = f.readlines()
for line in lines:
    splited = line[:-1].split(',')
    ref_dict_high_level[splited[1]] = ','.join(splited[2:3] + splited[4:])

    
f.close()

f = open("../tax_id.txt", "r")
lines = f.readlines()
f.close()
species_dict = {}
for i in range(len(lines)):
    spilited = lines[i][:-1].split(' ')
    if len(spilited) == 2:
        species_dict[spilited[1]] = spilited[0]
print(species_dict)
print(len(species_dict))


results_all = {}

for i in range(len(output_files)):
    output_file = output_files[i]
    correct_species = output_file.split('.')[2]
    species_name = correct_species
    result_dict = {}
    key_id_set = set()
    cou = 0
    query_set = set()
    with open(output_file_path + output_file, 'r') as file:
        for line in file:
            # Split the line into components
            components = line.strip().split('\t')
            if components[1] not in query_set:
                query_set.add(components[1])
                if components[2] in species_dict.keys():
                    name = species_dict[components[2]]
                    if name in key_id_set:
                        result_dict[ref_dict[name]] += 1
                    else:
                        result_dict[ref_dict[name]] = 1
                        key_id_set.add(name)
                elif components[2] in ref_dict_high_level.keys():
                    if ref_dict_high_level[components[2]] in key_id_set:
                        result_dict[ref_dict_high_level[components[2]]] += 1
                        key_id_set.add(ref_dict_high_level[components[2]])
                    else:
                        result_dict[ref_dict_high_level[components[2]]] = 1
                else:
                    if not components[2] == '0':
                        print(components[2])
    results_all[species_name] = result_dict



with open(output_file_name, 'w') as json_file:
    json.dump(results_all, json_file)
print("done")



