import sys
import numpy as np
import json
import pandas as pd
import math

classified_file = sys.argv[1]
output_file_name = sys.argv[2]


labels = ["Candida_metapsilosis", "Candida_parapsilosis", "Pichia_membranifaciens"]
print(labels)
print(len(labels))

ref_dict = {}
with open("../taxonomy.json", 'r') as json_file:
    ref_dict_base = json.load(json_file)

for key in ref_dict_base.keys():
    ref_dict[key] = ','.join(["Fungi"] + ref_dict_base[key].split('\t') + [key])

f = open(classified_file, "r")
lines = f.readlines()
f.close()
print(len(lines))
all_predicted = []
for i in range(1, len(lines)):
    all_predicted.append(lines[i].split('\t')[2])




df = pd.read_csv('../Mock1Abundance_test_overfit.csv')
species_list = df['organism'].to_list()
sequence_abundances = df['Mock1'].to_list()
assert len(species_list) == len(sequence_abundances)
species_sequence_number_mock1 = {}
for i in range(len(species_list)):
    species = species_list[i]
    species = "_".join(species.split(' '))
    species_sequence_number_mock1[species] = sequence_abundances[i]

    
    
species_sequence_number_mock = species_sequence_number_mock1


results_all = {}
total_index = 0
for j in range(len(labels)):
    result_dict = {}
    key_id_set = set()
    for i in range(species_sequence_number_mock[labels[j]]):
        

        name = all_predicted[i+total_index]
        if name in key_id_set:
            result_dict[ref_dict[name]] += 1
        else:
            result_dict[ref_dict[name]] = 1
            key_id_set.add(name)
    total_index += species_sequence_number_mock[labels[j]]
    results_all[labels[j]] = result_dict



with open(output_file_name, 'w') as json_file:
    json.dump(results_all, json_file)

