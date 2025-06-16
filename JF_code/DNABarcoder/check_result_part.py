import sys
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import math
import json




input_file_name = sys.argv[1]

output_file_name = sys.argv[2]

mock_ind = int(sys.argv[3]) - 1

df = pd.read_csv(input_file_name, delimiter = '\t')

# Output the column with the attribute 'species'
species_column = df['ReferenceID'].to_list()
qid_column = df['ID'].to_list()

print(species_column[:10])
assert len(qid_column) == len(species_column)



df = pd.read_csv('../Mock1Abundance.csv')
species_list = df['organism'].to_list()
sequence_abundances = df['Mock1'].to_list()
assert len(species_list) == len(sequence_abundances)
species_sequence_number_mock1 = {}
for i in range(len(species_list)):
    species = species_list[i]
    species = "_".join(species.split(' '))
    species_sequence_number_mock1[species] = math.ceil(float(sequence_abundances[i]))

    
df = pd.read_csv('../Mock2Abundance.tsv')
species_list = df['organism'].to_list()
sequence_abundances = df['sequence_abundance'].to_list()
assert len(species_list) == len(sequence_abundances)
species_sequence_number_mock2 = {}
for i in range(len(species_list)):
    species = species_list[i]
    species = "_".join(species.split(' '))
    species_sequence_number_mock2[species] = math.ceil(float(sequence_abundances[i]))
    
    
df = pd.read_csv('../Mock3Abundance.csv')
species_list = df['organism'].to_list()
sequence_abundances = df['abundance'].to_list()
assert len(species_list) == len(sequence_abundances)
species_sequence_number_mock3 = {}
for i in range(len(species_list)):
    species = species_list[i]
    species = "_".join(species.split(' '))
    species_sequence_number_mock3[species] = math.ceil(float(sequence_abundances[i]))
    

    
species_sequence_number_mocks_all = [species_sequence_number_mock1, species_sequence_number_mock2, species_sequence_number_mock3]
species_sequence_number_mock = species_sequence_number_mocks_all[mock_ind]


ref_dict = {}
with open("../taxonomy.json", 'r') as json_file:
    ref_dict_base = json.load(json_file)

for key in ref_dict_base.keys():
    ref_dict[key] = ','.join(["Fungi"] + ref_dict_base[key].split('\t') + [key])

results_all = {}
species_list.sort()



total_index = 0
for j in range(len(species_list)):
    result_dict = {}
    key_id_set = set()
    cou = 0
    query_set = set()
    correct = "_".join(species_list[j].split(' '))
    for i in range(species_sequence_number_mock[correct]):
        if not isinstance(species_column[i+total_index], float):
            if qid_column[i+total_index] not in query_set:
                query_set.add(qid_column[i+total_index])
                if species_column[i+total_index] in key_id_set:
                    result_dict[ref_dict[species_column[i+total_index]]] += 1
                else:
                    key_id_set.add(species_column[i+total_index])
                    result_dict[ref_dict[species_column[i+total_index]]] = 1
    total_index += species_sequence_number_mock[correct]
    results_all[correct] = result_dict



with open(output_file_name, 'w') as json_file:
    json.dump(results_all, json_file)
print("done")
