import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json
import pandas as pd
import os
import sys




common_species_json = sys.argv[1]
ncbi_species_json = sys.argv[2]
taxonomy_json = sys.argv[3]
result_paths_file = sys.argv[4]
dataset_index = [sys.argv[5], sys.argv[6]]
rank_level = sys.argv[7]
out_path = sys.argv[8]
ref_set = sys.argv[9]

f_out_data = open(out_path + rank_level + "_level_" + dataset_index[0] + "_" + dataset_index[1] + "_" + ref_set + "_precision.csv", 'w')

f_out_data2 = open(out_path + rank_level + "_level_" + dataset_index[0] + "_" + dataset_index[1] + "_" + ref_set + "_base_precision.csv", 'w')

f_out_data_rate = open(out_path + rank_level + "_level_" + dataset_index[0] + "_" + dataset_index[1] + "_" + ref_set + "_precision_rate.csv", 'w')


with open(common_species_json) as json0:
    common_species = json.load(json0)

with open(ncbi_species_json) as json0:
    species_ncbi_dict = json.load(json0)


species_to_phylum = {}
species_to_genus = {}
phylum_set = set()
genus_set = set()
with open(taxonomy_json) as json0:
    species_rank = json.load(json0)


with open("Mock1.json") as json0:
    mock1_base = json.load(json0)
with open("Mock2.json") as json0:
    mock2_base = json.load(json0)
with open("Mock3.json") as json0:
    mock3_base = json.load(json0)

common_species = [key for key, value in mock1_base.items()]
common_species.sort()
print(len(common_species))
mock1_base_filtered = {key: value for key, value in mock1_base.items() if key in common_species}
mock2_base_filtered = {key: value for key, value in mock2_base.items() if key in common_species}
mock3_base_filtered = {key: value for key, value in mock3_base.items() if key in common_species}
mock1_base = mock1_base_filtered
mock2_base = mock2_base_filtered
mock3_base = mock3_base_filtered

species_list = list(mock1_base.keys())
species_list.sort()
assert set(mock1_base.keys()) == set(mock2_base.keys())
assert set(mock2_base.keys()) == set(mock3_base.keys())
species_keys = list((mock1_base.keys()))
species_keys.sort()
mock1_base_sequence_number = []
mock2_base_sequence_number = []
mock3_base_sequence_number = []
for i in range(len(species_keys)):
    mock1_base_sequence_number.append(mock1_base[species_keys[i]])
    mock2_base_sequence_number.append(mock2_base[species_keys[i]])
    mock3_base_sequence_number.append(mock3_base[species_keys[i]])

result_paths = []
with open(result_paths_file, 'r') as file:
        for line in file:
            # Strip the trailing newline and append to the list
            result_paths.append(line.rstrip('\n'))


mock_index = int(dataset_index[0][-1:])

mock_bases = [np.array(mock1_base_sequence_number), np.array(mock2_base_sequence_number),\
              np.array(mock3_base_sequence_number)]
mock_base = mock_bases[mock_index - 1]


file_list = []
for result_path in result_paths:
    for root, dirs, files in os.walk(result_path):
        for file in files:
            if not result_path[-8:] == "kraken2/":
                if dataset_index[0] in file and dataset_index[1] in file and file[-5:] == '.json':
                    file_list.append(os.path.join(result_path, file))
            else:
                if dataset_index[0] in file and dataset_index[1] in file and rank_level in file:
                    file_list.append(os.path.join(result_path, file))

print(file_list)

for species1 in common_species:
    genus_rank = species_rank[species1].split('\t')[-1]
    if genus_rank[-1] == '_':
        genus_rank = genus_rank[:-1]
    genus_set.add(genus_rank)
    phylum_set.add(species_rank[species1].split('\t')[0])
    species_to_genus[species1] = genus_rank
    species_to_phylum[species1] = species_rank[species1].split('\t')[0]


species_to_genus["Leptosphaeria_maculans"] = "Leptosphaeria"



genus_set = list(genus_set)
phylum_set = list(phylum_set)
genus_set.sort()
phylum_set.sort()

print(len(genus_set))
print(len(phylum_set))
print(genus_set)
print(phylum_set)

mock_base_genus = []
for i in range(len(genus_set)):
    mock_base_genus.append(0)
assert len(common_species) == len(mock_base)
for i in range(len(common_species)):
    mock_base_genus[genus_set.index(species_to_genus[common_species[i]])] += mock_base[i]

mock_base_phylum = []
for i in range(len(phylum_set)):
    mock_base_phylum.append(0)
assert len(common_species) == len(mock_base)
for i in range(len(common_species)):
    mock_base_phylum[phylum_set.index(species_to_phylum[common_species[i]])] += mock_base[i]



plt.figure(figsize=(18, 6))
methods = []


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





for i in range(len(file_list)):


    mock_base_species = []
    for j in range(len(common_species)):
        mock_base_species.append(0)
    assert len(common_species) == len(mock_base)



    file = file_list[i]
    print(file)
    print(result_paths[i].split('/')[6])
    f_out_data.write("_".join(result_paths[i].split('/')[6:-1]))
    f_out_data.write(",")
    f_out_data2.write("_".join(result_paths[i].split('/')[6:-1]))
    f_out_data2.write(",")
    f_out_data_rate.write("_".join(result_paths[i].split('/')[6:-1]))
    f_out_data_rate.write(",")


    correct_classified = []
    correct_classified_genus = []
    for j in range(len(genus_set)):
        correct_classified_genus.append(0)
    correct_classified_phylum = []
    for j in range(len(phylum_set)):
        correct_classified_phylum.append(0)
    if result_paths[i].split('/')[6] == "emu_mix":
        # print(results.keys())
        for key in common_species:
            correct_count = 0
            for key2 in results.keys():
                if key2.split(',')[-2] == species_to_genus[key]:
                    correct_classified_genus[genus_set.index(species_to_genus[key])] += 0
            correct_classified.append(correct_count)
    
    elif result_paths[i].split('/')[6] == "kraken2":
        print("here")
        df = pd.read_excel(file)
        df_sorted = df.sort_values(by=df.columns[0])
        correct_classified = np.array(list(df_sorted['correctly mapped'])[:54])
        mock_base_species = np.array(list(df_sorted['mapped to this species'])[:54])


    else:
        with open(file, 'r') as json0:
            results = json.load(json0)
        
        results_filtered = {key: value for key, value in results.items() if key in common_species}
        results = results_filtered
        for key in common_species:
            correct_count = 0
            classify_result = results[key]
            for key2 in classify_result:
                if key2.split(',')[-1] == key:
                    correct_count += classify_result[key2]
                    mock_base_species[common_species.index(key)] += classify_result[key2]
                elif key2.split(',')[-1] in species_ncbi_dict.keys() and species_ncbi_dict[key2.split(',')[-1]] == key:
                    correct_count += classify_result[key2]
                    mock_base_species[common_species.index(key)] += classify_result[key2]
                elif key2.split(',')[-1] in common_species:
                    mock_base_species[common_species.index(key2.split(',')[-1])] += classify_result[key2]
                elif key2.split(',')[-1] in species_ncbi_dict.keys():
                    mock_base_species[common_species.index(species_ncbi_dict[key2.split(',')[-1]])] += classify_result[key2]
            correct_classified.append(correct_count)
    print(correct_classified)
    y_axis = []
    
    methods.append(' '.join(result_paths[i].split('/')[6:]))

    arr_correct_classified = np.array(correct_classified)
    arr_mock_base_species  = np.array(mock_base_species)

    result_correct = np.where(
        arr_mock_base_species == 0,
        0, 
        arr_correct_classified / arr_mock_base_species
    )
    f_out_data.write(",".join(map(str, correct_classified)))
    f_out_data.write("\n")
    f_out_data2.write(",".join(map(str, mock_base_species)))
    f_out_data2.write("\n")
    f_out_data_rate.write(",".join(map(str, result_correct.tolist())))
    f_out_data_rate.write("\n")

    plt.scatter(range(i * 65 + 1, i * 65 + len(common_species) + 1), result_correct, \
             color=colors[i], s = 15)
             
    plt.scatter(np.mean(np.array(range(i * 65 + 1, i * 65 + len(common_species) + 1))), np.mean(result_correct), marker='s', \
     s=30, edgecolors='black', linewidths=1, color=colors[-1])
    # Plot the data for the first subplot (Order)
    
plt.title(rank_level + ' level on ' + dataset_index[0] + ' dataset with Q' + dataset_index[1] + ' mapped to ' + ref_set, fontsize=18)
plt.ylabel('Precision Rate', fontsize=16)
plt.ylim(0, 1)
plt.xlim(0, 520)

print(methods)
plt.xticks(np.arange(1, len(file_list) + 1) * (len(common_species) + 12) - len(common_species)//2 - 20, methods, fontsize=16)
plt.legend(loc='best')
# Display the plot
# plt.show()
save_filepath = out_path + rank_level + "_level_" + dataset_index[0] + "_" + dataset_index[1] + "_" + ref_set + "_precision.png"
plt.savefig(save_filepath)

f_out_data.close()
f_out_data2.close()
f_out_data_rate.close()

