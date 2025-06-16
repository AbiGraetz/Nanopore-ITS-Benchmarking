import os
import json
import pandas as pd
from Bio import Entrez
from Bio import Phylo

output_files = []
output_file_paths = ["kraken2_output_folder"]
for output_file_path in output_file_paths:
    for filename in os.listdir(output_file_path):
        file_path = os.path.join(output_file_path, filename)
        if os.path.isfile(file_path):
            output_files.append(file_path)
output_files.sort()
print(output_files[:5])
print(len(output_files))

f = open("../tax_id.txt", "r")
lines = f.readlines()
f.close()
species_dict = {}
for i in range(len(lines)):
    spilited = lines[i][:-1].split(' ')
    if len(spilited) == 2:
        species_dict[spilited[1]] = spilited[0]
# print(species_dict)
print(len(species_dict))

keys_need_to_find = set()
for i in range(len(output_files)):
    output_file = output_files[i]
    correct_species = output_file.split('.')[2]
    species_name = correct_species
    result_dict = {}
    key_id_set = set()
    cou = 0
    query_set = set()
    with open(output_file, 'r') as file:
        for line in file:
            # Split the line into components
            components = line.strip().split('\t')
            if components[2] not in species_dict.keys():
                keys_need_to_find.add(components[2])
                

Entrez.email = "your.email@example.com"

def get_taxonomy_rank(taxid):
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    
    taxonomy_info = records[0]
    scientific_name = taxonomy_info["ScientificName"]
    rank = taxonomy_info["Rank"]
    lineage = taxonomy_info["LineageEx"]
    
    return scientific_name, rank, lineage

ranks = []
not_find_id = []
tax_id_to_scientificname = {}
for key in keys_need_to_find:
    species_rank = {}
    print("processing " + key)
    taxid = key
    species_rank['taxid'] = taxid
    try:
        scientific_name, rank, lineage = get_taxonomy_rank(taxid)
        tax_id_to_scientificname[key] = scientific_name
        for taxon in lineage:
            if taxon['Rank'] == 'kingdom':
                species_rank['kingdom'] = taxon['ScientificName']
            if taxon['Rank'] == 'phylum':
                species_rank['phylum'] = taxon['ScientificName']
            if taxon['Rank'] == 'class':
                species_rank['class'] = taxon['ScientificName']
            if taxon['Rank'] == 'order':
                species_rank['order'] = taxon['ScientificName']
            if taxon['Rank'] == 'family':
                species_rank['family'] = taxon['ScientificName']
            if taxon['Rank'] == 'genus':
                species_rank['genus'] = taxon['ScientificName']
        species_rank[rank] = scientific_name
        ranks.append(species_rank)
    except:
#         scientific_name, rank, lineage = get_taxonomy_rank(taxid)
        not_find_id.append(taxid)


with open('ranks_app.json', 'w') as json_file:
    json.dump(tax_id_to_scientificname, json_file)
