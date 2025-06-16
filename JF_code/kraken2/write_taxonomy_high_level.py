from Bio import SeqIO
import gzip
import os
import json


with open("ranks_app.json", 'r') as file:
    ranks = json.load(file)

levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
def calculate_next(rank0, taxo_level):
    ind0 = levels.index(taxo_level)
    next_keys = rank0.keys() & set(levels[ind0+1:])
    return len(next_keys) == 0

accessions_to_rank = {}
f2 = open("taxonomy3_high_level.txt", 'w')
ind = 0
for rank in ranks:
    line = []
    line.append(str(ind))
    line.append(rank['taxid'])
    if 'kingdom' in rank.keys():
        line.append(rank['kingdom'])
    elif not calculate_next(rank, 'kingdom'):
        line.append('NA')
    line.append('subkingdom')
    if 'phylum' in rank.keys():
        line.append(rank['phylum'])
    elif not calculate_next(rank, 'phylum'):
        line.append('NA')
    if 'class' in rank.keys():
        line.append(rank['class'])
    elif not calculate_next(rank, 'class'):
        line.append('NA')
    if 'order' in rank.keys():
        line.append(rank['order'])
    elif not calculate_next(rank, 'order'):
        line.append('NA')
    if 'family' in rank.keys():
        line.append(rank['family'])
    elif not calculate_next(rank, 'family'):
        line.append('NA')
    if 'genus' in rank.keys():
        line.append(rank['genus'])
    elif not calculate_next(rank, 'genus'):
        line.append('NA')
    if 'species' in rank.keys():
        line.append(rank['species'])
    elif not calculate_next(rank, 'species'):
        line.append('NA')
    ind += 1
    f2.write(','.join(line))
    f2.write('\n')

f2.close()

