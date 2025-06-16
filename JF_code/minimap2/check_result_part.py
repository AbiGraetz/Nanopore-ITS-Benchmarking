import sys
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import json

def pull_mapping_results_v3(fn):
    min_header = ['qseqid', 'qlen', 'qstart', 'qstop', 'strand', 'tname', 'tlen', 'tstart', \
                  'tend', 'nmatch', 'alen', 'mquality']
    tmp_df = pd.read_csv(fn, sep='\t', header = None, usecols=[x for x in range(0,12)], \
                         names=min_header)
    tmp_df['cscore'] = tmp_df['alen']/(tmp_df['alen']-tmp_df['nmatch'])
    sub_df = tmp_df[tmp_df['cscore'] == tmp_df.groupby('qseqid')['cscore'].transform(max)]\
    .reset_index(drop=True)
    hit_df = pd.DataFrame(sub_df.groupby('tname')['cscore'].count().tolist(), \
                          sub_df.groupby('tname')['cscore'].count().index, columns=['count'])
    hit_df.sort_values(by='count', ascending=False, inplace=True)
    for key in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
        hit_df[key] = False
        tmp_df[key] = False
    return hit_df, sub_df


ref_dict = {}
f = open("../taxonomy3.txt", 'r')

lines = f.readlines()
for line in lines:
    splited = line[:-1].split(',')
    ref_dict[splited[1]] = ','.join(splited[2:3] + splited[4:])

    
f.close()


paf_files = []
output_file_name = sys.argv[1]
paf_file_path =  sys.argv[2]
for filename in os.listdir(paf_file_path):
    file_path = os.path.join(paf_file_path, filename)
    if os.path.isfile(file_path):
        paf_files.append(filename)
paf_files.sort()
print(paf_files)


cou1 = {}
cou2 = {}
for i in range(len(paf_files)):
    file = paf_files[i]
    cou1[file.split('.')[2]] = 0
    cou2[file.split('.')[2]] = 0






results_all = {}
for i in range(len(paf_files)):
    file = paf_files[i]
    result_dict = {}
    key_id_set = set()
    species_name = file.split('.')[2]
    
    key0 = file.split('.')[2]
    mapping_results , full_results_df = pull_mapping_results_v3(paf_file_path + file)
    
    names = full_results_df['tname']
    query_seq_ids = full_results_df['qseqid']
    query_set = set()
    for i in range(len(names)):
#         if int(name[-1]) > 2:
#             print(file[:-4])
        name = names[i]
    
        if query_seq_ids[i] not in query_set:
            if name in key_id_set:
                result_dict[ref_dict[name]] += 1
            else:
                result_dict[ref_dict[name]] = 1
                key_id_set.add(name)
            query_set.add(query_seq_ids[i])
    results_all[species_name] = result_dict



with open(output_file_name, 'w') as json_file:
    json.dump(results_all, json_file)





