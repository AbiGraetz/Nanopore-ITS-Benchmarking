import sys
from Bio import SeqIO
import os

work_dirs = ["sequences/folders/"]

ranks_level = ["p__", "c__", "o__", "f__", "g__"]

for work_dir in work_dirs:
    print("processing " + work_dir)
    fastq_file_path = work_dir + "_train/"
    fastq_files = set()
    for filename in os.listdir(fastq_file_path):
        file_path = os.path.join(fastq_file_path, filename)
        if os.path.isfile(file_path):
            fastq_files.add(filename)
    fastq_files = list(fastq_files)
    fastq_files.sort()



    taxa_dict = {}
    f0 = open("../taxonomy3.txt","r")
    lines0 = f0.readlines()
    f0.close()
    for line in lines0:
        ranks = line[:-1].split(',')
    #     ranks.reverse()
        for i in range(len(ranks)):
            ranks[i] = '_'.join(ranks[i].split(' '))
        taxa_dict[ranks[1]] = "\t".join(ranks[4:-1])
    
    selected_records = []
    for file in fastq_files:
        species_name = file.split('.')[2]
        print("processing " + species_name)
        records = list(SeqIO.parse(fastq_file_path + file, "fastq"))
        for record in records:
            record_id_ori = record.id
            record.id = record_id_ori + "|" + "k__fungi;"
            for i in range(len(taxa_dict[species_name].split('\t'))):
                record.id += ranks_level[i]
                record.id += taxa_dict[species_name].split('\t')[i]
                record.id += ";"
            record.id += "s__"
            record.id += species_name
            record.id += "|"
            record.id += record_id_ori
            
            selected_records.append(record)
    with open(work_dir + "_train_MycoAI/train_sequences.fasta", "w") as output_handle:
        SeqIO.write(selected_records, output_handle, "fasta")
    
    
    

