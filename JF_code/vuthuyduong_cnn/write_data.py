import sys
from Bio import SeqIO
import os

work_dir = "/work/dir"


fastq_file_path = "fasta/folder/path"
fastq_files = set()
for filename in os.listdir(fastq_file_path):
    file_path = os.path.join(fastq_file_path, filename)
    if os.path.isfile(file_path):
        fastq_files.add(filename)
fastq_files = list(fastq_files)
fastq_files.sort()
print(fastq_files)



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
print(taxa_dict)



f1 = open("/path/to/train_taxa.txt", "w")




f1.write("#Sequence ID\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
selected_records = []
for file in fastq_files:
    print("processing " + file[:-6])
    species_name = file.split('.')[2]
    records = list(SeqIO.parse(fastq_file_path + file, "fastq"))
    for record in records:
        selected_records.append(record)
        record_id = record.id
        f1.write(record_id)
        f1.write("\t")
        f1.write(taxa_dict[species_name])
        f1.write("\t")
        f1.write(species_name)
        f1.write("\n")
with open("/path/to/train_sequences.fasta", "w") as output_handle:
    SeqIO.write(selected_records, output_handle, "fasta")
f1.close()
    
    
    

