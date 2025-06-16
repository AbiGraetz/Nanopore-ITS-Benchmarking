import os
from Bio import SeqIO

taxa_dict = {}
f0 = open("../taxonomy3.txt","r")
lines0 = f0.readlines()
f0.close()
for line in lines0:
    ranks = line[:-1].split(',')
#     ranks.reverse()
    for i in range(len(ranks)):
        ranks[i] = '_'.join(ranks[i].split(' '))
    taxa_dict[ranks[1]] = "\t".join(ranks[2:3] + ranks[4:] + ['N/A', 'N/A'])
print(taxa_dict)

f1 = open("/media/largeData/Jinghang/benchmark/ref.classification", "w")
f1.write("id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain number/type\tsynonym\n")
selected_records = []
records = list(SeqIO.parse("ref.fna", \
                           "fasta"))
for record in records:
#     print("processing " + record.id)
    species_name = record.id
    
    record_id = record.id
    f1.write(record_id)
    f1.write("\t")
    f1.write(taxa_dict[species_name])
    f1.write("\n")
    
f1.close()
