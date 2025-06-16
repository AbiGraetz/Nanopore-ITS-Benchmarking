import pandas as pd
import json
import csv

tax_id_dict = {}

f = open("tax_id.txt", 'r')
lines = f.readlines()
f.close()

f1 = open("seq2tax.map", 'w')
for line in lines:
    line = line[:-1]
    line_spt = line.split(' ')
    
    tax_id_dict[line_spt[0]] = line_spt[1]
    f1.write(line_spt[0])
    f1.write('\t')
    f1.write(line_spt[1])
    f1.write('\n')
    
f1.close()
print(tax_id_dict)


data = []

f = open("taxonomy3.txt", 'r')
lines = f.readlines()
f.close()


for line in lines:
    line_spt = line.split(',')
    species_name = '_'.join(line_spt[1].split(' '))
    data0 = {}
    if species_name in tax_id_dict.keys():
        data0['tax_id']=tax_id_dict[species_name]
        data0['species']=species_name
        data0['genus']=line_spt[8]
        data0['family']=line_spt[7]
        data0['order']=line_spt[6]
        data0['class']=line_spt[5]
        data0['phylum']=line_spt[4]
        data0['kingdom']=line_spt[2]
    data.append(data0)


fieldnames = ['tax_id', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

# Open a file to write to, newline='' is recommended in Python 3
with open('taxonomy.tsv', 'w', newline='') as out_file:
    # Create a DictWriter object with tab as delimiter
    writer = csv.DictWriter(out_file, fieldnames=fieldnames, delimiter='\t')

    # Write the header (optional, but typically very useful)
    writer.writeheader()

    # Write data rows
    for row in data:
        writer.writerow(row)
