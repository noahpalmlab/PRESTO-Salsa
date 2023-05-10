import re
import csv
import collections
import os
import difflib
import gzip

directory = os.getcwd()

all_seq = list()
all_files = os.listdir(directory)
for f in all_files:
    if f.endswith("R1_001.fastq.gz"):
        all_seq.append(f)

reference = csv.reader(open('barcode_dictionary_2022.11.04.csv','r'))
mapping = {}
for row in reference:
    bc, prot = row
    mapping[bc] = prot

final_list = list()
all_bcs = list(set(mapping.keys()))
all_prots = list()
for item in all_bcs:
    all_prots.append(mapping[item])

for seq in all_seq:
    print("seq", seq)
    seq_name = seq.split('_')[0]
    with gzip.open(seq,'rt') as start_file:
    	raw_data = start_file.readlines()
    sequence_data = raw_data[1::4]
    start_file.close()
    barcodes = list()
    for sequence in sequence_data:
        barcode_find = re.search('CTGTACAAGTAA.*GAATTCGATATC', sequence)
        if not barcode_find:
            pass
        else:
            barcode_found = barcode_find.group(0)
            barcode = barcode_found[12:-12]
            barcodes.append(barcode)

    all_counts = list()
    all_counts.append(seq_name)
    bc_count = collections.Counter(barcodes)
    total_count = float(sum(bc_count.values()))
    for item in all_bcs:
        if item in bc_count.keys():
            all_counts.append(bc_count[item])
        else:
            all_counts.append(0)
    all_counts_tuple = tuple(all_counts)
    final_list.append(all_counts_tuple)


all_bcs.insert(0,'Barcode')
all_bcs_tuple = tuple(all_bcs)
all_prots.insert(0,'Protein')
all_proteins = tuple(all_prots)
final_list.insert(0,all_bcs_tuple)
final_list.insert(0,all_proteins)

with open('Mapped_bcs.csv','w', newline='') as f:
    w = csv.writer(f)
    w.writerows(final_list)
