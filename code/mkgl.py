#!/usr/bin/env python

import sys
import pandas as pd
import csv
import scipy.io

ensembls = open("/home/zz2565/scHPF/resources/gencode.v31.annotation.gene_l1l2.pc_TRC_IGC.stripped.txt")
ensembls_dict = {}
for line in ensembls:
    ensembl, name = line.strip().split("\t")
    ensembls_dict[name] = ensembl

sc_data = pd.read_csv(sys.argv[1], index_col = 'Unnamed: 0')
data_names = sc_data.columns
drop_names = []
# delete genename without ensembl id
for name in data_names:
    if name not in ensembls_dict:
        drop_names.append(name)
sc_data = sc_data.drop(columns=drop_names)
data_names = sc_data.columns
ensembl_names = []

# get gene list
for name in data_names:
    ensembl_names.append(ensembls_dict[name])

# write a new csv for schpf 
sc_data = sc_data.transpose()
with open("test_cells", 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(sc_data.columns)
if len(ensembl_names) != len(data_names):
    print(len(ensembl_names), len(data_names))
sc_data.insert(0, "ensemble_id", ensembl_names, True)
sc_data.insert(1, "gene_name", data_names, True)
print(sc_data.head())
sc_data.to_csv(sys.argv[1] + ".umi", sep=' ', header=False, index=False)

# write into a sep
# with open(sys.argv[0]+".genes", 'w') as f:
#    writer = csv.writer(f, delimiter='\t')
#    writer.writerows(zip(ensembl_names, data_names))

