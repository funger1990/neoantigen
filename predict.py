#!/usr/bin/env python

import sys
import re
import os
import gzip
import pandas as pd
import neoantigen

file_hla = 'hla'
file_var = 'variant'
outfile = 'neo'

thread = 1


# --------------------------------------------------------------------
# HLA
with open(file_hla) as f:
    hla_list = f.read().splitlines()

# variant
df_var = pd.read_csv(file_var, sep='\t')

# mutant protein
df_var = neoantigen.mutation_to_protein(df_var)

# chop kmer
df_pep = neoantigen.chop_kmer(df_var[['uid', 'old', 'new']])
pep_list = df_pep['pep_new'].tolist()

# run netmhc
df_net = neoantigen.predict_netmhc(pep_list, hla_list, thread)
df_net.rename(columns={'pep': 'pep_new'}, inplace=True)

# merge
df = pd.merge(df_var, df_pep, how='outer')
df = pd.merge(df, df_net, how='left')

# choose the smallest IC50 for each (uid, hla)
df.sort_values(['hla', 'ic50'], inplace=True)
df.drop_duplicates(['uid', 'hla'], inplace=True)

df = df.loc[df['ic50'] < 500, :]

# output
df.to_csv(outfile, sep='\t', index=False)

