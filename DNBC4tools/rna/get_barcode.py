#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os,argparse
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='trans hex to ATCG')
parser.add_argument('-r','--raw', type=str, help='The cDNA_barcode_counts_raw.txt')
parser.add_argument('-c','--select',type=str, help='The cell_calling beads_barcodes_hex.txt')
parser.add_argument('-o','--outdir',type=str, help='Result output dir')
args = parser.parse_args()

def seq_comp(seq):
    nt_comp = {'A':'0', 'C':'1', 'G':'2', 'T':'3'}
    length = len(seq)-1
    sum = 0
    for k,v in enumerate(seq.upper()):
        sum += int(nt_comp[v])*(4**(length-k))
    return str('%010x'%sum).upper()

barcode_all = pd.read_table(args.raw,sep = '\t',header=None)
barcode_all.columns = ['barcode','count']
barcode_all["hex"]= barcode_all["barcode"].map(seq_comp)

select_barcode = []
with open(args.select,'r') as select:
    for line in select:
        line = line.strip()
        select_barcode.append(line)

select_df = barcode_all.loc[barcode_all['hex'].isin(select_barcode)]

barcode_all['barcode'].to_csv(os.path.join(args.outdir,'beads_barcode_all.txt'),index=False,header=False)
select_df['barcode'].to_csv(os.path.join(args.outdir,'beads_barcodes.txt'),index=False,header=False)
