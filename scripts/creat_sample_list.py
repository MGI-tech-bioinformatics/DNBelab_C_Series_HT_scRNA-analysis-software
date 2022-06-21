#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import pandas as pd
import argparse,os

parser = argparse.ArgumentParser(description='Convert tmp.txt to sample.txt')
parser.add_argument('-i','--infile', metavar='FILE', type=str,help='tmp.txt for analysis',default='tmp.txt')
parser.add_argument('-d','--indir',help='storage outfile',default=os.getcwd())
parser.add_argument('--mix',action='store_true',help='Whether to mix and match lane and barcode')
args = parser.parse_args()

infile = args.infile
indir = args.indir

test = pd.read_table(os.path.join(indir,infile),sep='\t',index_col=None)

df_for_analysis = pd.DataFrame(columns=['sampleID', 'cDNA_PATH','oligo_PATH','species'])
df_for_analysis['sampleID'] = test['sample']
df_for_analysis['species'] = test['species']

for i in range(test.shape[0]):
    cDNA = []
    for j in test.iloc[i,1].split(','):
        chip = j.split('/')[-1]
        a = test.iloc[i,2].split(',')
        b = test.iloc[i,3].split(',')
        if not args.mix:
            for k in range(len(test.iloc[i,2].split(','))):
                cDNA.append(j+'/'+ a[k] +'/'+ chip+'_'+a[k] +'_' + b[k] + '_1.fq.gz')
                cDNA.append(j+'/'+ a[k] +'/'+ chip+'_'+a[k] +'_' + b[k] + '_2.fq.gz')
        else:
            for k1 in range(len(test.iloc[i,2].split(','))):
                for k2 in range(len(test.iloc[i,3].split(','))):
                    cDNA.append(j+'/'+ a[k1] +'/'+ chip+'_'+a[k1] +'_' + b[k2] + '_1.fq.gz')
                    cDNA.append(j+'/'+ a[k1] +'/'+ chip+'_'+a[k1] +'_' + b[k2] + '_2.fq.gz')

    Oligo = []
    for j in test.iloc[i,4].split(','):
        chip = j.split('/')[-1]
        a = test.iloc[i,5].split(',')
        b = test.iloc[i,6].split(',')
        if not args.mix:
            for k in range(len(test.iloc[i,5].split(','))):
                Oligo.append(j+'/'+ a[k] +'/'+ chip+'_'+a[k] +'_' + b[k] + '_1.fq.gz')
                Oligo.append(j+'/'+ a[k] +'/'+ chip+'_'+a[k] +'_' + b[k] + '_2.fq.gz')
        else:
            for k1 in range(len(test.iloc[i,5].split(','))):
                for k2 in range(len(test.iloc[i,6].split(','))):
                    Oligo.append(j+'/'+ a[k1] +'/'+ chip+'_'+a[k1] +'_' + b[k2] + '_1.fq.gz')
                    Oligo.append(j+'/'+ a[k1] +'/'+ chip+'_'+a[k1] +'_' + b[k2] + '_2.fq.gz')
    
    cDNA_fq1 = ','.join(cDNA[::2])
    cDNA_fq2 = ','.join(cDNA[1::2])
    df_for_analysis.iloc[i,1] = str(cDNA_fq1)+';'+str(cDNA_fq2)
    Oligo_fq1 = ','.join(Oligo[::2])
    Oligo_fq2 = ','.join(Oligo[1::2])
    df_for_analysis.iloc[i,2] = str(Oligo_fq1)+';'+str(Oligo_fq2)

df_for_analysis.to_csv(os.path.join(indir,'sample.txt'),sep='\t',index=None,header=False)
