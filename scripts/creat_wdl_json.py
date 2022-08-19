#!/usr/bin/env python
# -*- encoding: utf-8 -*-

'''
Generate "config.json" and "run.sh" for each sample from the table corresponding to the sample.txt
The format of "sample.txt" is as follows, no header required:
----------------------------------------------------------
test1 cDNA1_L01_1.fq.gz;cDNA1_L01_2.fq.gz    oligo1_L01_1.fq.gz,oligo1_L02_1.fq.gz;oligo1_L01_2.fq.gz,oligo1_L02_2.fq.gz Mouse
test2 cDNA2_L01_1.fq.gz,cDNA2_L02_1.fq.gz;cDNA1_L01_2.fq.gz,cDNA2_L02_2.fq.gz   oligo2_L01_1.fq.gz;oligo2_L01_2.fq.gz  Mouse
test3 cDNA3_L01_1.fq.gz;cDNA3_L01_2.fq.gz    oligo3_L01_1.fq.gz,oligo3_L02_1.fq.gz;oligo3_L01_2.fq.gz,oligo3_L02_2.fq.gz Mouse
-----------------------------------------------------------
The list file includes four columns, and the columns are separated by a delimiter.
The first column is the sample name, the second column is the fastq file of cDNA, the third column is the fastq file of oligo, and the fourth column is the species name of the sample.
A semicolon is used to separate fq1 and fq2. Multiple fastqs are separated by commas. The order of fq1 and fq2 needs to be in one-to-one correspondence.

"wdl.json" needs to be in the same directory as "creat_wdl_json.py". 
Before using it for the first time, please modify "wdl.json" according to your actual value. 
Species information can be added to the DB in json format, which needs to correspond to the "species" column name of "sample.txt".
'''

from collections import defaultdict
import json,os,argparse

parser = argparse.ArgumentParser(description='create run for all samples')
parser.add_argument('-i','--infile', metavar='FILE', type=str,help='sample.txt for analysis, default: sample.txt',default='sample.txt')
parser.add_argument('-m','--mixseq', action='store_true', help='cDNA and oligo sequeencing in one chip')
parser.add_argument('-o','--outdir',help='storage outfile',default=os.getcwd())
args = parser.parse_args()

file_dir = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(file_dir,'wdl.json'),'r') as f:
        configlist = json.load(f)

with open(args.infile, 'r') as f:
    for line in f:
        line = line.strip()
        if len(line.split('\t')) == 4:
            adict = defaultdict()
            lst = line.split('\t')
            os.makedirs("%s/%s"%(args.outdir,lst[0]), exist_ok=True)
            adict["main.Outdir"] = "%s/%s"%(args.outdir,lst[0])
            adict["main.SampleName"] = "%s"%lst[0]
            adict["main.cDNA_Fastq1"] = "%s"%lst[1].split(';')[0]
            adict["main.cDNA_Fastq2"] = "%s"%lst[1].split(';')[1]
            adict["main.Oligo_Fastq1"] = "%s"%lst[2].split(';')[0]
            adict["main.Oligo_Fastq2"] = "%s"%lst[2].split(';')[1]
            adict["main.BeadsBarcode"] = "%s"%configlist['BeadsBarcode']
            if args.mixseq:
                adict["main.OligoBarcode"] = "%s"%configlist['OligoBarcodeMix']
            else:
                adict["main.OligoBarcode"] = "%s"%configlist['OligoBarcode']
            adict["main.Root"] = "%s"%configlist['Root']
            adict["main.Refdir"] = "%s"%configlist['DB'][lst[3]]['refdir']
            adict["main.Gtf"] = "%s"%configlist['DB'][lst[3]]['gtf']
            adict["main.Oligo_type8"] = "%s"%configlist['Oligo_type8']
            adict["main.Species"] = "%s"%lst[3]
            adict["main.expectCellNum"] = configlist['expectCellNum']
            adict["main.calling_method"] = "%s"%configlist['calling_method']
            adict["main.forceCellNum"] = configlist['forceCellNum']
            adict["main.Intron"] = configlist['Intron']
            json_str = json.dumps(adict, indent=4)
            with open('%s/%s/config.json'%(args.outdir,lst[0]), 'w') as json_file:
                json_file.write(json_str)
                json_file.write('\n')

            with open('%s/%s/work.sh'%(args.outdir,lst[0]), 'w') as runfile:
                runfile.write('export PATH=%s/bin:$PATH'%configlist['Env']+'\n')
                runfile.write('export LD_LIBRARY_PATH=%s/lib:$LD_LIBRARY_PATH'%configlist['Env']+'\n')
                runfile.write('java -jar %s/wdl/cromwell-35.jar run -i config.json %s/wdl/DNBC4_scRNA.wdl'%(configlist['Root'],configlist['Root'])+'\n')