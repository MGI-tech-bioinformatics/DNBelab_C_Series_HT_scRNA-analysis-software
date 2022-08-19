import argparse,os
from multiprocessing import Pool
from DNBC4tools.tools.utils import logging_call
from DNBC4tools.__init__ import _root_dir

parser = argparse.ArgumentParser(description='QC,star and anno')
parser.add_argument('--name',help='sample name', type=str)
parser.add_argument('--outdir',help='output dir, default is current directory', default=os.getcwd())
parser.add_argument('--cDNAfastq1',help='cDNAR1 fastq file, Multiple files are separated by commas.', required=True)
parser.add_argument('--cDNAfastq2',help='cDNAR2 fastq file, Multiple files are separated by commas.', required=True)
parser.add_argument('--cDNAconfig',
    help='whitelist file in JSON format for cDNA fastq.',
    default='%s/config/DNBelabC4_scRNA_beads_readStructure.json'%_root_dir)
parser.add_argument('--oligofastq1',help='oligoR1 fastq file, Multiple files are separated by commas.',required=True)
parser.add_argument('--oligofastq2',help='oligoR2 fastq file, Multiple files are separated by commas.',required=True)
parser.add_argument('--oligoconfig',
    help='whitelist file in JSON format for oligo fastq.',
    default='%s/config/DNBelabC4_scRNA_oligo_readStructure.json'%_root_dir)
parser.add_argument('--starIndexDir',type=str, help='star index dir.')
parser.add_argument('--gtf',type=str, help='gtf file')
parser.add_argument('--thread',type=int, default=4,help='Analysis threads.')
parser.add_argument('--mixseq',action='store_true',help='If cDNA and oligo sequence in one chip, add this parameter.')
parser.add_argument('--no_introns',action='store_true',help='Not include intronic reads in count.')
args = parser.parse_args()

def cDNA_para():
    cDNA_in1 = open('%s/01.data/cDNAin1'%args.outdir,'w')
    cDNA_in2 = open('%s/01.data/cDNAin2'%args.outdir,'w')
    for list in args.cDNAfastq1.strip().split(','):
        cDNA_in1.write(list+'\n')
    cDNA_in1.close()
    for list in args.cDNAfastq2.strip().split(','):
        cDNA_in2.write(list+'\n')
    cDNA_in2.close()
    cDNA_conf = open('%s/01.data/cDNA_para'%args.outdir,'w')
    cDNA_conf.write('in1=%s/01.data/cDNAin1'%args.outdir+'\n')
    cDNA_conf.write('in2=%s/01.data/cDNAin2'%args.outdir+'\n')
    cDNA_conf.write('config=%s'%args.cDNAconfig+'\n')
    cDNA_conf.write('cbdis=%s/01.data/cDNA_barcode_counts_raw.txt'%args.outdir+'\n')
    cDNA_conf.write('report=%s/01.data/cDNA_sequencing_report.csv'%args.outdir+'\n')
    cDNA_conf.close()
    
def oligo_para():
    oligo_conf = open('%s/01.data/oligo_para'%args.outdir,'w')
    oligo_conf.write('in1=%s'%args.oligofastq1+'\n')
    oligo_conf.write('in2=%s'%args.oligofastq2+'\n')
    if args.mixseq:
        oligo_conf.write('config=%s/config/DNBelabC4_scRNA_oligomix_readStructure.json'%_root_dir+'\n')
    else:
        oligo_conf.write('config=%s'%args.oligoconfig+'\n')
    oligo_conf.write('cbdis=%s/01.data/Index_barcode_counts_raw.txt'%args.outdir+'\n')
    oligo_conf.write('report=%s/01.data/Index_sequencing_report.csv'%args.outdir+'\n')
    oligo_conf.write('outFq=%s/01.data/Index_reads.fq.gz'%args.outdir+'\n')
    oligo_conf.write('threads=%s'%args.thread+'\n')
    oligo_conf.close()

def cDNA_qcstaranno():
    cDNA_para()
    cDNA_star_cmd = '%s/soft/scStar --outSAMattributes singleCell --outSAMtype BAM Unsorted --genomeDir %s --outFileNamePrefix %s/01.data/ --stParaFile %s/01.data/cDNA_para --outSAMmode NoQS --runThreadN %s --limitOutSJcollapsed 10000000 --limitIObufferSize 350000000'\
        %(_root_dir,args.starIndexDir,args.outdir,args.outdir,args.thread)
    cDNA_anno_cmd = ['%s/soft/Anno -I %s/01.data/Aligned.out.bam -a %s -L %s/01.data/cDNA_barcode_counts_raw.txt -o %s/01.data -c %s -m chrM -B %s --anno 1 '\
        %(_root_dir,args.outdir,args.gtf,args.outdir,args.outdir,args.thread,args.cDNAconfig)]
    if args.no_introns:
        pass
    else:
        cDNA_anno_cmd += ['--intron']
    cDNA_anno_cmd = ' '.join(cDNA_anno_cmd)
    return cDNA_star_cmd,cDNA_anno_cmd

def oligo_qc():
    oligo_para()
    oligo_qc_cmd = '%s/soft/parseFq %s/01.data/oligo_para'%(_root_dir,args.outdir)
    return oligo_qc_cmd

if __name__ == '__main__':
    cDNA_star_cmd,cDNA_anno_cmd = cDNA_qcstaranno()
    oligo_qc_cmd = oligo_qc()
    mission = [oligo_qc_cmd,cDNA_star_cmd]
    pool = Pool(2)
    for i in mission:
        pool.apply_async(logging_call,(i,'data',args.outdir,))
    pool.close()
    pool.join()
    logging_call(cDNA_anno_cmd,'data',args.outdir)


