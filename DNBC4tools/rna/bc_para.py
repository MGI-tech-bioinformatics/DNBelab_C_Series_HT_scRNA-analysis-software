import argparse,os,sys

parser = argparse.ArgumentParser(description='process para')
parser.add_argument('--outdir',help='output dir, default is current directory', default=os.getcwd())
parser.add_argument('--cDNAfastq1',help='cDNAR1 fastq file, Multiple files are separated by commas.', required=True)
parser.add_argument('--cDNAfastq2',help='cDNAR2 fastq file, Multiple files are separated by commas.', required=True)
parser.add_argument('--cDNAconfig',help='whitelist file in JSON format for cDNA fastq.')
parser.add_argument('--oligofastq1',help='oligoR1 fastq file, Multiple files are separated by commas.',required=True)
parser.add_argument('--oligofastq2',help='oligoR2 fastq file, Multiple files are separated by commas.',required=True)
parser.add_argument('--oligoconfig', help='whitelist file in JSON format for oligo fastq.')
parser.add_argument('--genomeDir', help='star index dir.')
parser.add_argument('--gtf', help='gtf file.')
parser.add_argument('--adapter', help='adapter file.')
parser.add_argument('--oligotype', help='oligotype file.')
parser.add_argument('--thread',type=int, default=4,help='Analysis threads.')
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
    cDNA_conf.write('adapter=%s'%args.adapter+'\n')
    cDNA_conf.write('cbdis=%s/01.data/cDNA_barcode_counts_raw.txt'%args.outdir+'\n')
    cDNA_conf.write('report=%s/01.data/cDNA_sequencing_report.csv'%args.outdir+'\n')
    cDNA_conf.close()
    
def oligo_para():
    oligo_conf = open('%s/01.data/oligo_para'%args.outdir,'w')
    oligo_conf.write('in1=%s'%args.oligofastq1+'\n')
    oligo_conf.write('in2=%s'%args.oligofastq2+'\n')
    oligo_conf.write('config=%s'%args.oligoconfig+'\n')
    oligo_conf.write('cbdis=%s/01.data/Index_barcode_counts_raw.txt'%args.outdir+'\n')
    oligo_conf.write('report=%s/01.data/Index_sequencing_report.csv'%args.outdir+'\n')
    oligo_conf.write('outFq=%s/01.data/Index_reads.fq.gz'%args.outdir+'\n')
    oligo_conf.write('threads=%s'%args.thread+'\n')
    oligo_conf.close()

def judgeFilexits(*args):
    for input_files in args:
        for input_file in input_files.split(','):
            if not os.path.exists(input_file): 
                print(" ------------------------------------------------") 
                print("Error: Cannot find input file or dir %s."%(str(input_file))) 
                print(" ------------------------------------------------") 
                sys.exit()
            else:
                pass

if __name__=='__main__':
    judgeFilexits(args.cDNAfastq1,args.cDNAfastq2,args.oligofastq1,args.oligofastq2,\
        args.cDNAconfig,args.oligoconfig,args.genomeDir,args.gtf,args.oligotype,args.adapter)
    cDNA_para()
    oligo_para()
