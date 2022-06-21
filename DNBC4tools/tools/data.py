import os
from .utils import rm_temp, str_mkdir,start_print_cmd,judgeFilexits,change_path,python_path,logging_call
from DNBC4tools.__init__ import _root_dir

class Data:
    def __init__(self, args):
        self.cDNAr1 = args.cDNAfastq1
        self.cDNAr2 = args.cDNAfastq2
        self.oligor1 = args.oligofastq1
        self.oligor2 = args.oligofastq2
        self.thread = args.thread
        self.name = args.name
        self.cDNAconfig = args.cDNAconfig
        self.oligoconfig = args.oligoconfig
        self.outdir = os.path.join(args.outdir,args.name)
        self.starindex = args.starIndexDir
        self.gtf = args.gtf
        self.no_introns = args.no_introns
        self.mixseq = args.mixseq

    def run(self):
        judgeFilexits(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.cDNAconfig,self.oligoconfig,self.starindex,self.gtf)
        str_mkdir('%s/01.data'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        change_path()
        new_python = python_path()

        scRNA_parse_cmd = ['%s %s/rna/star_anno.py --name %s --cDNAfastq1 %s --cDNAfastq2 %s --oligofastq1 %s --oligofastq2 %s --thread %s --cDNAconfig %s --oligoconfig %s --outdir %s --starIndexDir %s --gtf %s'\
            %(new_python,_root_dir,self.name,self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.thread,self.cDNAconfig,self.oligoconfig,self.outdir,self.starindex,self.gtf)]
        if self.no_introns:
            scRNA_parse_cmd += ['--no_introns']
        if self.mixseq:
            scRNA_parse_cmd += ['--mixseq']
        scRNA_parse_cmd = ' '.join(scRNA_parse_cmd)
        start_print_cmd(scRNA_parse_cmd)
        str_mkdir('%s/01.data/raw_matrix'%self.outdir)
        raw_matrix_cmd = '%s/soft/PISA count -@ %s -cb CB -anno-tag GN -umi UB -outdir %s/01.data/raw_matrix %s/01.data/final.bam'\
            %(_root_dir,self.thread,self.outdir,self.outdir)
        logging_call(raw_matrix_cmd,'data',self.outdir)
        rm_temp('%s/01.data/Aligned.out.bam'%self.outdir)

def data(args):
    Data(args).run()

def parse_data(parser):
    parser.add_argument('--name',metavar='NAME',help='sample name.', type=str)
    parser.add_argument('--outdir',metavar='DIR',help='output dir, [default is current directory].', default=os.getcwd())
    parser.add_argument('--cDNAfastq1',metavar='FASTQ',help='cDNAR1 fastq file, Multiple files are separated by commas.', required=True)
    parser.add_argument('--cDNAfastq2',metavar='FASTQ',help='cDNAR2 fastq file, Multiple files are separated by commas, the files order needs to be consistent with cDNAfastq1.', required=True)
    parser.add_argument('--cDNAconfig',metavar='JSON',help='whitelist of cell barcode and structure file in JSON format for cDNA fastq.',default='%s/config/DNBelabC4_scRNA_beads_readStructure.json'%_root_dir)
    parser.add_argument('--oligofastq1',metavar='FASTQ',help='oligoR1 fastq file, Multiple files are separated by commas.',required=True)
    parser.add_argument('--oligofastq2',metavar='FASTQ',help='oligoR2 fastq file, Multiple files are separated by commas, the files order needs to be consistent with oligofastq1.',required=True)
    parser.add_argument('--oligoconfig',metavar='JSON',help='whitelist of cell barcode and structure file in JSON format for oligo fastq.',default='%s/config/DNBelabC4_scRNA_oligo_readStructure.json'%_root_dir)
    parser.add_argument('--thread',metavar='INT',type=int, default=4,help='Analysis threads.')
    parser.add_argument('--starIndexDir',metavar='DIR',type=str, help='star index dir.')
    parser.add_argument('--gtf',metavar='GTF',type=str, help='gtf file.')
    parser.add_argument('--mixseq',action='store_true',help='If cDNA and oligo sequence in one chip, add this parameter.')
    parser.add_argument('--no_introns', action='store_true',help='Not include intronic reads in count.')
    return parser
