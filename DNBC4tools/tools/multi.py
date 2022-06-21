import os
from DNBC4tools.__init__ import _root_dir

class Multi_list:
    def __init__(self, args):
        self.list = args.list
        self.starIndexDir = args.starIndexDir
        self.gtf = args.gtf
        self.outdir = args.outdir
        self.thread = args.thread
        self.cDNAconfig = args.cDNAconfig
        self.oligoconfig = args.oligoconfig
        self.oligotype = args.oligotype
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        self.process = args.process
        self.mtgenes = args.mtgenes
        self.no_introns = args.no_introns
        self.no_bam = args.no_bam
        self.mixseq = args.mixseq
    
    def run(self):
        with open(self.list) as samplelist:
            for line in samplelist:
                lst = line.strip().split('\t')
                name = lst[0]
                cDNAr1 = lst[1].split(';')[0]
                cDNAr2 = lst[1].split(';')[-1]
                oligor1 = lst[2].split(';')[0]
                oligor2 = lst[2].split(';')[-1]
                species = lst[-1]
                shelllist = open('%s.sh'%name,'w')
                path = '/'.join(str(_root_dir).split('/')[0:-4])+ '/bin'
                cmd_line = ['%s/DNBC4tools run --name %s --cDNAfastq1 %s --cDNAfastq2 %s --oligofastq1 %s --oligofastq2 %s --starIndexDir %s --gtf %s --species %s'
                %(path,name,cDNAr1,cDNAr2,oligor1,oligor2,self.starIndexDir,self.gtf,species)]
                if self.thread:
                    cmd_line += ['--thread %s'%self.thread]
                if self.outdir:
                    cmd_line += ['--outdir %s'%self.outdir]
                if self.cDNAconfig:
                    cmd_line += ['--cDNAconfig %s'%self.cDNAconfig]
                if self.oligoconfig:
                    cmd_line += ['--oligoconfig %s'%self.oligoconfig]
                if self.oligotype:
                    cmd_line += ['--oligotype %s'%self.oligotype]
                if self.calling_method:
                    cmd_line += ['--calling_method %s'%self.calling_method]
                if self.expectcells:
                    cmd_line += ['--expectcells %s'%self.expectcells]
                if self.forcecells:
                    cmd_line += ['--forcecells %s'%self.forcecells]
                if self.process:
                    cmd_line += ['--process %s'%self.process]
                if self.mtgenes:
                    cmd_line += ['--mtgenes %s'%self.mtgenes]
                if self.no_introns:
                    cmd_line += ['--no_introns']
                if self.no_bam:
                    cmd_line += ['--no_bam']
                if self.mixseq:
                    cmd_line += ['--mixseq']
                cmd_line = ' '.join(cmd_line)
                shelllist.write(cmd_line + '\n')
                
def multi(args):
    Multi_list(args).run()

def parse_multi(parser):
    parser.add_argument('--list', metavar='FILE',help='sample list.', type=str,required=True)
    parser.add_argument('--starIndexDir',type=str, metavar='PATH',help='Star index dir path.',required=True)
    parser.add_argument('--gtf',type=str, metavar='GTF',help='GTF file.',required=True)
    parser.add_argument('--outdir', metavar='PATH',help='output dir, [default: current directory].')
    parser.add_argument('--thread',type=int, metavar='INT',help='Analysis threads, [defult: 4].')
    parser.add_argument('--cDNAconfig', metavar='JASON',help='whitelist of cell barcode and structure file in JSON format for cDNA fastq, [defalut: DNBelabC4_scRNA_beads_readStructure.json].')
    parser.add_argument('--oligoconfig', metavar='JASON',help='whitelist of cell barcode and structure file in JSON format for oligo fastq, [defalut: DNBelabC4_scRNA_oligo_readStructure.json].')
    parser.add_argument('--oligotype', metavar='FILE',help='Whitelist of oligo index, [default: oligo_type8.txt].')
    parser.add_argument('--calling_method',metavar='STR',help='Cell calling method, Choose from barcoderanks and emptydrops, [default: emptydrops].')
    parser.add_argument('--expectcells',metavar='INT',help='Expected number of recovered beads, [default: 3000].')
    parser.add_argument('--forcecells',metavar='INT',help='Force pipeline to use this number of beads.')
    parser.add_argument('--mtgenes',metavar='LIST',help='Set mitochondrial genes(mtgene list file path) or auto, [default: auto].')
    parser.add_argument('--process', metavar='TEXT',help='Custom your analysis steps, steps are separated by commas, [default: data,count,analysis,report].')
    parser.add_argument('--no_introns', action='store_true',help='Not include intronic reads in count.')
    parser.add_argument('--mixseq', action='store_true',help='cDNA and oligo sequence in same chip.')
    parser.add_argument('--no_bam', action='store_true',help='Do not move filter bam file to output dir.')
    return parser
    
