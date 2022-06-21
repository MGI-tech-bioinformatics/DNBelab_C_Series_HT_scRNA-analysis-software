import os
from .utils import str_mkdir,judgeFilexits,start_print_cmd,change_path
from DNBC4tools.__init__ import _root_dir

class Runpipe:
    def __init__(self, args):
        self.name = args.name
        self.cDNAr1 = args.cDNAfastq1
        self.cDNAr2 = args.cDNAfastq2
        self.oligor1 = args.oligofastq1
        self.oligor2 = args.oligofastq2
        self.starIndexDir = args.starIndexDir
        self.gtf = args.gtf
        self.species = args.species
        
        self.outdir = args.outdir
        self.thread = args.thread
        self.cDNAconfig = args.cDNAconfig
        self.oligoconfig = args.oligoconfig
        self.oligotype = args.oligotype
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        self.mtgenes = args.mtgenes
        
        self.process = args.process
        self.no_introns = args.no_introns
        self.no_bam = args.no_bam
        self.mixseq = args.mixseq
        self.dry = args.dry
        
    def runpipe(self):
        change_path()
        judgeFilexits(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.cDNAconfig,self.oligoconfig,self.starIndexDir,self.gtf,self.oligotype)
        data_cmd = ['DNBC4tools data --cDNAfastq1 %s --cDNAfastq2 %s --oligofastq1 %s --oligofastq2 %s --thread %s --name %s --cDNAconfig %s --oligoconfig %s --outdir %s --starIndexDir %s --gtf %s'
        %(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.thread,self.name,self.cDNAconfig,self.oligoconfig,self.outdir,self.starIndexDir,self.gtf)]
        if self.mixseq:
            data_cmd += ['--mixseq']
        if self.no_introns:
            data_cmd += ['--no_introns']
        data_cmd = ' '.join(data_cmd)

        count_cmd = 'DNBC4tools count --name %s --raw_matrix %s/%s/01.data/raw_matrix --bam %s/%s/01.data/final.bam --calling_method %s --expectcells %s --forcecells %s --minumi 1000 --cDNAbarcodeCount %s/%s/01.data/cDNA_barcode_counts_raw.txt --Indexreads %s/%s/01.data/Index_reads.fq.gz --oligobarcodeCount %s/%s/01.data/Index_barcode_counts_raw.txt --thread %s --oligotype %s --outdir %s'\
        %(self.name,self.outdir,self.name,self.outdir,self.name,self.calling_method,self.expectcells,self.forcecells,self.outdir,self.name,self.outdir,self.name,self.outdir,self.name,self.thread,self.oligotype,self.outdir)

        analysis_cmd = 'DNBC4tools analysis --name %s --matrix %s/%s/02.count/filter_matrix --species %s --outdir %s --mtgenes %s'\
            %(self.name,self.outdir,self.name,self.species,self.outdir,self.mtgenes) 
        
        report_cmd = ['DNBC4tools report --name %s --species %s --outdir %s --thread %s'
        %(self.name,self.species,self.outdir,self.thread)]
        if self.no_bam:
            report_cmd += ['--no_bam']
        if self.no_introns:
            report_cmd += ['--no_introns']
        report_cmd = ' '.join(report_cmd)

       
        pipelist = str(self.process).split(',')
        cmdlist = []
        if 'data' in pipelist:
            cmdlist.append(data_cmd)
        if 'count' in pipelist:
            cmdlist.append(count_cmd)
        if 'analysis' in pipelist:
            cmdlist.append(analysis_cmd)
        if 'report' in pipelist:
            cmdlist.append(report_cmd)
        if not self.dry:
            str_mkdir('%s'%os.path.join(self.outdir,self.name))
            for pipecmd in cmdlist:
                start_print_cmd(pipecmd)
        else:
            shelllist = open('%s/%s.shell'%(self.outdir,self.name),'w')
            for pipecmd in cmdlist:
                shelllist.write(pipecmd+'\n')
            shelllist.close()
            
def run(args):
    Runpipe(args).runpipe()

def parse_run(parser):
    parser.add_argument('--name', metavar='NAME',help='sample name.', type=str,required=True)
    parser.add_argument('--cDNAfastq1', metavar='FASTQ',help='cDNA R1 fastq file, Multiple files are separated by commas.', required=True)
    parser.add_argument('--cDNAfastq2', metavar='FASTQ',help='cDNA R2 fastq file, Multiple files are separated by commas, the files order needs to be consistent with cDNAfastq1.', required=True)
    parser.add_argument('--oligofastq1', metavar='FASTQ',help='oligo R1 fastq file, Multiple files are separated by commas.',required=True)
    parser.add_argument('--oligofastq2', metavar='FASTQ',help='oligo R2 fastq file, Multiple files are separated by commas, the files order needs to be consistent with oligofastq1.',required=True)
    parser.add_argument('--starIndexDir',type=str, metavar='PATH',help='Star index dir path.',required=True)
    parser.add_argument('--gtf',type=str, metavar='GTF',help='GTF file.',required=True)
    parser.add_argument('--species',type=str, metavar='STR',default='NA',help='Species name. Only [Human] and [Mouse] can analysis cell annotation, [default: NA].')
    parser.add_argument('--outdir', metavar='PATH',help='output dir, [default: current directory].', default=os.getcwd())
    parser.add_argument('--thread',type=int, metavar='INT',default=4,help='Analysis threads, [default: 4].')
    parser.add_argument('--cDNAconfig', metavar='JASON',help='whitelist of cell barcode and structure file in JSON format for cDNA fastq, [defalut: DNBelabC4_scRNA_beads_readStructure.json].',default='%s/config/DNBelabC4_scRNA_beads_readStructure.json'%_root_dir)
    parser.add_argument('--oligoconfig', metavar='JASON',help='whitelist of cell barcode and structure file in JSON format for oligo fastq, [defalut: DNBelabC4_scRNA_oligo_readStructure.json].',default='%s/config/DNBelabC4_scRNA_oligo_readStructure.json'%_root_dir)
    parser.add_argument('--oligotype', metavar='FILE',help='Whitelist of oligo index, [default: oligo_type8.txt].',default='%s/config/oligo_type8.txt'%_root_dir)
    parser.add_argument('--calling_method',metavar='STR',help='Cell calling method, Choose from barcoderanks and emptydrops, [default: emptydrops].', default='emptydrops')
    parser.add_argument('--expectcells',metavar='INT',help='Expected number of recovered beads, [default: 3000].', default=3000)
    parser.add_argument('--forcecells',metavar='INT',help='Force pipeline to use this number of beads.', default=0)
    parser.add_argument('--mtgenes',metavar='LIST',default='auto',help='Set mitochondrial genes(mtgene list file path) or auto, [default: auto].')
    parser.add_argument('--process', metavar='TEXT',help='Custom your analysis steps, steps are separated by commas, [default: data,count,analysis,report].',type=str,default='data,count,analysis,report')
    parser.add_argument('--no_introns', action='store_true',help='Not include intronic reads in count.')
    parser.add_argument('--mixseq', action='store_true',help='cDNA and oligo sequence in same chip.')
    parser.add_argument('--no_bam', action='store_true',help='Do not move filter bam file to output dir.')
    parser.add_argument('--dry', help=' Do not execute the pipeline. Generate a pipeline shell file.',action='store_true')
    return parser