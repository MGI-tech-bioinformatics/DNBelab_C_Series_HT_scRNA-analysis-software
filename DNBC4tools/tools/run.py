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
        self.genomeDir = args.genomeDir
        self.gtf = args.gtf
        self.species = args.species
        
        self.outdir = args.outdir
        self.thread = args.thread
        self.chemistry = args.chemistry
        self.darkreaction = args.darkreaction
        self.customize = args.customize
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        self.mtgenes = args.mtgenes
        
        self.process = args.process
        self.no_introns = args.no_introns
        self.no_bam = args.no_bam
        self.dry = args.dry
        
    def runpipe(self):
        change_path()
        judgeFilexits(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.genomeDir,self.gtf)
        data_cmd = ['DNBC4tools data --cDNAfastq1 %s --cDNAfastq2 %s --oligofastq1 %s --oligofastq2 %s --thread %s --name %s --chemistry %s --darkreaction %s --outdir %s --genomeDir %s --gtf %s'
        %(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.thread,self.name,self.chemistry,self.darkreaction,self.outdir,self.genomeDir,self.gtf)]
        if self.customize:
            data_cmd += ['--customize %s'%self.customize]
        if self.no_introns:
            data_cmd += ['--no_introns']
        data_cmd = ' '.join(data_cmd)

        count_cmd = 'DNBC4tools count --name %s --raw_matrix %s/%s/01.data/raw_matrix --bam %s/%s/01.data/final_sorted.bam --calling_method %s --expectcells %s --forcecells %s --minumi 1000 --cDNAbarcodeCount %s/%s/01.data/cDNA_barcode_counts_raw.txt --Indexreads %s/%s/01.data/Index_reads.fq.gz --oligobarcodeCount %s/%s/01.data/Index_barcode_counts_raw.txt --thread %s --outdir %s'\
        %(self.name,self.outdir,self.name,self.outdir,self.name,self.calling_method,self.expectcells,self.forcecells,self.outdir,self.name,self.outdir,self.name,self.outdir,self.name,self.thread,self.outdir)

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
    parser.add_argument('--name', metavar='NAME',help='Sample name.', type=str,required=True)
    parser.add_argument('--cDNAfastq1', metavar='FASTQ',help='cDNA R1 fastq file, use commas to separate multiple files.', required=True)
    parser.add_argument('--cDNAfastq2', metavar='FASTQ',help='cDNA R2 fastq file, use commas to separate multiple files, the files order needs to be consistent with cDNAfastq1.', required=True)
    parser.add_argument('--oligofastq1', metavar='FASTQ',help='oligo R1 fastq file, use commas to separate multiple files.',required=True)
    parser.add_argument('--oligofastq2', metavar='FASTQ',help='oligo R2 fastq file, use commas to separate multiple files, the files order needs to be consistent with oligofastq1.',required=True)
    parser.add_argument('--genomeDir',type=str, metavar='PATH',help='Path of folder containing reference index.',required=True)
    parser.add_argument('--gtf',type=str, metavar='GTF',help='Path of gtf file.',required=True)
    parser.add_argument('--species',type=str, metavar='STR',default='undefined',help='Species name. Only "Homo_sapiens","Human","Mus_musculus" and "Mouse" can perform cell annotation analysis, [default: undefined].')
    parser.add_argument('--outdir', metavar='PATH',help='Output diretory, [default: current directory].', default=os.getcwd())
    parser.add_argument('--thread',type=int, metavar='INT',default=4,help='Number of threads to use, [default: 4].')
    parser.add_argument('--calling_method',metavar='STR',choices=["barcoderanks","emptydrops"],help='Cell calling method, Choose from barcoderanks and emptydrops, [default: emptydrops].', default='emptydrops')
    parser.add_argument('--expectcells',metavar='INT',help='Expected number of recovered beads, [default: 3000].', default=3000)
    parser.add_argument('--forcecells',metavar='INT',help='Force pipeline to use this number of beads.', default=0)
    parser.add_argument('--chemistry',metavar='STR',choices=["scRNAv1HT","scRNAv2HT","auto"],help='Chemistry version. Automatic detection is recommended. If setting, needs to be used with --darkreaction, can be "scRNAv1HT", "scRNAv2HT", [default: auto].',default='auto')
    parser.add_argument('--darkreaction',metavar='STR',help='Sequencing dark reaction. Automatic detection is recommended. If setting, needs to be used with --chemistry, use comma to separate cDNA and oligo, can be "R1,R1R2", "R1,R1", "unset,unset", etc, [default: auto].', default='auto')
    parser.add_argument('--customize',metavar='STR',help='Customize files for whitelist and readstructure in JSON format for cDNA and oligo, use comma to separate cDNA and oligo.')
    parser.add_argument('--process', metavar='STR',help='Custom your analysis steps, use comma to separate steps, [default: data,count,analysis,report].',type=str,default='data,count,analysis,report')
    parser.add_argument('--mtgenes',metavar='FILE',default='auto',help='Path of file with mitochondrial genes, [default: auto].')
    parser.add_argument('--no_introns', action='store_true',help='Not include intronic reads in count.')
    parser.add_argument('--no_bam', action='store_true',help='Do not move filter bam file to the output dir.')
    parser.add_argument('--dry', help=' Do not execute the pipeline. Generate a pipeline shell file.',action='store_true')
    return parser