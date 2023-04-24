import os,collections
from dnbc4tools.tools.utils import str_mkdir,judgeFilexits,change_path,read_json,logging_call
from dnbc4tools.__init__ import __root_dir__

class Runpipe:
    def __init__(self, args):
        self.name = args.name
        self.cDNAr1 = args.cDNAfastq1
        self.cDNAr2 = args.cDNAfastq2
        self.oligor1 = args.oligofastq1
        self.oligor2 = args.oligofastq2
        self.genomeDir = args.genomeDir
        
        self.outdir = args.outdir
        self.threads = args.threads
        self.chemistry = args.chemistry
        self.darkreaction = args.darkreaction
        self.customize = args.customize
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        
        self.process = args.process
        self.no_introns = args.no_introns
        self.minumi = 1000
        
    def runpipe(self):
        change_path()
        genomeDir = os.path.abspath(self.genomeDir)
        judgeFilexits('%s/ref.json'%genomeDir)
        indexConfig = read_json('%s/ref.json'%genomeDir)
        gtf = indexConfig['gtf']

        judgeFilexits(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.genomeDir,gtf)
        data_cmd = ['dnbc4rna data --cDNAfastq1 %s --cDNAfastq2 %s --oligofastq1 %s --oligofastq2 %s --threads %s --name %s --chemistry %s --darkreaction %s --outdir %s --genomeDir %s'
        %(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,self.threads,self.name,self.chemistry,self.darkreaction,self.outdir,self.genomeDir)]
        if self.customize:
            data_cmd += ['--customize %s'%self.customize]
        if self.no_introns:
            data_cmd += ['--no_introns']
        data_cmd = ' '.join(data_cmd)

        count_cmd = 'dnbc4rna count --name %s --calling_method %s --expectcells %s --forcecells %s --minumi %s --threads %s --outdir %s'\
        %(self.name,self.calling_method,self.expectcells,self.forcecells,self.minumi,self.threads,self.outdir)

        analysis_cmd = 'dnbc4rna analysis --name %s  --outdir %s --genomeDir %s'\
            %(self.name,self.outdir,self.genomeDir) 
        
        report_cmd = ['dnbc4rna report --name %s --genomeDir %s --outdir %s --threads %s'
        %(self.name,self.genomeDir,self.outdir,self.threads)]
        if self.no_introns:
            report_cmd += ['--no_introns']
        report_cmd = ' '.join(report_cmd)

       
        pipelist = str(self.process).split(',')
        for pipe in pipelist:
            if pipe not in ['data','count','analysis','report','']:
                print('\033[0;31;40mUnable to recognize pipe!\033[0m')
                raise Exception('Unable to recognize pipe!')
        
        cmdlist = collections.OrderedDict()
        if 'data' in pipelist:
            cmdlist['data'] = data_cmd
        if 'count' in pipelist:
            cmdlist['count'] = count_cmd
        if 'analysis' in pipelist:
            cmdlist['analysis'] = analysis_cmd
        if 'report' in pipelist:
            cmdlist['report'] = report_cmd

        str_mkdir('%s/log'%os.path.join(self.outdir,self.name))
        for pipe,pipecmd in cmdlist.items():
            logging_call(pipecmd,pipe,os.path.join(self.outdir,self.name))
            
def run(args):
    Runpipe(args).runpipe()

def helpInfo_run(parser):
    parser.add_argument(
        '--name', 
        metavar='NAME',
        help='Sample name.', 
        type=str,
        required=True
        )
    parser.add_argument(
        '--cDNAfastq1', 
        metavar='FASTQ',
        help='Paths to the raw R1 fastq files of cDNA library.', 
        required=True
        )
    parser.add_argument(
        '--cDNAfastq2', 
        metavar='FASTQ',
        help='Paths to the raw R2 fastq files of cDNA library.', 
        required=True
        )
    parser.add_argument(
        '--oligofastq1', 
        metavar='FASTQ',
        help='Paths to the raw R1 fastq files of oligo library.',
        required=True
        )
    parser.add_argument(
        '--oligofastq2', 
        metavar='FASTQ',
        help='Paths to the raw R2 fastq files of oligo library.',
        required=True
        )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '--outdir', 
        metavar='PATH',
        help='Output directory, [default: current directory].', 
        default=os.getcwd()
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads used during analysis, [default: 4].'
        )
    parser.add_argument(
        '--calling_method',
        metavar='STR',
        choices=["barcoderanks","emptydrops"],
        help='Cell calling method, choose from barcoderanks and emptydrops, [default: emptydrops].', 
        default='emptydrops'
        )
    parser.add_argument(
        '--expectcells',
        metavar='INT',
        help='Expected number of recovered beads, [default: 3000].', 
        default=3000
        )
    parser.add_argument(
        '--forcecells',
        metavar='INT',
        help='Force pipeline to use this number of beads.', 
        default=0
        )
    parser.add_argument(
        '--chemistry',
        metavar='STR',
        choices=["scRNAv1HT","scRNAv2HT","auto"],
        help='Chemistry version. Automatic detection is recommended , [default: auto].',
        default='auto'
        )
    parser.add_argument(
        '--darkreaction',
        metavar='STR',
        help='Sequencing dark cycles. Automatic detection is recommended, [default: auto].', 
        default='auto'
        )
    parser.add_argument(
        '--customize',
        metavar='STR',
        help='Customize files for whitelist and readstructure in JSON format for cDNA and oligo.'
        )
    parser.add_argument(
        '--process', 
        metavar='STR',
        help='Custom analysis steps allow skipping of certain analysis steps that are not needed to be re-run, [default: data,count,analysis,report].',
        type=str,
        default='data,count,analysis,report'
        )
    parser.add_argument(
        '--no_introns', 
        action='store_true',
        help='Intron reads are not included in the expression matrix.'
        )
    return parser