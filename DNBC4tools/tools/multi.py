import os
from DNBC4tools.__init__ import _root_dir

class Multi_list:
    def __init__(self, args):
        self.list = args.list
        self.genomeDir = args.genomeDir
        self.gtf = args.gtf
        self.outdir = args.outdir
        self.thread = args.thread
        self.chemistry = args.chemistry
        self.darkreaction = args.darkreaction
        self.customize = args.customize
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        self.process = args.process
        self.mtgenes = args.mtgenes
        self.no_introns = args.no_introns
        self.no_bam = args.no_bam
    
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
                cmd_line = ['%s/DNBC4tools run --name %s --cDNAfastq1 %s --cDNAfastq2 %s --oligofastq1 %s --oligofastq2 %s --genomeDir %s --gtf %s --species %s'
                %(path,name,cDNAr1,cDNAr2,oligor1,oligor2,self.genomeDir,self.gtf,species)]
                if self.thread:
                    cmd_line += ['--thread %s'%self.thread]
                if self.outdir:
                    cmd_line += ['--outdir %s'%self.outdir]
                if self.chemistry:
                    cmd_line += ['--chemistry %s'%self.chemistry]
                if self.darkreaction:
                    cmd_line += ['--darkreaction %s'%self.darkreaction]
                if self.customize:
                    cmd_line += ['--customize %s'%self.customize]
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
                cmd_line = ' '.join(cmd_line)
                shelllist.write(cmd_line + '\n')
                
def multi(args):
    Multi_list(args).run()

def parse_multi(parser):
    parser.add_argument('--list', metavar='FILE',help='sample list.', type=str,required=True)
    parser.add_argument('--genomeDir',type=str, metavar='PATH',help='Path of folder containing reference index.',required=True)
    parser.add_argument('--gtf',type=str, metavar='GTF',help='Path of gtf file.',required=True)
    parser.add_argument('--outdir', metavar='PATH',help='Output diretory, [default: current directory].')
    parser.add_argument('--thread',type=int, metavar='INT',default=4,help='Number of threads to use, [default: 4].')
    parser.add_argument('--calling_method',metavar='STR',choices=["barcoderanks","emptydrops"],help='Cell calling method, Choose from barcoderanks and emptydrops, [default: emptydrops].')
    parser.add_argument('--expectcells',metavar='INT',help='Expected number of recovered beads, [default: 3000].')
    parser.add_argument('--forcecells',metavar='INT',help='Force pipeline to use this number of beads.')
    parser.add_argument('--chemistry',metavar='STR',choices=["scRNAv1HT","scRNAv2HT","auto"],help='Chemistry version. Automatic detection is recommended. If setting, needs to be used with --darkreaction, can be "scRNAv1HT", "scRNAv2HT", [default: auto].')
    parser.add_argument('--darkreaction',metavar='STR',help='Sequencing dark reaction. Automatic detection is recommended. If setting, needs to be used with --chemistry, use comma to separate cDNA and oligo, can be "R1,R1R2", "R1,R1", "unset,unset", etc, [default: auto].')
    parser.add_argument('--customize',metavar='STR',help='Customize files for whitelist and readstructure in JSON format for cDNA and oligo, use comma to separate cDNA and oligo.')
    parser.add_argument('--process', metavar='STR',help='Custom your analysis steps, use comma to separate steps, [default: data,count,analysis,report].',type=str)
    parser.add_argument('--mtgenes',metavar='FILE',help='Path of file with mitochondrial genes, [default: auto].')
    parser.add_argument('--no_introns', action='store_true',help='Not include intronic reads in count.')
    parser.add_argument('--no_bam', action='store_true',help='Do not move filter bam file to the output dir.')
    return parser
    
