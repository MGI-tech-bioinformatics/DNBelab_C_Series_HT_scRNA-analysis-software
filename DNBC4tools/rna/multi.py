import os
from dnbc4tools.__init__ import __root_dir__

class Multi_list:
    def __init__(self, args):
        self.list = args.list
        self.genomeDir = args.genomeDir
        self.outdir = args.outdir
        self.threads = args.threads
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.no_introns = args.no_introns
    
    def run(self):
        with open(self.list) as samplelist:
            for line in samplelist:
                lst = line.strip().split('\t')
                name = lst[0]
                cDNAr1 = lst[1].split(';')[0]
                cDNAr2 = lst[1].split(';')[-1]
                oligor1 = lst[2].split(';')[0]
                oligor2 = lst[2].split(';')[-1]
                shelllist = open('%s.sh'%name,'w')
                path = '/'.join(str(__root_dir__).split('/')[0:-4])+ '/bin'
                cmd_line = ['%s/dnbc4rna run --name %s --cDNAfastq1 %s --cDNAfastq2 %s --oligofastq1 %s --oligofastq2 %s --genomeDir %s'
                %(path,name,cDNAr1,cDNAr2,oligor1,oligor2,self.genomeDir)]
                if self.threads:
                    cmd_line += ['--threads %s'%self.threads]
                if self.outdir:
                    cmd_line += ['--outdir %s'%self.outdir]
                if self.calling_method:
                    cmd_line += ['--calling_method %s'%self.calling_method]
                if self.expectcells:
                    cmd_line += ['--expectcells %s'%self.expectcells]
                if self.no_introns:
                    cmd_line += ['--no_introns']
                cmd_line = ' '.join(cmd_line)
                shelllist.write(cmd_line + '\n')
                
def multi(args):
    Multi_list(args).run()

def helpInfo_multi(parser):
    parser.add_argument(
        '--list', 
        metavar='FILE',
        help='sample list.', 
        type=str,
        required=True
        )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path of folder containing reference database.',
        required=True
        )
    parser.add_argument(
        '--outdir', 
        metavar='PATH',
        help='Output directory, [default: current directory].'
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads used for the analysis, [default: 4].'
        )
    parser.add_argument(
        '--calling_method',
        metavar='STR',
        choices=["barcoderanks","emptydrops"],
        help='Cell calling method, choose from barcoderanks and emptydrops, [default: emptydrops].'
        )
    parser.add_argument(
        '--expectcells',
        metavar='INT',
        help='Expected number of recovered beads, [default: 3000].'
        )
    parser.add_argument(
        '--no_introns', 
        action='store_true',
        help='Intron reads are not included in the expression matrix.'
        )
    return parser
    
