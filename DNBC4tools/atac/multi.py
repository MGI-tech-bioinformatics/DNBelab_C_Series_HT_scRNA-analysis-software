from typing import Optional
from dnbc4tools.__init__ import __root_dir__

class MultiList:
    def __init__(self, args):
        """
        Initialize the MultiList class.

        :param list_file: str, the path to the sample list file.
        :param genome_dir: str, the path to the genome directory.
        :param outdir: Optional[str], the path to the output directory.
        :param threads: Optional[int], the number of threads to use.
        """
        self.list = args.list
        self.genomeDir = args.genomeDir
        self.outdir = args.outdir
        self.threads = args.threads
    
    def run(self) -> None:
        """
        Run the dnbc4atac pipeline for all samples listed in the sample list file.
        """
        with open(self.list) as samplelist:
            for line in samplelist:
                lst = line.strip().split('\t')
                name = lst[0]
                fastqr1 = lst[1].split(';')[0]
                fastqr2 = lst[1].split(';')[-1]
                shelllist = open('%s.sh'%name,'w')
                path = '/'.join(str(__root_dir__).split('/')[0:-4])+ '/bin'
                cmd_line = ['%s/dnbc4atac run --name %s --fastq1 %s --fastq2 %s --genomeDir %s'
                %(path,name,fastqr1,fastqr2,self.genomeDir)]
                if self.threads:
                    cmd_line += ['--threads %s'%self.threads]
                if self.outdir:
                    cmd_line += ['--outdir %s'%self.outdir]
                cmd_line = ' '.join(cmd_line)
                shelllist.write(cmd_line + '\n')
                
def multi(args):
    MultiList(args).run()

def helpInfo_multi(parser):
    parser.add_argument(
        '--list', 
        metavar='FILE',
        help='sample list.', 
        type=str,
        required=True
        )
    parser.add_argument(
        '--outdir', 
        metavar='PATH',
        help='Output diretory, [default: current directory].', 
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads used for the analysis, [default: 4].'
        )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    return parser