import os
from typing import List, Optional
from dnbc4tools.tools.utils import str_mkdir, judgeFilexits, change_path, start_print_cmd
from dnbc4tools.__init__ import __root_dir__

class Runpipe:
    def __init__(self, args):
        self.name: str = args.name
        self.fastq1: str = args.fastq1 
        self.fastq2: str = args.fastq2
        self.threads: int = args.threads
        self.darkreaction: Optional[str] = args.darkreaction
        self.customize: Optional[str] = args.customize
        self.outdir: str = args.outdir
        self.genomeDir: str = args.genomeDir
        self.forcebeads: int = args.forcebeads
        self.forcefrags: int = args.forcefrags
        self.process: List[str] = args.process.split(",")
        
    def runpipe(self) -> None:
        change_path()
        judgeFilexits(self.fastq1, self.fastq2, self.genomeDir)

        data_cmd: List[str] = [
            f'dnbc4atac data --fastq1 %s --fastq2 %s --threads %s --name %s --darkreaction %s --outdir %s --genomeDir %s --bcerror 1'
            %(self.fastq1, self.fastq2, self.threads, self.name, self.darkreaction, self.outdir, self.genomeDir)
            ]
        if self.customize:
            data_cmd += [
                f'--customize %s'
                % self.customize
                ]
        data_cmd = ' '.join(data_cmd)

        decon_cmd: List[str] = [
            f'dnbc4atac decon --name %s --threads %s --outdir %s --genomeDir %s' 
            %(self.name, self.threads, self.outdir, self.genomeDir)
            ]
        if self.forcefrags:
            decon_cmd += [
                f'--forcefrags %s'
                % self.forcefrags
                ]
        if self.forcebeads:
            decon_cmd += [
                f'--forcebeads %s'
                % self.forcebeads
                ]
        decon_cmd = ' '.join(decon_cmd)

        analysis_cmd: str = f'dnbc4atac analysis --name %s --outdir %s --genomeDir %s'\
            %(self.name, self.outdir, self.genomeDir)
        report_cmd: str = f'dnbc4atac report --name %s --outdir %s --genomeDir %s'\
            %(self.name, self.outdir, self.genomeDir)

        cmdlist: List[str] = []
        if 'data' in self.process:
            cmdlist.append(data_cmd)
        if 'decon' in self.process:
            cmdlist.append(decon_cmd)
        if 'analysis' in self.process:
            cmdlist.append(analysis_cmd)
        if 'report' in self.process:
            cmdlist.append(report_cmd)

        str_mkdir(os.path.join(self.outdir, self.name))
        for pipecmd in cmdlist:
            start_print_cmd(pipecmd)
            
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
        '--fastq1', 
        metavar='FASTQ',
        help='The input R1 fastq files.', 
        required=True
        )
    parser.add_argument(
        '--fastq2', 
        metavar='FASTQ',
        help='The input R2 fastq files.', 
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
        help='Output diretory, [default: current directory].', 
        default=os.getcwd()
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads used for the analysis, [default: 4].'
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
        help='Customize readstructure.'
        )
    parser.add_argument(
        '--forcebeads', 
        type=int,
        metavar='INT',
        help='Top N number of beads to be thresholded.'
        )
    parser.add_argument(
        '--forcefrags', 
        type=int, 
        metavar='INT',
        help='Minimum number of fragments to be thresholded.'
        )
    parser.add_argument(
        '--process', 
        metavar='STR',
        help='Custom analysis steps allow skipping of certain analysis steps that are not needed to be re-run, [default: data,decon,analysis,report].',
        type=str,
        default='data,decon,analysis,report')
    return parser