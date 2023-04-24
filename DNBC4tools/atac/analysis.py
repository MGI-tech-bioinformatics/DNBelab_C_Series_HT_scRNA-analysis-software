import os,argparse
from dnbc4tools.tools.utils import str_mkdir,judgeFilexits,change_path,logging_call,read_json,bin_path
from dnbc4tools.__init__ import __root_dir__

class Analysis:
    def __init__(self, args: argparse.Namespace) -> None:
        """
        Initialize the Analysis class.

        Args:
        - args (argparse.Namespace): parsed command-line arguments
        """
        self.name = args.name
        self.outdir = os.path.join(args.outdir, args.name)
        self.genomeDir = args.genomeDir

    def run(self) -> None:
        """
        Run the analysis.
        """
        judgeFilexits(self.genomeDir)
        str_mkdir('%s/03.analysis/peak'%self.outdir)
        # str_mkdir('%s/03.analysis/promoter'%self.outdir)
        str_mkdir('%s/03.analysis/images'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        change_path()
        bin_command = bin_path()

        genomeDir = os.path.abspath(self.genomeDir)
        indexConfig = read_json('%s/ref.json'%genomeDir)
        tss = indexConfig['tss']
        chrmt = indexConfig['chrmt']
        genomesize = indexConfig['genomesize']
        species = indexConfig['species']
        # promoter = indexConfig['promoter']

        macs2_cmd = (
            '%s/macs2 callpeak -t %s/02.decon/%s.fragments.tsv.gz -f BED -g %s -n %s -B -q 0.001 --nomodel --outdir %s/03.analysis'
            %(bin_command,self.outdir,self.name,genomesize,self.name,self.outdir)
        )

        cluster_cmd = (
            '%s/Rscript %s/atac/src/Cluster_Annotation.R -I %s/03.analysis/%s_peaks.narrowPeak -F %s/02.decon/%s.fragments.tsv.gz -T %s -MT %s -Q %s/02.decon/%s.Metadata.tsv -O %s/03.analysis -S %s'
            %(bin_command,__root_dir__,self.outdir,self.name,self.outdir,self.name,tss,chrmt,self.outdir,self.name,self.outdir,species)
        )
        
        print('\nPeak calling.')
        logging_call(macs2_cmd, 'analysis', self.outdir)
        print('\nDimensionality reduction, clustering')
        logging_call(cluster_cmd, 'analysis', self.outdir)

def analysis(args: argparse.Namespace) -> None:
    """
    Run the analysis.

    Args:
    - args (argparse.Namespace): parsed command-line arguments
    """
    Analysis(args).run()


def helpInfo_analysis(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add command-line arguments for the analysis subcommand.

    Args:
    - parser (argparse.ArgumentParser): argparse parser

    Returns:
    - argparse.ArgumentParser: argparse parser
    """
    parser.add_argument(
        '--name', 
        metavar='NAME',
        help='Sample name.', 
        type=str,
        required=True
    )
    parser.add_argument(
        '--outdir', 
        metavar='PATH',
        help='Output directory, [default: current directory].', 
        default=os.getcwd()
    )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path of folder containing reference database.',
        required=True
    )
    return parser