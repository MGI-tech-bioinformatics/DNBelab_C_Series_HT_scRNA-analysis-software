import sys,argparse
import textwrap,importlib

from DNBC4tools.__init__ import _version, _pipelist

def pipeline_package(pipe):
    package = importlib.import_module(f"DNBC4tools.tools.{pipe}")
    return package

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
        DNBC4tools contain the following pipelines and functions:
        --------------------------------
        \033[0;32;40mpipeline\033[0m:
        \033[0;32;40mDNBC4tools run\033[0m      Run data, count, analysis, report for a complete pipeline.
        \033[0;32;40mDNBC4tools data\033[0m     Run quality control and filtering on the raw fastq, 
                            align the cDNA data to the reference genome and annotate it with gtf file.
        \033[0;32;40mDNBC4tools count\033[0m    Determine the inflection point and judge the empty droplets, 
                            merge multiple beads in the same droplet, calculate the cell * gene expression matrix.
        \033[0;32;40mDNBC4tools analysis\033[0m Run quality control on the cell expression matrix, filter low-quality cells and genes,
                            perform cell clustering analysis and marker gene screening based on the expression matrix.
        \033[0;32;40mDNBC4tools report\033[0m   Summary result files and generate html report.
        
        \033[0;32;40mfunction\033[0m:
        \033[0;32;40mDNBC4tools multi\033[0m    Creat shell list for multi samples.
        \033[0;32;40mDNBC4tools mkref\033[0m    Create a genome reference directory.
        \033[0;32;40mDNBC4tools clean\033[0m    Delete the temp files.'''))
    parser.add_argument('-v', '--version', action='version', version=_version)
    subparsers = parser.add_subparsers(dest='parser_step')

    for _pipe in _pipelist:
        package= pipeline_package(_pipe)
        steps = getattr(package, _pipe)
        steps_help = getattr(package, f"parse_{_pipe}")
        parser_step = subparsers.add_parser(_pipe, description='DNBC4tools %s'%_pipe,formatter_class=argparse.RawDescriptionHelpFormatter)
        steps_help(parser_step)
        parser_step.set_defaults(func=steps)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()