import sys,argparse
import textwrap,importlib

from dnbc4tools.__init__ import __version__
from dnbc4tools.rna.__init__ import _pipe
from dnbc4tools.tools.text import help_text,help_sub_text

def pipeline_package(pipe):
    package = importlib.import_module("dnbc4tools.rna.%s"%pipe)
    return package

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent(help_text('rna','dnbc4rna')))
    parser.add_argument(
        '-v', 
        '--version', 
        action='version', 
        version=__version__
        )
    subparsers = parser.add_subparsers(dest='parser_step')

    for pipe in _pipe:
        package= pipeline_package(pipe)
        steps = getattr(package, pipe)
        steps_help = getattr(package, "helpInfo_%s"%pipe)
        parser_step = subparsers.add_parser(
            pipe, 
            description=textwrap.dedent(help_sub_text('rna','dnbc4rna',pipe)),
            formatter_class=argparse.RawDescriptionHelpFormatter
            )
        steps_help(parser_step)
        parser_step.set_defaults(func=steps)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()