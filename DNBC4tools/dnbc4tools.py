import sys,argparse
import textwrap,importlib

from dnbc4tools.__init__ import __version__,__category__
from dnbc4tools.tools.text import help_text,sum_help,help_sub_text

def category_pipe(pipe):
    package = importlib.import_module("dnbc4tools.%s.__init__"%pipe)
    pipelist = package._pipe
    return pipelist

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(sum_help)
        )
    parser.add_argument(
        '-v', 
        '--version', 
        action='version', 
        version=__version__
        )
    subparsers = parser.add_subparsers(dest='parser_step')

    for _category in __category__:
        sub_subparser = subparsers.add_parser(
            _category, 
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(help_text(_category,'dnbc4tools %s'%_category))
            )
        sub_subparser.add_argument(
            '-v', 
            '--version', 
            action='version', 
            version=__version__
            )
        subsub_subparser= sub_subparser.add_subparsers()
        _pipeList = category_pipe(_category)

        for _pipe in _pipeList:
            package = importlib.import_module(
                "dnbc4tools.%s.%s"%(_category,_pipe)
                )
            pipes = getattr(package, _pipe)
            pipes_help = getattr(package, "helpInfo_%s"%_pipe)
            parser_step = subsub_subparser.add_parser(
                _pipe, 
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=textwrap.dedent(help_sub_text(_category,'dnbc4tools %s'%_category,_pipe))
                )
            pipes_help(parser_step)
            parser_step.set_defaults(func=pipes)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()