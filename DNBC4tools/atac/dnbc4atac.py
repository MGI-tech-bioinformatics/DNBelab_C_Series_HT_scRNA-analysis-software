import sys
import argparse
import textwrap
import importlib
from typing import Any

from dnbc4tools.__init__ import __version__
from dnbc4tools.atac.__init__ import _pipe
from dnbc4tools.tools.text import help_text, help_sub_text


def pipeline_package(pipe: str) -> Any:
    package = importlib.import_module("dnbc4tools.atac.%s" % pipe)
    return package


def main() -> None:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(help_text('atac', 'dnbc4atac'))
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=__version__
    )
    subparsers = parser.add_subparsers(dest='parser_step')

    for pipe in _pipe:
        package = pipeline_package(pipe)
        steps = getattr(package, pipe)
        steps_help = getattr(package, "helpInfo_%s" % pipe)
        parser_step = subparsers.add_parser(
            pipe,
            description=textwrap.dedent(help_sub_text('atac', 'dnbc4atac', pipe)),
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
        steps_help(parser_step)
        parser_step.set_defaults(func=steps)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()