#!/usr/bin/env python

"""
A helper script to generate synthetic data for a project in PEP format.
"""

import sys

from argparse import ArgumentParser

from ngs_toolkit.demo import generate_project
from ngs_toolkit.utils import filter_kwargs_by_callable


def parse_arguments():
    """
    Argument Parsing.
    """
    import inspect

    parser = ArgumentParser(
        prog="python -m ngs_toolkit.recipes.generate_project", description=__doc__)

    sig = inspect.signature(generate_project)
    for arg in sig.parameters:
        if arg in ["kwargs", "initialize"]:
            continue
        d = sig.parameters[arg].default
        if d is None:
            parser.add_argument("--" + arg.replace("_", "-"))
        else:
            parser.add_argument(
                "--" + arg.replace("_", "-"),
                default=d, type=type(d))
    parser.add_argument("--debug", action="store_true")
    return parser


def main(cli=None):
    """Generate synthetic data for a project in PEP format."""
    args = parse_arguments().parse_args(cli)
    if args.debug:
        print(args)
    kwargs = {k: v for k, v in args.__dict__.items() if v is not None}
    kwargs = filter_kwargs_by_callable(kwargs, generate_project)
    if args.debug:
        print(kwargs)
    pep = generate_project(**kwargs, initialize=False)
    sys.stdout.write(pep + "\n")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
