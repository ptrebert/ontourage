#!/usr/bin/env python3

import sys
import argparse as argp
import logging
import traceback

# JUST FOR DEV/DBG
sys.path.insert(0, '/home/local/work/code/github/ontourage')

from ontourage import LOG_MESSAGE_FORMAT as logging_format
from ontourage.process_graph import process_gfa_cli_parser
from ontourage.process_align import process_gaf_cli_parser


def parse_command_line():

    parser = argp.ArgumentParser(prog='ONTourage', add_help=True)

    # parser.add_argument("--version", "-v", action="version", version=__version__)

    noise_level = parser.add_mutually_exclusive_group(required=False)
    noise_level.add_argument(
        "--debug",
        "-dbg",
        action="store_true",
        default=False,
        help="Print debug log messages to stderr.",
    )
    noise_level.add_argument(
        "--verbose",
        "-vrb",
        action="store_true",
        default=False,
        help="Print progress messages to stdout.",
    )

    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        dest='num_jobs',
        default=1,
        help="Number of CPU cores to use. Default: 1 (no sanity checks!)",
    )

    parser = add_sub_parsers(parser)

    return parser


def add_sub_parsers(main_parser):
    subparsers = main_parser.add_subparsers(
        title="Commands",
        dest="command_name",
        description="Valid ONTourage commands"
    )
    subparsers = process_gfa_cli_parser(subparsers)
    subparsers = process_gaf_cli_parser(subparsers)
    return main_parser


if __name__ == "__main__":
    parser = parse_command_line()
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(stream=sys.stderr, level=logging.DEBUG, format=logging_format)
    elif args.verbose:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=logging_format)
    else:
        logging.basicConfig(stream=sys.stderr, level=logging.WARNING, format=logging_format)
    logger = logging.getLogger(None)
    logger.info("Logging system initialized")
    try:
        args.execute(args)
    except AttributeError as attr_error:
        if 'execute' in str(attr_error):
            parser.print_help()
            sys.stderr.write('\nPlease select a valid command to execute.\n')
        else:
            raise
    except Exception:
        traceback.print_exc()
        raise
