"""
TODO
"""

import argparse


def parse_command_line_args() -> argparse.Namespace:
    """
    TODO
    """
    parser = argparse.ArgumentParser()

    # set up some subparsers
    subparsers = parser.add_subparsers(title="Subcommands", dest="subcommands")
    fetch = subparsers.add_parser("fetch", help="Fetch SRA accession dataset(s).")
    _haplotype = subparsers.add_parser(
        "haplo",
        help="Identify unique haplotypes within SRA reads.",
    )
    _extract = subparsers.add_parser(
        "extract",
        help="Extract mutations present in each SRA haplotype.",
    )

    # add arguments for the fetch subcommand
    fetch.add_argument(
        "--SRA",
        type=str,
        help="SRA number must be case-sensitive perfect match",
        required=True,
    )

    return parser.parse_args()
