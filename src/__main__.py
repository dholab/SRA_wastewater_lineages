#!/usr/bin/env python3

"""
TODO
"""

from loguru import logger

from .cli import parse_command_line_args


def main() -> None:
    """
    TODO
    """
    logger.debug("Hi mom!")
    _ = parse_command_line_args()


if __name__ == "__main__":
    main()
