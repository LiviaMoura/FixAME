##!/usr/bin/env python3

# **********************************
# **********************************
# *** A mente é seu pior inimigo ***
# **********************************
# **********************************

__author__ = "Livia Moura"
__copyright__ = "Copyright 2021"
__maintainer__ = "Livia Moura, Rohan Sachdeva"
__email__ = "liviam.moura@gmail.com, rohansach@berkeley.edu"
__status__ = "Development"

import sys
from fixame.fixame_control import Controller
from fixame.fixame_all_menus import parse_args
from fixame.fixame_logging import logger, create_fixame_logging_folders, fixame_logging


if sys.version_info[0] < 3:
    print(
        """
    ***********************************************************
    *** You're using Python 2, but FixAME needs Python 3 ;D ***
    ***********************************************************"""
    )
    sys.exit(1)

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    full_path = create_fixame_logging_folders(args.output_dir, create_folders=["tmp", "new_fastas"])
    logger = fixame_logging(logger, full_path)
    args.full_path = full_path

    control = Controller()
    control.parseArguments(args)
