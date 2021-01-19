##!/usr/bin/env python3

#***************************
#***************************
#*** Nosso ponto inicial ***
#***************************
#***************************

__author__ = "Livia Moura"
__copyright__ = "Copyright 2021"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"

import argparse
import sys
import os
from fixame.fixame_control import Controller

if sys.version_info[0] < 3:
    print('''
    ***********************************************************
    *** You're using Python 2, but FixAME needs Python 3 ;D ***
    ***********************************************************''')
    sys.exit(1)

import fixame.fixame_all_menus

if __name__ == '__main__':
    args = fixame.fixame_all_menus.parse_args(sys.argv[1:])
    print(args)
    # do what you came here to do :)
    control = Controller()
    control.parseArguments(args)
