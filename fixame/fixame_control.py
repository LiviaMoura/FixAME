import argparse
import os
import logging
import sys
import fixame
import fixame.fixame_curation

__author__ = "Livia Moura"
__copyright__ = "Copyright 2019"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"

class Controller():
    def __init__(self):
        self.logger = logging.getLogger()

    def curation_operation(self, **kwargs):
        logging.debug("Start teste")
        fixame.fixame_curation.fixame_curation_validate(**kwargs)
        logging.debug("!!! Finished teste!!!")

    #def error_finder_operation(self,**kwargs):

    #def merge_operation(self,**kwargs):
    
    def parseArguments(self, args):
        ''' Organize the arguments'''

        # Load the workDirectory
        #wd_loc = str(os.path.abspath(args.work_directory))
        #wd = WorkDirectory(wd_loc)

        # Set up the logger
        #self.setup_logger(wd.get_loc('log'))
        #logging.debug(str(args))

        # Call the appropriate workflow
        if args.operation == "curation":
            self.curation_operation(**vars(args))
        if args.operation == "error_finder":
            self.error_finder_operation(**vars(args))
        if args.operation == "merge":
            self.merge_operation(**vars(args))
        #if args.operation == "find_circular":
        #    self.find_circular_operation(**vars(args))

    def loadDefaultArgs(self):
        pass
