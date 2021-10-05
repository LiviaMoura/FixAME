from fixame.fixame_curation import main as curation_main
from fixame.fixame_error_finder import main as error_finder_main

__author__ = "Livia Moura"
__copyright__ = "Copyright 2021"
__maintainer__ = "Livia Moura, Rohan Sachdeva"
__email__ = "liviam.moura@gmail.com, rohansach@berkeley.edu"
__status__ = "Development"


class Controller:
    def __init__(self):
        pass

    def curation_operation(self, **kwargs):
        curation_main(**kwargs)

    def error_finder_operation(self, **kwargs):
        error_finder_main(**kwargs)

    def merge_operation(self, **kwargs):
        pass

    def parseArguments(self, args):
        """Organize the arguments"""

        # Call the appropriate workflow
        if args.operation == "curation":
            self.curation_operation(**vars(args))
        if args.operation == "error_finder":
            self.error_finder_operation(**vars(args))
        # if args.operation == "merge":
        #    self.merge_operation(**vars(args))
        # if args.operation == "find_circular":
        #    self.find_circular_operation(**vars(args))

    def loadDefaultArgs(self):
        pass
