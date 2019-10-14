__author__ = "Livia Moura"
__copyright__ = "Copyright 2019"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"

import argparse,textwrap
import sys
import os
import fixame

class CustomHelpFormatter(argparse.HelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string

fmt = lambda prog: CustomHelpFormatter(prog)

def version():
    versionFile = open(os.path.join(fixame.__path__[0], 'VERSION'))
    return versionFile.read().strip()

VERSION = version()

def printHelp():
    print('')
    print('               -|╞ FixAME - Love to Fix ╡|-')
    print('                       v. '+VERSION)
    print('''\


 ► Workflow      
 curation        Locate and try to fix potentials locals assembly errors. 
                 If desired, it'll merge the curated contigs into larger contigs or scaffold them.                 

 ► Simple process 
 error_finder    Only locate the position of the errors from the 
                    original file and report them 
 merge           Locate overlaps between contigs edges and connect them and single reads that connect contigs and join them creating a longer
                    contig or a scaffold

 Example: fixame [curation|error_finder|merge] -h 
 
 Further details can be found at: Manuscript
 If you use it, please cite us :) 

    ''')

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=fmt)
    subparsers = parser.add_subparsers(help='Desired operation',dest='operation')

    #SYSTEM submenu
    system_parser = argparse.ArgumentParser(add_help=False)
    sysflags = system_parser.add_argument_group('SYSTEM')
    sysflags.add_argument('-t','--threads',help='Number of threads [6]',default=6,type=int)
    sysflags.add_argument("--keep", action="store_true", default=False, help="Keep intermediate files [False] (It uses a lot of disk space)")
    sysflags.add_argument("-h", "--help", action="help", help="show this help message and exit")

    #REQUERIMENTS submenu
    require_parser = argparse.ArgumentParser(add_help=False)
    reqflags = require_parser.add_argument_group('REQUIRED')
    reqflags.add_argument('-f','--fasta',help='fasta file for genome|metagenome [.fasta|.fa|.fna]')
    reqflags.add_argument('-b','--bins',help='folder cointaining bins [.fasta|.fa|.fna]')        
    reqflags.add_argument('-r12',help='Interlaced SYNCED forward and reverse paired-end reads')
    reqflags.add_argument('-r1',metavar="R1", help='Forward paired-end reads')
    reqflags.add_argument('-r2',metavar="R2", help='Reverse paired-end reads')
    reqflags.add_argument('-o','--output_dir', help="output directory")
    
    #CURATION submenu
    
 #   curation_parser = subparsers.add_parser("curation",formatter_class=fmt,parents = [curation_parent,require_parser,system_parser], add_help=False)

    curation_parent = argparse.ArgumentParser(add_help=False)
    curflags = curation_parent.add_argument_group('CURATION PARAMETERS')
    curflags.add_argument("-x","--xtimes", help= "Number of alignments during the curation [10]",
                            default=10,metavar="INT", type = int)
    curflags.add_argument("-minid",metavar="FLOAT", help="Minumum identity for the first alignment [0.9]", 
                            default = 0.9, type = float)
    curflags.add_argument("-cov","--fasta_cov", metavar="INT", help="Local errors will be called on regions with coverage >= [5]",
                            default = 5, type = int)
    curflags.add_argument("--dp_cov", metavar="INT", help="Number of reads needed to extend the gaps/curate [1]", 
                            default = 1, type = int)
    curflags.add_argument("-nm","--num_mismatch", metavar="INT", help="Number of mismatches allowed to filter out the initial reads [2]", 
                            default = 2, type = int)
    curflags.add_argument("--ext_multifasta", help="Execute the merge between curated contigs [True]",
                            choices=['True','False'], default = "True")

    #CHAMANDO OS MENUS
    curation_parser = subparsers.add_parser("curation", description=textwrap.dedent('''\
            
                               -o| CURATION workflow |o-

            The complete workflow to curate your genome, bins or metagenome
            Here, It'll look for assembly local errors and try to fix them. 
            
            If possible, it'll merge contigs edges and create scaffolds with
            the curated version (default).
            If bins, the merge will only considers the contigs from that BIN.


     '''),formatter_class=argparse.RawDescriptionHelpFormatter, parents=[require_parser,curation_parent,system_parser], add_help=False,
    epilog=textwrap.dedent('''\
            Usage: fixame curation -f {fasta_file} -r1 {R1.fastq} -r2 {R2.fastq} -o {path} 
                   fixame curation -b {bin_folder} -r12 {R12.fastq} -o {path}
   
                   '''))
    
    error_finder_parser = subparsers.add_parser("error_finder", description=textwrap.dedent('''\
            
                       -o| Error_finder script |o-

            Reports the assembly local errors positions


     '''),formatter_class=argparse.RawDescriptionHelpFormatter, parents=[require_parser,system_parser], add_help=False,
    epilog=textwrap.dedent('''\
            Usage: fixame error_finder -f {fasta_file} -r1 {R1.fastq} -r2 {R2.fastq} -o {path} 
                   fixame error_finder -b {bin_folder} -r12 {R12.fastq} -o {path}
   
                   ''')) 
     
    
    merge_parser = subparsers.add_parser("merge", description=textwrap.dedent('''\
            
                                    -o| Merge script |o-

            If possible, it'll merge contigs edges and create scaffolds based on aligned reads. 
            If bins, the merge will only considers the contigs from that BIN

     '''),formatter_class=argparse.RawDescriptionHelpFormatter, parents=[require_parser,system_parser], add_help=False,
     epilog=textwrap.dedent('''\
            Usage: fixame merge -f {fasta_file} -r1 {R1.fastq} -r2 {R2.fastq} -o {path}
                   fixame merge -b {bin_folder} -r12 {R12.fastq} -o {path} 
   
                   '''))

    # Handle the situation where the user wants the raw help
    #args = None
    #if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv == '--help'):
    if (len(args) == 0 or args[0] == '-h' or args[0] == '--help'):
        printHelp()
        sys.exit(0)
    else:
        return parser.parse_args(args)


  
