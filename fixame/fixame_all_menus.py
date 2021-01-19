__author__ = "Livia Moura"
__copyright__ = "Copyright 2020"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"

import argparse,textwrap
import sys
import os
import fixame

import argparse
import textwrap
import sys


class CustomFormatter(argparse.HelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            # change to
            #    -s, --long ARGS
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    #parts.append('%s %s' % (option_string, args_string))
                    parts.append('%s' % option_string)
                parts[-1] += ' %s' % args_string
            return ', '.join(parts)


def printHelp(version):
    print('')
    print('--- FixAME - Love to Fix ---')
    print('        v. {}'.format(version))
    print('''\

Workflow
  curation         Locate and try to fix potentials locals assembly errors.

Single process
  error_finder     Return the positional of assembly errors
  merge            Locate overlaps between contigs edges and connect them
  
Example: fixame [curation|error_finder|merge|find_circular] --help
''')


def cf(prog): return CustomFormatter(prog)


parser = argparse.ArgumentParser(
    prog='fixame',
    epilog="See '<command> --help'",
    formatter_class=cf
)

system_parser = argparse.ArgumentParser(add_help=False)
system_parser_opt = system_parser.add_argument_group('System Options')
system_parser_opt.add_argument(
    '-t', '--threads', help='Number of threads [6]', default=6, type=int)
system_parser_opt.add_argument("--keep", action="store_true", default=False,
                               help="Keep intermediate files [False] (It uses a lot of disk space)")
# system_parser_opt.add_argument(
#    "-h", "--help", action="help", help="show this help message and exit")

req_parser = argparse.ArgumentParser(add_help=False)
req_parser_opt = req_parser.add_argument_group('Required Inputs')
req_parser_opt.add_argument(
    '-f', '--fasta', help='fasta file for genome|metagenome [.fasta|.fa|.fna]')
req_parser_opt.add_argument(
    '-b', '--bins', help='folder cointaining bins [.fasta|.fa|.fna]')
req_parser_opt.add_argument(
    '-r12', help='Interlaced SYNCED forward and reverse paired-end reads')
req_parser_opt.add_argument('-r1', help='Forward paired-end reads')
req_parser_opt.add_argument('-r2', help='Reverse paired-end reads')
req_parser_opt.add_argument('-o', '--output_dir', help="output directory")
req_parser_opt.add_argument(
    "-l", "--min_assembly_len", action="store", help="Minimum fasta3 length")

curation_parser = argparse.ArgumentParser(add_help=False)
curation_parser_opt = curation_parser.add_argument_group('Curation')
curation_parser_opt.add_argument(
    "-x", "--xtimes", help="Number of alignments during the curation [10]", default=1, type=int)
curation_parser_opt.add_argument(
    "-minid", help="Minumum identity for the first alignment [0.9]", default=0.9, type=float)
curation_parser_opt.add_argument(
    "-cov", "--fasta_cov", help="Local errors will be called on regions with coverage >= [5]", default=5, type=int)
curation_parser_opt.add_argument(
    "--dp_cov", help="Number of reads needed to extend the gaps/curate [1]", default=1, type=int)
curation_parser_opt.add_argument(
    "-nm", "--num_mismatch", help="Number of mismatches allowed to filter out the initial reads [2]", default=2, type=int)
curation_parser_opt.add_argument("--ext_multifasta", help="Execute the merge between curated contigs [True]",
                                 choices=['True', 'False'], default="True")

# Creating the subparse
subparsers = parser.add_subparsers(dest='act')

curation_call = subparsers.add_parser(
    'curation', help='curation option', parents=[req_parser, curation_parser, system_parser], formatter_class=cf)

merge_call = subparsers.add_parser(
    'merge', help='merge command', parents=[req_parser])


# Calling menus
args = parser.parse_args()


if (len(sys.argv[1:]) == 0 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
    version = '1.0.0'
    printHelp(version)
    sys.exit(0)
else:
    if args.act == 'curation':
        pass
    elif args.act == 'error_finder':
        pass
    elif args.act == 'merge':
        pass

print(args)



  
