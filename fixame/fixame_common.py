import argparse
import os
import logging
import sys
import fixame
from xopen import xopen

__author__ = "Livia Moura"
__copyright__ = "Copyright 2019"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"


def common_validate(**kwargs):
    '''
    *common args*
    fasta:          fasta file for genome|metagenome [.fasta|.fa|.fna]
    output_dir:     output directory where it'll created a fixame_[date] folder
    bins:           folder cointaining bins [.fasta|.fa|.fna]
    #r12:            Interlaced SYNCED forward and reverse paired-end reads
    r1:             Forward paired-end reads
    r2:             Reverse paired-end reads
    '''
    
    if kwargs.get('fasta') and kwargs.get('bins'):
        logging.info('Checking chosen mode [genome/metagenome|bin]')
        logging.error("If you want to curate --bins, do not set --fasta, and vice-versa")
        sys.exit()
        
    if not os.path.exists(kwargs.get('output_dir')):
        logging.info('Checking output folder')
        logging.error("The given path {} does not exist!".format(kwargs.get('output_dir')))
        sys.exit()
    
    method = 0
    if kwargs.get('fasta'):
        if os.path.isfile(kwargs.get('fasta')):
            logging.info('Checking if fasta file exists')
        else:
            logging.error('There is no fasta file {}'.format(kwargs.get('fasta')))
            sys.exit()
        
        if not os.path.splitext(kwargs.get('fasta'))[1][1:].strip().lower() in ("fasta","fa","fna"): 
            logging.error('The fasta file should be [.fasta|.fa|.fna]')
            sys.exit()


    if kwargs.get('bins'):
        if not os.path.exists(kwargs.get('bins')):
            logging.info('Checking bins folder')
            logging.error("The given path for the bins {} does not exist!".format(kwargs.get('bins')))
            sys.exit()
        method = 1    
    
    return method

def script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
