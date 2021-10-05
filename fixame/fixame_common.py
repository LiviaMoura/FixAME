import os
import sys

from fixame.fixame_logging import logger


__author__ = "Livia Moura"
__copyright__ = "Copyright 2021"
__maintainer__ = "Livia Moura, Rohan Sachdeva"
__email__ = "liviam.moura@gmail.com, rohansach@berkeley.edu"
__status__ = "Development"


def common_validate(**kwargs):
    """
    *common args*
    fasta:          fasta file for genome|metagenome [.fasta|.fa|.fna]
    output_dir:     output directory where it'll created a fixame_[date] folder
    bins:           folder cointaining bins [.fasta|.fa|.fna]
    #r12:            Interlaced SYNCED forward and reverse paired-end reads
    r1:             Forward paired-end reads
    r2:             Reverse paired-end reads
    """

    if kwargs.get("fasta") and kwargs.get("bins"):
        logger.info("Checking chosen mode [genome/metagenome|bin]")
        logger.error(
            "If you want to curate --bins, do not set --fasta, and vice-versa"
        )
        sys.exit()

    # if not os.path.exists(kwargs.get("output_dir")):
    #     logger.info("Checking output folder")
    #     logger.error(
    #         "The given path {} does not exist!".format(kwargs.get("output_dir"))
    #     )
    #     sys.exit()
    
    #     if not os.path.exists(loc):
    #             os.makedirs(loc)

    if (kwargs.get("minid") < 0.76) or (kwargs.get("minid") > 1.00):
        logger.info("Checking minimum identify for the first alignment")
        logger.error("Please, verify if -minid is >= 0.76 and <= 1.00 ")
        sys.exit()

    if (kwargs.get("num_mismatch") < 0) or (kwargs.get("num_mismatch") > 5):
        logger.info("Checking number of mismatches allowed for the filtering reads")
        logger.error("-num_mismatch must be 0>=x>=5")
        sys.exit()

    method = 0
    if kwargs.get("fasta"):
        if os.path.isfile(kwargs.get("fasta")):
            logger.info("Checking if fasta file exists")
        else:
            logger.error("There is no fasta file {}".format(kwargs.get("fasta")))
            sys.exit()

        if not os.path.splitext(kwargs.get("fasta"))[1][1:].strip().lower() in (
            "fasta",
            "fa",
            "fna",
        ):
            logger.error("The fasta file should be [.fasta|.fa|.fna]")
            sys.exit()

    if kwargs.get("bins"):
        if not os.path.exists(kwargs.get("bins")):
            logger.info("Checking bins folder")
            logger.error(
                "The given path for the bins {} does not exist!".format(
                    kwargs.get("bins")
                )
            )
            sys.exit()
        method = 1

    return method


def script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
