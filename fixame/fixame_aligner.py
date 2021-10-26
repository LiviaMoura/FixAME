import os, sys
import subprocess
from fixame.fixame_logging import logger

__author__ = "Livia Moura"
__copyright__ = "Copyright 2021"
__maintainer__ = "Livia Moura, Rohan Sachdeva"
__email__ = "liviam.moura@gmail.com, rohansach@berkeley.edu"
__status__ = "Development"


def aligner(
    output_dir, thread, minid, fasta_path, r1, r2, r12, bam_out, semi=False, last=False
):
    subprocess.run(
        [
            "bbmap.sh",
            "ref=" + str(fasta_path),
            "path=" + output_dir + "/new_fastas",
            "overwrite=t",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if not semi:
        if not r12:
            result = subprocess.run(
                [
                    os.path.join("bbmap.sh"),
                    "ref=" + str(fasta_path),
                    "path=" + output_dir + "/new_fastas",
                    "in1=" + r1,
                    "in2=" + r2,
                    "outm=" + output_dir + "/tmp/" + bam_out + ".bam",
                    "threads=" + str(thread),
                    "minid=" + str(minid),
                    "ambiguous=random",
                    "overwrite=t",
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
            if 'Duplicated sequence' in str(result.stderr):
                logger.error('There are duplicate fastas names in your fasta/bins, please fix it before running Fixame')
                sys.exit()

        else:
            result = subprocess.run(
                [
                    os.path.join("bbmap.sh"),
                    "ref=" + str(fasta_path),
                    "path=" + output_dir + "/new_fastas",
                    "in=" + r12,
                    "outm=" + output_dir + "/tmp/" + bam_out + ".bam",
                    "threads=" + str(thread),
                    "minid=" + str(minid),
                    "ambiguous=random",
                    "overwrite=t",
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
            if 'Duplicated sequence' in str(result.stderr):
                logger.error('There are duplicate fastas names in your fasta/bins, please fix it before running Fixame')
                sys.exit()
    else:
        result = subprocess.run(
            [
                os.path.join("bbmap.sh"),
                "ref=" + str(fasta_path),
                "path=" + output_dir + "/new_fastas",
                "in1=" + r1,
                "in2=" + r2,
                "outm=" + output_dir + "/tmp/" + bam_out + ".bam",
                "threads=" + str(thread),
                "ambiguous=random",
                "killbadpairs=t",
                "semiperfectmode=t",
                "overwrite=t",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        if 'Duplicated sequence' in str(result.stderr):
            logger.error('There are duplicate fastas names in your fasta/bins, \nplease fix it before running Fixame')
            sys.exit()

    subprocess.run(
        [
            "samtools",
            "sort",
            output_dir + "/tmp/" + bam_out + ".bam",
            "-o",
            output_dir + "/tmp/" + bam_out + "_sorted.bam",
            "--threads",
            str(thread),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        ["samtools", "index", output_dir + "/tmp/" + bam_out + "_sorted.bam"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    os.remove(output_dir + "/tmp/" + bam_out + ".bam")
