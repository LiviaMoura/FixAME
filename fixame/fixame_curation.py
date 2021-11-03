from logging import error
import os
import sys
import re
import subprocess
import shutil
import pysam as ps
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2

from fixame.fixame_aligner import aligner
from fixame.fixame_common import common_validate
from fixame.fixame_logging import logger
from fixame.fixame_error_finder import (
    calculate_reference_lengths,
    parse_map,
    check_local_assembly_errors_parallel,
    organizing_found_errors,
    check_direct_features_parallel,
)

from functools import partial
from concurrent.futures import ProcessPoolExecutor
from xopen import xopen


__author__ = "Livia Moura"
__copyright__ = "Copyright 2021"
__maintainer__ = "Livia Moura, Rohan Sachdeva"
__email__ = "liviam.moura@gmail.com, rohansach@berkeley.edu"
__status__ = "Development"


def main(**kwargs):
    """
    *Fixame curation args:
    xtimes:         Number of alignments during the curation [10]
    minid:          Minumum identity for the first alignment [0.9]
    fasta_cov:      Local errors will be called on regions with coverage >= [5]
    dp_cov:         Number of reads needed to extend the gaps/curate [1]
    num_mismatch:   Number of mismatches allowed to filter out the initial reads [2]
    ext_multifasta: Execute the merge between curated contigs [True]

    *common args*
    fasta:          fasta file for genome|metagenome [.fasta|.fa|.fna]
    output_dir:     output directory where it'll created a fixame_[date] folder
    bins:           folder cointaining bins [.fasta|.fa|.fna]
    r12:            Interlaced SYNCED forward and reverse paired-end reads
    r1:             Forward paired-end reads
    r2:             Reverse paired-end reads

    *Returns:
    """
    mydir = kwargs.get("full_path")

    logger.info("Starting Curation")

    if kwargs.get("min_ctg_len") < 800:
        logger.info("Checking minimum contig length")
        logger.error("--min_ctg_len must be >= 800")
        sys.exit()
    else:
        minimum_assembly_length = kwargs.get("min_ctg_len")

    if kwargs.get("r12"):
        read12_in = os.path.realpath(os.path.expanduser(kwargs.get("r12")))
        read1_in = ""
        read2_in = ""
    else:
        read1_in = os.path.realpath(os.path.expanduser(kwargs.get("r1")))
        read2_in = os.path.realpath(os.path.expanduser(kwargs.get("r2")))
        read12_in = ""

    av_readlen = temp_average_read(kwargs.get("r1"))

    # Checking the pipeline - genome/metagenoms vs Bins
    method = common_validate(**kwargs)

    if method == 0:
        fasta_in = os.path.realpath(os.path.expanduser(kwargs.get("fasta")))
        name_fasta = os.path.splitext(os.path.basename(fasta_in))[0]

        logger.info("\n --- Analysing the file {} ---\n".format(name_fasta))

        try:
            logger.info(
                "Checking overlaping at N regions on {} and fix them".format(fasta_in)
            )
            os.mkdir(os.path.join(mydir, "fixing_log"))

            closed_N = open(
                os.path.join(mydir, "fixing_log", "fixame_initial_Ns_closed.txt"),
                "w+",
            )
            check_overlap(
                mydir, fasta_in, av_readlen, kwargs.get("threads"), True, fixed=closed_N
            )
            logger.info(
                "A new reference fasta {} was created".format(
                    mydir + "/new_fastas/" + name_fasta + "_renewed.fasta"
                )
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        try:
            logger.info("Mapping reads against the new reference")
            aligner(
                mydir,
                kwargs.get("threads"),
                kwargs.get("minid"),
                mydir + "/new_fastas/" + name_fasta + "_renewed.fasta",
                r1=read1_in,
                r2=read2_in,
                r12=read12_in,
                bam_out=name_fasta + "_renewed",
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        fasta_cov, num_mm = kwargs.get("fasta_cov"), kwargs.get("num_mismatch")

        # Filtering the fastq - Make curation process faster
        try:
            logger.info("Filtering out reads that didn't map the new ref")
            filtering_bam(
                mydir,
                kwargs.get("threads"),
                num_mm,
                mydir + "/tmp/" + name_fasta + "_renewed",
                read1_in,
                read2_in,
                read12_in,
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        try:
            logger.info("Generating some metrics to keep running")
            reference_to_length = calculate_reference_lengths(
                mydir + "/new_fastas/" + name_fasta + "_renewed.fasta",
                minimum_assembly_length,
            )

            (
                bam_dict,
                reference_read_lengths,
                average_template_length,
                average_read_length,
                average_gap_length,
                template_length_min,
                template_length_max,
                average_gap_std,
            ) = parse_map(
                mydir + "/tmp/" + name_fasta + "_renewed_sorted.bam",
                num_mm,
                kwargs.get("threads"),
                minimum_assembly_length,
                reference_to_length,
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        features_list_dict = check_direct_features_parallel(
            mydir + "/new_fastas/" + name_fasta + "_renewed.fasta",
            kwargs.get("threads"),
        )

        try:
            logger.info("Trying to find regions with local assembly errors")
            (
                reference_to_error_regions,
                coverage_dict,
                reference_to_high_mismatch_positions,
            ) = check_local_assembly_errors_parallel(
                reference_to_length.keys(),
                kwargs.get("threads"),
                reference_read_lengths,
                reference_to_length,
                fasta_cov,
                bam_dict,
                num_mm,
                template_length_max,
            )

        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        try:
            logger.info("Selecting the errors regions")
            organized_errors = organizing_found_errors(
                average_read_length, reference_to_error_regions
            )

            fasta_len, error_df = build_N(
                mydir,
                mydir + "/new_fastas/" + name_fasta + "_renewed.fasta",
                average_read_length,
                organized_errors,
            )

        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        logger.info("\nStarting to fix sample {}\n".format(name_fasta))

        for count, r in enumerate(range(kwargs.get("xtimes")), 1):
            fixed = open(
                os.path.join(mydir, "fixing_log", "fixame_loop_" + str(count) + ".txt"),
                "w+",
            )
            try:
                logger.info("Loop {} from {}".format(count, kwargs.get("xtimes")))
                var_cal_fix(
                    mydir,
                    count,
                    fixed,
                    kwargs.get("threads"),
                    kwargs.get("dp_cov"),
                    av_readlen,
                    error_df,
                )
            except Exception as e:
                logger.exception("Something went wrong")
                raise e

            fixed.close()

        logger.info("Errors fixing complete")
        logger.info(
            "You can find the fixing log at {}".format(
                os.path.join(mydir, "fixing_log")
            )
        )
        os.mkdir(os.path.join(mydir, "fixame_results"))
        logger.info("Polishing the sequences...")
        remove_N(
            mydir,
            name_fasta,
            os.path.join(mydir, "tmp", "v_" + str(kwargs.get("xtimes")) + ".fasta"),
            organized_errors,
            fasta_len,
            av_readlen,
            average_gap_length,
            average_gap_std,
            kwargs.get("threads"),
        )

        error_df["sample_name"] = name_fasta

        final_output(mydir, error_df, features_list_dict)

        if kwargs.get("keep") is False:
            try:
                logger.info("Removing temporary files")
                shutil.rmtree(os.path.join(mydir, "tmp"))
            except:
                logger.info("It wasn't possible to remove the /tmp folder")
        logger.info("\n\nFixame proccess done!\n")

    else:  # BINS MODE
        fasta_array = []
        name_sample = "bins"
        contigs_bins = open(os.path.join(mydir, "tmp", "bin_contigs.txt"), "w+")
        # contigs_bins = defaultdict(list)
        merged = open(os.path.join(mydir, "tmp", "bins.fasta"), "w+")
        for sample in os.listdir(
            os.path.realpath(os.path.expanduser(kwargs.get("bins")))
        ):
            name = sample.split(".")[0]
            if (
                sample.lower().endswith(".fasta")
                or sample.lower().endswith(".fa")
                or sample.lower().endswith(".fna")
            ):
                fasta_array.append(sample)
                for seq_record in SeqIO.parse(
                    os.path.join(kwargs.get("bins"), sample), "fasta"
                ):
                    contigs_bins.write("{}\t{}\n".format(name, seq_record.id))
                with open(os.path.join(kwargs.get("bins"), sample), "r") as readfile:
                    shutil.copyfileobj(readfile, merged)

        contigs_bins.close()
        merged.close()

        fasta_in = os.path.join(mydir, "tmp", "bins.fasta")

        logger.info(
            "\n --- Analysing a metagenome sample with {} bins ---\n".format(
                len(fasta_array)
            )
        )
        try:
            logger.info(
                "Checking overlaping at N regions on {} and fix them".format(fasta_in)
            )
            os.mkdir(os.path.join(mydir, "fixing_log"))
            closed_N = open(
                os.path.join(mydir, "fixing_log", "fixame_initial_Ns_closed.txt"),
                "w+",
            )
            check_overlap(
                mydir, fasta_in, av_readlen, kwargs.get("threads"), True, fixed=closed_N
            )
            logger.info(
                "A new reference fasta {} was created".format(
                    mydir + "/new_fastas/" + name_sample + "_renewed.fasta"
                )
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        try:
            logger.info("Mapping reads against the new reference")
            aligner(
                mydir,
                kwargs.get("threads"),
                kwargs.get("minid"),
                mydir + "/new_fastas/" + name_sample + "_renewed.fasta",
                r1=read1_in,
                r2=read2_in,
                r12=read12_in,
                bam_out=name_sample + "_renewed",
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        fasta_cov, num_mm = kwargs.get("fasta_cov"), kwargs.get("num_mismatch")

        # Filtering the fastq - Make curation process faster
        try:
            logger.info("Filtering out reads that didn't map the new ref")
            filtering_bam(
                mydir,
                kwargs.get("threads"),
                num_mm,
                mydir + "/tmp/" + name_sample + "_renewed",
                read1_in,
                read2_in,
                read12_in,
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        try:
            logger.info("Generating some metrics to keep running")
            reference_to_length = calculate_reference_lengths(
                mydir + "/new_fastas/" + name_sample + "_renewed.fasta",
                minimum_assembly_length,
            )
            (
                bam_dict,
                reference_read_lengths,
                average_template_length,
                average_read_length,
                average_gap_length,
                template_length_min,
                template_length_max,
                average_gap_std,
            ) = parse_map(
                mydir + "/tmp/" + name_sample + "_renewed_sorted.bam",
                num_mm,
                kwargs.get("threads"),
                minimum_assembly_length,
                reference_to_length,
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        features_list_dict = check_direct_features_parallel(
            mydir + "/new_fastas/" + name_sample + "_renewed.fasta",
            kwargs.get("threads"),
        )

        try:
            logger.info("Trying to find regions with local assembly errors")
            (
                reference_to_error_regions,
                coverage_dict,
                reference_to_high_mismatch_positions,
            ) = check_local_assembly_errors_parallel(
                reference_to_length.keys(),
                kwargs.get("threads"),
                reference_read_lengths,
                reference_to_length,
                fasta_cov,
                bam_dict,
                num_mm,
                template_length_max,
            )
        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        try:
            logger.info("Selecting the errors regions")
            organized_errors = organizing_found_errors(
                average_read_length, reference_to_error_regions
            )

            fasta_len, error_df = build_N(
                mydir,
                mydir + "/new_fastas/" + name_sample + "_renewed.fasta",
                average_read_length,
                organized_errors,
            )

        except Exception as e:
            logger.exception("Something went wrong")
            raise e

        logger.info("Starting to fix all bins\n")

        for count, r in enumerate(range(kwargs.get("xtimes")), 1):
            fixed = open(
                os.path.join(mydir, "fixing_log", "fixame_loop_" + str(count) + ".txt"),
                "w+",
            )
            try:
                logger.info("Loop {} from {}".format(count, kwargs.get("xtimes")))
                var_cal_fix(
                    mydir,
                    count,
                    fixed,
                    kwargs.get("threads"),
                    kwargs.get("dp_cov"),
                    av_readlen,
                    error_df,
                )
            except Exception as e:
                logger.exception("Something went wrong")
                raise e

            fixed.close()

        logger.info("Errors fixing complete")
        logger.info(
            "You can find the fixing log at {}".format(
                os.path.join(mydir, "fixing_log")
            )
        )
        os.mkdir(os.path.join(mydir, "fixame_results"))
        logger.info("Polishing the sequences...")
        remove_N(
            mydir,
            name_sample,
            os.path.join(mydir, "tmp", "v_" + str(kwargs.get("xtimes")) + ".fasta"),
            organized_errors,
            fasta_len,
            av_readlen,
            average_gap_length,
            average_gap_std,
            kwargs.get("threads"),
        )

        ## Spliting the bins
        df = pd.read_table(
            os.path.join(mydir, "tmp", "bin_contigs.txt"),
            names=["sample_name", "contig"],
        )
        error_df = pd.merge(error_df, df, on=["contig"], how="left")
        df = df.groupby("sample_name").agg({"contig": lambda x: list(x)}).reset_index()

        final_output(mydir, error_df, features_list_dict)

        for index, (sample_name, fasta_id) in df.iterrows():
            fasta_ids = "\n".join(fasta_id)
            tmp_id = open(os.path.join(mydir, "tmp", "tmp_fasta"), "w")
            tmp_id.write(fasta_ids)
            tmp_id.close()
            cmd = "filterbyname.sh in={} out={} names={} include=t".format(
                os.path.join(mydir, "tmp", "bins_unordered.fasta"),
                os.path.join(mydir, "fixame_results", sample_name + "_fixame.fasta"),
                os.path.join(mydir, "tmp", "tmp_fasta"),
            )
            subprocess.run(
                cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )

            cmd = "filterbyname.sh in={} out={} names={} include=t".format(
                os.path.join(mydir, "new_fastas", "bins_renewed.fasta"),
                os.path.join(mydir, "new_fastas", sample_name + "_renewed.fasta"),
                os.path.join(mydir, "tmp", "tmp_fasta"),
            )
            subprocess.run(
                cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )

        os.remove(os.path.join(mydir, "tmp", "bins_unordered.fasta"))
        os.remove(os.path.join(mydir, "new_fastas", "bins_renewed.fasta"))
        os.remove(os.path.join(mydir, "tmp", "tmp_fasta"))

        if not kwargs.get("keep"):
            try:
                logger.info("Removing temporary files")
                shutil.rmtree(os.path.join(mydir, "tmp"))
            except:
                logger.info("It wasn't possible to remove the /tmp folder")
        logger.info("\n\nFixame proccess done!\n")


def temp_average_read(r1_fastq):
    """Initial fast avg read length - It'll be recalculated with more precision later"""
    sum_read_len = 0
    for num, record in enumerate(SeqIO.parse(xopen(r1_fastq), "fastq")):
        sum_read_len += len(record.seq)
        if num == 9999:
            av_read_len = int(sum_read_len / 10000)
            break
    return av_read_len


def check_overlap(
    output_dir,
    fasta,
    av_readlen,
    thread,
    user_file=False,
    error_df=pd.DataFrame(),
    fixed="",
    count="",
):
    """Check overlap from the border of N's regions"""
    if user_file:
        name_fasta = os.path.splitext(os.path.basename(fasta))[0]
        new_fasta_temp = open(
            os.path.join(output_dir, "new_fastas", name_fasta + "_renewed.fasta"), "w"
        )
    else:
        if os.path.exists(
            os.path.join(output_dir, "tmp", "target")
        ):  ### remove target file if it exists
            # copied = open(os.path.join(output_dir, "tmp", "target_"+str(count)), "w+")
            shutil.copyfile(
                os.path.join(output_dir, "tmp", "target"),
                os.path.join(output_dir, "tmp", "target_" + str(count)),
            )
            # copied.close()
            os.remove(os.path.join(output_dir, "tmp", "target"))

        new_fasta_temp = open(
            os.path.join(output_dir, "tmp", "v_" + str(count) + ".fasta"), "w"
        )

    items = [seq_record for seq_record in SeqIO.parse(fasta, "fasta")]

    with ProcessPoolExecutor(thread) as executor:
        check_npos_result = executor.map(partial(check_npos, av_readlen), items)

    npos_list = list()
    ntarget_list = list()

    for npos_unit in check_npos_result:
        npos_list.append(npos_unit[0])
        ntarget_list.append(npos_unit[1])

    N_pos = {}
    N_target = {}
    for i in npos_list:
        N_pos.update(i)
    for i in ntarget_list:
        N_target.update(i)

    if not user_file:
        for k, v in N_pos.items():
            new_value = v[1:-1]
            N_pos[k] = new_value

    for item in items:
        contig_name = item.id
        new_target = list()
        control_index = list()
        new_subfasta = ""

        if not N_pos.get(item.id):
            new_subfasta = item.seq

        if contig_name in N_pos.keys():
            temp_dif_pos = 0

            if not error_df.empty:
                for (start, end, number) in N_pos.get(contig_name):
                    idx = error_df.loc[
                        (error_df["contig"] == contig_name)
                        & ((start - error_df["N_build_start"]) <= av_readlen)
                        & ((start - error_df["N_build_start"]) >= -av_readlen)
                        & ((error_df["N_build_end"] - end) <= av_readlen)
                        & ((error_df["N_build_end"] - end) >= -av_readlen)
                    ].index.tolist()
                    if idx:
                        error_df.at[idx[0], "N_build_start"] = start
                        error_df.at[idx[0], "N_build_end"] = end

            for j, (start, end, number) in enumerate(N_pos.get(contig_name)):

                if j == 0:  ### first time doesnt have index modification
                    temp_left = item[
                        0:start
                    ].seq  ### start = N pos; but here would be Last base pos
                    temp_right = item[
                        end + 1 : -1
                    ].seq  ### end = N pos; but here would be first base pos
                    comp_left = temp_left[-(av_readlen * 2) :]
                    comp_right = temp_right[0:12]

                    ## Self remember alignment list 0 to N
                    if comp_right in comp_left:
                        alignments = pairwise2.align.localms(
                            comp_left, comp_right, 2, -1, -0.5, -0.1
                        )
                        variat = (av_readlen * 2) - alignments[0][3]
                        new_subfasta = temp_left[:-variat] + temp_right
                        control_index.append(
                            tuple([j + 1, N_pos.get(contig_name)[j][2] + variat])
                        )

                        if fixed != "":
                            fixed.write(
                                "{}\t{}\t{}\tfixed\n".format(
                                    contig_name, start - variat, end + 1 + 12
                                ),
                            )

                        if not error_df.empty:
                            idx = list()
                            idx = error_df.loc[
                                (error_df["contig"] == contig_name)
                                & (error_df["status"] != "fixed")
                                & (
                                    (start - variat)
                                    > (error_df["N_build_start"] - 2 * av_readlen)
                                )
                                & ((start - variat) < (error_df["N_build_end"] - 1))
                            ].index.tolist()
                            temp_dif_pos += number + variat

                            if idx:
                                error_df.at[idx[0], "status"] = "fixed"
                                error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (
                                        error_df["order"] > error_df.at[idx[0], "order"]
                                    ),
                                    "N_build_start",
                                ] = (
                                    error_df["N_build_start"] - temp_dif_pos
                                )  # + (start - error_df['N_build_start'])
                                error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (
                                        error_df["order"] > error_df.at[idx[0], "order"]
                                    ),
                                    "N_build_end",
                                ] = (
                                    error_df["N_build_end"] - temp_dif_pos
                                )  # - (error_df['N_build_end'] - end)

                            else:
                                error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (error_df["N_build_start"] > (start - variat)),
                                    "N_build_start",
                                ] = (
                                    error_df["N_build_start"] - temp_dif_pos
                                )
                                error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (error_df["N_build_start"] > (start - variat)),
                                    "N_build_end",
                                ] = (
                                    error_df["N_build_end"] - temp_dif_pos
                                )
                        else:
                            temp_dif_pos += number + variat

                    ####### REVERSE CHECK #################
                    else:
                        comp_left = temp_left[-12:]
                        comp_right = temp_right[: (av_readlen * 2)]

                        if comp_left in comp_right:
                            alignments = pairwise2.align.localms(
                                comp_right, comp_left, 2, -1, -0.5, -0.1
                            )

                            new_subfasta = temp_left + temp_right[alignments[0][4] :]
                            control_index.append(
                                tuple(
                                    [
                                        j + 1,
                                        N_pos.get(contig_name)[j][2] + alignments[0][4],
                                    ]
                                )
                            )
                            if fixed != "":
                                fixed.write(
                                    "{}\t{}\t{}\tfixed\n".format(
                                        contig_name,
                                        (start - 12 - 1),
                                        (end + 1 + alignments[0][4]),
                                    ),
                                )

                            if not error_df.empty:
                                idx = list()
                                idx = error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (error_df["status"] != "fixed")
                                    & (
                                        (end + 1 + alignments[0][4])
                                        > (error_df["N_build_start"] + 1)
                                    )
                                    & (
                                        (end + 1 + alignments[0][4])
                                        < (error_df["N_build_end"] + 2 * av_readlen)
                                    )
                                ].index.tolist()
                                temp_dif_pos += number + alignments[0][4]

                                if idx:
                                    error_df.at[idx[0], "status"] = "fixed"
                                    error_df.loc[
                                        (error_df["contig"] == contig_name)
                                        & (
                                            error_df["order"]
                                            > error_df.at[idx[0], "order"]
                                        ),
                                        "N_build_start",
                                    ] = (
                                        error_df["N_build_start"] - temp_dif_pos
                                    )
                                    error_df.loc[
                                        (error_df["contig"] == contig_name)
                                        & (
                                            error_df["order"]
                                            > error_df.at[idx[0], "order"]
                                        ),
                                        "N_build_end",
                                    ] = (
                                        error_df["N_build_end"] - temp_dif_pos
                                    )
                                else:
                                    error_df.loc[
                                        (error_df["contig"] == contig_name)
                                        & (
                                            error_df["N_build_start"] > (start - 12 - 1)
                                        ),
                                        "N_build_start",
                                    ] = (
                                        error_df["N_build_start"] - temp_dif_pos
                                    )
                                    error_df.loc[
                                        (error_df["contig"] == contig_name)
                                        & (
                                            error_df["N_build_start"] > (start - 12 - 1)
                                        ),
                                        "N_build_end",
                                    ] = (
                                        error_df["N_build_end"] - temp_dif_pos
                                    )

                            else:
                                temp_dif_pos += number + alignments[0][4]
                        else:
                            new_subfasta = item.seq
                else:
                    temp_left = new_subfasta[
                        0 : start - temp_dif_pos
                    ]  ### start = N pos; but here would be Last base pos
                    temp_right = new_subfasta[
                        end + 1 - temp_dif_pos : -1
                    ]  ### end = N pos; but here would be first base pos
                    comp_left = temp_left[-(av_readlen * 2) :]
                    comp_right = temp_right[0:12]

                    ## Self remember alignment list 0-N
                    if comp_right in comp_left:
                        alignments = pairwise2.align.localms(
                            comp_left, comp_right, 2, -1, -0.5, -0.1
                        )

                        variat = (av_readlen * 2) - alignments[0][3]
                        new_subfasta = temp_left[:-variat] + temp_right
                        control_index.append(
                            tuple([j + 1, N_pos.get(contig_name)[j][2] + variat])
                        )

                        if fixed != "":
                            fixed.write(
                                "{}\t{}\t{}\tfixed\n".format(
                                    contig_name,
                                    (start - temp_dif_pos - variat - 1),
                                    (end + 1 - temp_dif_pos + 12),
                                ),
                            )

                        if not error_df.empty:
                            idx = list()
                            idx = error_df.loc[
                                (error_df["contig"] == contig_name)
                                & (error_df["status"] != "fixed")
                                & (
                                    (start - temp_dif_pos)
                                    > (error_df["N_build_start"] - 2 * av_readlen)
                                )
                                & (
                                    (start - temp_dif_pos)
                                    < (error_df["N_build_end"] - 1)
                                )
                            ].index.tolist()
                            temp_dif_pos += number + variat

                            if idx:
                                error_df.at[idx[0], "status"] = "fixed"
                                error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (
                                        error_df["order"] > error_df.at[idx[0], "order"]
                                    ),
                                    "N_build_start",
                                ] = error_df["N_build_start"] - (
                                    number + variat
                                )  # temp_dif_pos
                                error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (
                                        error_df["order"] > error_df.at[idx[0], "order"]
                                    ),
                                    "N_build_end",
                                ] = error_df["N_build_end"] - (
                                    number + variat
                                )  # temp_dif_pos

                            else:
                                error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (
                                        error_df["N_build_start"] > (start - variat - 1)
                                    ),
                                    "N_build_start",
                                ] = error_df["N_build_start"] - (
                                    number + variat
                                )  # temp_dif_pos | (start - temp_dif_pos - variat)
                                error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (
                                        error_df["N_build_start"] > (start - variat - 1)
                                    ),
                                    "N_build_end",
                                ] = error_df["N_build_end"] - (
                                    number + variat
                                )  # temp_dif_pos

                        else:
                            temp_dif_pos += number + variat

                    #######  REVERSE CHECK  #################
                    else:
                        comp_left = temp_left[-12:]
                        comp_right = temp_right[: (av_readlen * 2)]

                        if comp_left in comp_right:
                            alignments = pairwise2.align.localms(
                                comp_right, comp_left, 2, -1, -0.5, -0.1
                            )
                            new_subfasta = temp_left + temp_right[alignments[0][4] :]
                            control_index.append(
                                tuple(
                                    [
                                        j + 1,
                                        N_pos.get(contig_name)[j][2] + alignments[0][4],
                                    ]
                                )
                            )
                            if fixed != "":
                                fixed.write(
                                    "{}\t{}\t{}\tfixed\n".format(
                                        contig_name,
                                        (start - temp_dif_pos - 12 - 1),
                                        (end + 1 - temp_dif_pos + alignments[0][4]),
                                    ),
                                )

                            if not error_df.empty:
                                idx = list()
                                idx = error_df.loc[
                                    (error_df["contig"] == contig_name)
                                    & (error_df["status"] != "fixed")
                                    & (
                                        (end + 1 - temp_dif_pos)
                                        > (error_df["N_build_start"] + 1)
                                    )
                                    & (
                                        (end + 1 - temp_dif_pos)
                                        < (error_df["N_build_end"] + 2 * av_readlen)
                                    )
                                ].index.tolist()
                                temp_dif_pos += number + alignments[0][4]
                                if idx:
                                    error_df.at[idx[0], "status"] = "fixed"
                                    error_df.loc[
                                        (error_df["contig"] == contig_name)
                                        & (
                                            error_df["order"]
                                            > error_df.at[idx[0], "order"]
                                        ),
                                        "N_build_start",
                                    ] = error_df["N_build_start"] - (
                                        number + alignments[0][4]
                                    )  # temp_dif_pos
                                    error_df.loc[
                                        (error_df["contig"] == contig_name)
                                        & (
                                            error_df["order"]
                                            > error_df.at[idx[0], "order"]
                                        ),
                                        "N_build_end",
                                    ] = error_df["N_build_end"] - (
                                        number + alignments[0][4]
                                    )  # temp_dif_pos
                                else:
                                    error_df.loc[
                                        (error_df["contig"] == contig_name)
                                        & (error_df["N_build_start"] > (start - 13)),
                                        "N_build_start",
                                    ] = error_df["N_build_start"] - (
                                        number + alignments[0][4]
                                    )  # temp_dif_pos | (start - temp_dif_pos - 12)
                                    error_df.loc[
                                        (error_df["contig"] == contig_name)
                                        & (error_df["N_build_start"] > (start - 13)),
                                        "N_build_end",
                                    ] = error_df["N_build_end"] - (
                                        number + alignments[0][4]
                                    )  # - temp_dif_pos

                            else:
                                temp_dif_pos += number + alignments[0][4]

        ### Creating new target based on control and N_target ######
        if control_index:
            for index, number in control_index:
                for i in range(0, index):
                    new_target.append(
                        list(
                            [
                                N_target.get(contig_name)[i][0],
                                N_target.get(contig_name)[i][1],
                            ]
                        )
                    )
                for i in range(index, len(N_target.get(contig_name))):
                    new_target.append(
                        list(
                            [
                                N_target.get(contig_name)[i][0] - number,
                                N_target.get(contig_name)[i][1] - number,
                            ]
                        )
                    )
        else:
            if (N_target.get(contig_name)):
                new_target = N_target.get(contig_name)

        if not user_file:
            
            with open(os.path.join(output_dir, "tmp", "target"), "a") as file_target:
                for start, end in new_target:
                    file_target.write(
                        contig_name + "\t" + str(start) + "\t" + str(end) + "\n"
                    )

        new_fasta_temp.write(">{}\n{}\n".format(contig_name, new_subfasta))

    new_fasta_temp.close()
    if not user_file:
        file_target.close()
    if fixed != "":
        return fixed


def check_npos(av_readlen, one_fasta, error=False):
    """Check Npos from contigs"""
    ### Personal remember - genome pos (1..N), list pos (0..N)
    ###taking all N index
    lst = list(one_fasta)
    indexes = []
    start = 0
    while True:
        try:
            start = lst.index("N", start)
            indexes.append(start)
            start += 1
        except ValueError:
            break

    check_new = False
    final_list = list()

    ###tracking all pos btw the first and last value from a gap (individuals N as well)
    if indexes:
        for index, pos in enumerate(indexes):
            if index == 0:
                final_list.append(pos)
                pos_temp = pos
                continue
            else:
                if pos == pos_temp + 1:
                    check_new = False
                else:
                    final_list.append(pos_temp)
                    if not check_new:
                        final_list.append(pos)
                        check_new = True
                    else:
                        final_list.append(pos)
                        check_new = False
                pos_temp = pos

        final_list.append(pos)
        ### merging regions if they are too close (80% of readlength)
        temp_list = list()
        for i in range(1, len(final_list) - 1, 2):
            if final_list[i + 1] - final_list[i] < av_readlen * 0.8:
                temp_list.extend([i, i + 1])

        for item in sorted(temp_list, reverse=True):
            del final_list[item]

        ###Final N_pos - list with startN, end N, Number of N in total
        N_pos = list()
        N_target = list()

        for i in range(0, len(final_list), 2):
            N_pos.append(
                tuple(
                    (
                        final_list[i],
                        final_list[i + 1],
                        final_list[i + 1] - final_list[i] + 1,
                    )
                )
            )  ## return the python index mode
            N_target.append(
                list([final_list[i] + 1, final_list[i + 1] + 1])
            )  ### return the position for calculate new target

        dict_npos = {one_fasta.id: N_pos}
        dict_ntarget = {one_fasta.id: N_target}

        return dict_npos, dict_ntarget
    else:
        return list(), list()


def filtering_bam(output_dir, thread, num_mm, bam_sorted, r1, r2, r12):
    """Filtering in reads with less than [2] mismatch.
    Step used to keep only interesting reads to reduce the fastq's size"""
    samfile = ps.AlignmentFile(bam_sorted + "_sorted.bam", "rb")
    match_reads = list()
    for read in samfile.fetch():
        if not read.has_tag("NM"):
            continue
        if not (
            read.get_tag("NM") > num_mm
        ):  ## if TAG NM is not greater than num_mm (default = 1) -> read_list to be filtered
            match_reads.append(read.qname)
    samfile.close()

    ### creating a file with unique reads names that fulfill the criteria above
    match_reads = list(set(match_reads))
    with open(output_dir + "/tmp/matched_reads", "w") as outfile:
        outfile.write("\n".join(str(item) for item in match_reads))
    outfile.close()

    ### Filtering original fastq with the high stringency reads
    subprocess.run(
        [
            os.path.join("filterbyname.sh"),
            "in=" + r1,
            "in2=" + r2,
            "out=" + output_dir + "/tmp/res_R1.fastq",
            "out2=" + output_dir + "/tmp/res_R2.fastq",
            "names=" + output_dir + "/tmp/matched_reads",
            "include=t",
            "overwrite=t",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def build_N(
    output_dir,
    fasta,
    av_readlen,
    organized_errors,
):
    """Replace error for N and extend the contig's edges for N - pre-curation step"""

    if os.path.exists(
        os.path.join(output_dir, "tmp", "target")
    ):  ### remove target file if it exists
        os.remove(os.path.join(output_dir, "tmp", "target"))

    # #---------------- Replace region for N ---------------------------
    dict_error_pos = defaultdict(list)
    dict_only_errors = defaultdict(list)
    dict_len = defaultdict()
    ext_size = av_readlen * 3
    set_keys = set(organized_errors.keys())

    fasta_N = open(os.path.join(output_dir, "tmp", "v_0.fasta"), "w")
    target = open(os.path.join(output_dir, "tmp", "target"), "w")

    for seq_record in SeqIO.parse(fasta, "fasta"):
        dict_len[seq_record.id] = len(seq_record)
        seq_mutable = ""

        if seq_record.id in set_keys:
            set_keys.remove(seq_record.id)
            seq_mutable = seq_record.seq.tomutable()
            count = 0
            target.write("{}\t{}\t{}\n".format(seq_record.id, "1", ext_size))

            for j, (start, end) in enumerate(
                organized_errors[seq_record.id]
            ):  ## Add N after the end position
                seq_mutable[start - 1 + count : end + count] = (end - start + 1) * "N"
                dict_error_pos[seq_record.id].append(
                    list([start + (ext_size) * (j + 1), end + (ext_size) * (j + 2)])
                )
                seq_mutable[end + count : end + count] = "N" * ext_size
                count += ext_size
                target.write(
                    "{}\t{}\t{}\n".format(
                        seq_record.id, start + count, end + count + ext_size
                    )
                )

                if seq_record.seq[round((start + end) / 2)] != "N":
                    dict_only_errors[seq_record.id].append(tuple([start, end]))

            ### This is part of extension script - It adds av_readlen*3 "N"s in the beginning and end from a file
            seq_mutable[0:0] = "N" * ext_size
            seq_mutable.extend(("N" * ext_size))

            len_mut = len(seq_mutable)
            target.write(
                "{}\t{}\t{}\n".format(seq_record.id, len_mut + 1 - ext_size, len_mut)
            )
            fasta_N.write(">{}\n{}\n".format(seq_record.id, seq_mutable))
        else:
            ### This is part of extension script - It adds av_readlen*3 "N"s in the beginning and ending from a file
            seq_mutable = seq_record.seq.tomutable()
            len_mut = len(seq_mutable)
            seq_mutable[0:0] = "N" * ext_size
            seq_mutable.extend(("N" * ext_size))
            target.write("{}\t{}\t{}\n".format(seq_record.id, "1", ext_size))
            target.write(
                "{}\t{}\t{}\n".format(seq_record.id, len_mut + 1 - ext_size, len_mut)
            )
            fasta_N.write(">{}\n{}\n".format(seq_record.id, seq_mutable))

        dict_error_pos[seq_record.id].append(list([1, ext_size]))
        dict_error_pos[seq_record.id].append(
            list([len(seq_mutable) + 1 - (ext_size), len(seq_mutable)])
        )
        dict_error_pos[seq_record.id] = sorted(
            dict_error_pos[seq_record.id], key=lambda i: i[0]
        )

    fasta_N.close()
    target.close()

    L = [(k, *t) for k, v in dict_only_errors.items() for t in v]
    error_df = pd.DataFrame(
        L,
        columns=[
            "contig",
            "start",
            "end",
        ],
    )
    error_df["order"] = error_df.groupby("contig").cumcount() + 1
    error_df["len"] = error_df["contig"].map(dict_len)
    error_df["type_of_error"] = "local_assembly_error"
    error_df.loc[
        (error_df["start"] == 1) | (error_df["end"] == error_df["len"]), "type_of_error"
    ] = "edge_reads_cov"
    error_df["N_build_start"] = error_df["start"] + (ext_size * error_df["order"])
    error_df["N_build_end"] = (
        error_df["end"] + (ext_size * error_df["order"]) + ext_size
    )
    error_df["status"] = ""
    logger.warning(
        "\n\nFixame could detect a total of {} errors in {} contig(s)\n * Local Assembly Errors - {}\n * Edge reads coverage - {}\n".format(
            len(error_df),
            len(error_df["contig"].unique()),
            len(error_df[error_df["type_of_error"] == "local_assembly_error"]),
            len(error_df[error_df["type_of_error"] == "edge_reads_cov"]),
        )
    )

    return dict_len, error_df


def var_cal_fix(output_dir, count, fixed, thread, dp_cov, av_readlen, error_df):
    """Find changes on Ns position and replace them"""

    aligner(
        output_dir,
        thread,
        0,
        os.path.join(output_dir, "tmp", "v_" + str(count - 1) + ".fasta"),
        r1=os.path.join(output_dir, "tmp", "res_R1.fastq"),
        r2=os.path.join(output_dir, "tmp", "res_R2.fastq"),
        r12="",
        bam_out="v_" + str(count - 1),
        semi=True,
    )

    # Run variant finder - bcftools

    comm = "bcftools mpileup -R {}/target -A -B -f {}/v_{}.fasta {}/v_{}_sorted.bam".format(
        os.path.join(output_dir, "tmp"),
        os.path.join(output_dir, "tmp"),
        count - 1,
        os.path.join(output_dir, "tmp"),
        count - 1,
    )
    varfind = subprocess.Popen(
        comm,
        shell=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        stderr=subprocess.DEVNULL,
    )
    output = varfind.communicate()[0]
    varfind.stdout.close()

    # Process output based on dp_cov -default[1]
    var_list = dict()
    for line in output.splitlines():
        if re.match("^#", line):
            continue
        line_list = line.split(
            "\t"
        )  # 0=header;1=pos;3=ref;4=var("," slipted);7=flags(;/= splited)

        if line_list[0] in var_list.keys():
            if int(re.search("DP=(.+?);", line_list[7]).group(1)) >= dp_cov:
                if (
                    (line_list[4].split(",")[0].upper() == "A")
                    or (line_list[4].split(",")[0].upper() == "C")
                    or (line_list[4].split(",")[0].upper() == "T")
                    or (line_list[4].split(",")[0].upper() == "G")
                ):
                    var_list[line_list[0]].append(
                        (
                            int(line_list[1]) - 1,
                            len(line_list[3]),
                            len(line_list[4].split(",")[0]),
                            line_list[4].split(",")[0],
                        )
                    )
        else:
            if int(re.search("DP=(.+?);", line_list[7]).group(1)) >= dp_cov:
                if (
                    (line_list[4].split(",")[0].upper() == "A")
                    or (line_list[4].split(",")[0].upper() == "C")
                    or (line_list[4].split(",")[0].upper() == "T")
                    or (line_list[4].split(",")[0].upper() == "G")
                ):
                    var_list[line_list[0]] = [
                        (
                            int(line_list[1]) - 1,
                            len(line_list[3]),
                            len(line_list[4].split(",")[0]),
                            line_list[4].split(",")[0],
                        )
                    ]

    new_fasta_temp = open(os.path.join(output_dir, "tmp", "snp_fasta.fasta"), "w")
    for seq_record in SeqIO.parse(
        os.path.join(output_dir, "tmp", "v_" + str(count - 1) + ".fasta"), "fasta"
    ):
        pos_temp = 0
        if seq_record.id in var_list.keys():
            seq_mutable = seq_record.seq.tomutable()
            for value in var_list[seq_record.id]:
                pos, len_ref, len_var, var = value
                seq_mutable[pos + pos_temp : pos + pos_temp + len_ref] = var
                pos_temp += len_var - len_ref
            new_fasta_temp.write(">{}\n{}\n".format(seq_record.id, seq_mutable))
        else:
            new_fasta_temp.write(
                ">{}\n{}\n".format(seq_record.id, seq_record.seq)
            )  ##removi seq_mutable e substitui por seq_record.seq

    new_fasta_temp.close()

    fixed = check_overlap(
        output_dir,
        os.path.join(output_dir, "tmp", "snp_fasta.fasta"),
        av_readlen,
        thread,
        error_df=error_df,
        fixed=fixed,
        count=count,
    )
    return fixed


def remove_N(
    output_dir,
    name_fasta,
    fasta_semifinal,
    organized_errors,
    fasta_len,
    av_readlen,
    mean_gap,
    mean_gap_std,
    thread,
):
    fasta_wo_N = ""
    ext_size = av_readlen * 3
    unordered_fasta = open(
        os.path.join(output_dir, "tmp", name_fasta + "_unordered.fasta"), "w"
    )

    # Alignment needed to selected reads #fasta_semifinal
    aligner(
        output_dir,
        thread,
        1,
        fasta_semifinal,
        r1=os.path.join(output_dir, "tmp", "res_R1.fastq"),
        r2=os.path.join(output_dir, "tmp", "res_R2.fastq"),
        r12="",
        bam_out="check_read",
        semi=True,
    )
    to_be_checked = []

    for seq_record in SeqIO.parse(fasta_semifinal, "fasta"):
        actual_start, actual_end = 0, 0
        if seq_record.id in organized_errors.keys():
            to_be_checked.append(seq_record)
        else:
            seq_mutable = seq_record.seq.tomutable()
            fasta_wo_N = remove_N_slave(
                seq_record.id, seq_mutable, actual_start, actual_end, av_readlen
            )
            unordered_fasta.write(str(fasta_wo_N))

    with ProcessPoolExecutor(thread) as executor:
        check_npos_result = executor.map(partial(check_npos, av_readlen), to_be_checked)

    npos_list = list()

    for npos_unit in check_npos_result:
        npos_list.append(npos_unit[0])

    N_pos = {}

    for i in npos_list:
        N_pos.update(i)

    for checked in to_be_checked:
        actual_start, actual_end = 0, 0
        contig = checked.id
        seq_mutable = checked.seq.tomutable()
        
        if N_pos.get(contig):
            check_edges = N_pos[contig][1:-1]
            if check_edges:
                seq_mutable = check_reads_N_edges(
                    output_dir,
                    contig,
                    seq_mutable,
                    contig,
                    av_readlen,
                    mean_gap,
                    mean_gap_std,
                    check_edges,
                )

            # Edges
            if organized_errors[contig][0][0] == 1:
                actual_start = organized_errors[contig][0][1] + ext_size
            if organized_errors[contig][-1][1] == fasta_len[contig]:
                actual_end = (
                    organized_errors[contig][-1][1]
                    - organized_errors[contig][-1][0]
                    + 1
                    + ext_size
                )
            fasta_wo_N = remove_N_slave(
                contig, seq_mutable, actual_start, actual_end, av_readlen
            )
        else: 
            fasta_wo_N = ">" + contig + "\n" + seq_mutable + "\n"
        unordered_fasta.write(str(fasta_wo_N))
    unordered_fasta.close()
    records = list(
        SeqIO.parse(
            os.path.join(output_dir, "tmp", name_fasta + "_unordered.fasta"), "fasta"
        )
    )
    records.sort(key=lambda r: -len(r))
    SeqIO.write(
        records,
        os.path.join(output_dir, "fixame_results", name_fasta + "_fixame.fasta"),
        "fasta",
    )


def check_reads_N_edges(
    output_dir,
    contig_name,
    seq_mutable,
    seq_name,
    av_readlen,
    mean_gap,
    mean_gap_std,
    N_pos,
):
    r_left = os.path.join(output_dir, "tmp", "r_left")
    r_right = os.path.join(output_dir, "tmp", "r_right")
    left_right = os.path.join(output_dir, "tmp", "left_right")
    count = 0

    mean_frag_len = 2 * av_readlen + mean_gap
    no_support = open(
        os.path.join(output_dir, "fixing_log", "fixame_without_read_support.txt"), "a"
    )
    for start, end, space in N_pos:
        check_valid = ""

        if (start - mean_gap) < 1:
            start_mgap = 1
        else:
            start_mgap = start - mean_gap
        if (end + mean_gap) > len(seq_mutable):
            end_mgap = len(seq_mutable)
        else:
            end_mgap = end + mean_gap

        cmd = """samtools view {}/check_read_sorted.bam {}:{}-{} | \
                grep -v "*" |cut -f 1 | sort | uniq -u > {}""".format(
            os.path.join(output_dir, "tmp"),
            seq_name,
            (start_mgap + count),
            (start + count),
            r_left,
        )
        subprocess.run(
            cmd,
            shell=True,
        )
        cmd = """samtools view {}/check_read_sorted.bam {}:{}-{} | \
                grep -v "*" |cut -f 1 | sort | uniq -u > {}""".format(
            os.path.join(output_dir, "tmp"),
            seq_name,
            (end + count),
            (end + mean_gap + count),
            r_right,
        )
        subprocess.run(
            cmd,
            shell=True,
        )
        cmd = """cat {} {}| sort | uniq -d > {}""".format(r_left, r_right, left_right)
        subprocess.run(
            cmd,
            shell=True,
        )

        if os.stat(left_right).st_size == 0:
            no_support.write(
                "No read support around:\t{}:{}-{}\n".format(
                    contig_name,
                    start + count - (3 * av_readlen),
                    end + count - (3 * av_readlen),
                )
            )
            continue
        else:
            # print (seq_name,(start-2*av_readlen), (end+2*av_readlen), left_right)
            cmd = """samtools view {}/check_read_sorted.bam {}:{}-{} | grep -F -f {}| awk -v FS="\\t" '$9 > 0 {{print}}' """.format(
                os.path.join(output_dir, "tmp"),
                seq_name,
                (start_mgap + count),
                (end_mgap + 1 + count),
                left_right,
            )
            check_valid = subprocess.check_output(
                cmd, universal_newlines=True, shell=True
            )

            if check_valid == "":
                continue

            cmd = """samtools view {}/check_read_sorted.bam {}:{}-{} | grep -F -f {}| awk -v FS="\\t" '$9 > 0 {{ sum += $9; n++ }} END {{print int(sum/n)}}' """.format(
                os.path.join(output_dir, "tmp"),
                seq_name,
                (start_mgap + count),
                (end_mgap + 1 + count),
                left_right,
            )

            distance = int(
                subprocess.check_output(
                    cmd, universal_newlines=True, shell=True
                ).split()[0]
            )
            if distance == 0:
                continue

            if distance < (mean_frag_len - mean_gap_std):  # rare but possible
                if space <= (0.1 * av_readlen):
                    seq_start = int(start + count)
                    seq_end = int(start + count)
                else:
                    seq_start = int(start + count + (0.1 * av_readlen))
                    seq_end = int(start + count + (0.1 * av_readlen))

                seq_mutable[seq_start:seq_end] = (mean_frag_len - distance) * "N"
                count += mean_frag_len - distance

            elif distance > (mean_frag_len + mean_gap_std):
                seq_start = int(start + count + (0.1 * av_readlen))
                seq_end = int(
                    start + count + (0.1 * av_readlen) + (distance - mean_frag_len)
                )
                if space <= (seq_end - seq_start):
                    continue
                seq_mutable[seq_start:seq_end] = ""
                count -= distance - mean_frag_len

            else:
                continue

    no_support.close()
    return seq_mutable


def remove_N_slave(fasta_header, seq_mutable, actual_start, actual_end, av_readlen):
    len_seq = len(seq_mutable)
    new_fasta = ""
    ext_size = av_readlen * 3
    for number in range(len_seq - (ext_size) - actual_end, len_seq):
        if seq_mutable[number] == "N":
            remov_end = number
            break
        else:
            remov_end = len_seq  ## this considers if all edges were fullfilled by bases (same as below)
    for number in range((ext_size) + actual_start - 1, 0, -1):
        if seq_mutable[number] == "N":
            remov_start = number + 1  ### +1?
            break
        else:
            remov_start = 0
    new_fasta += ">" + fasta_header + "\n" + seq_mutable[remov_start:remov_end] + "\n"

    return new_fasta


def final_output(output, df, features):
    df.to_csv(
        os.path.join(output, "Fixame_AssemblyErrors_report.txt"),
        columns=[
            "contig",
            "start",
            "end",
            "type_of_error",
            "status",
            "sample_name",
        ],
        index=None,
        sep="\t",
    )

    # grouped = df.groupby(['sample_name', 'contig', 'type_of_error', 'status'])[['type_of_error']].agg('count')
    # g_unstack = grouped.unstack(['status','type_of_error']).reset_index()
    # g_unstack.columns = ['sample_name','contig','Edges_events', 'Local_error_events', 'Fixed_local_error']
    # g_unstack[['Edges_events','Local_error_events','Fixed_local_error']] = g_unstack[['Edges_events','Local_error_events','Fixed_local_error']].fillna(0).astype(int)
    # g_unstack.to_csv(os.path.join(output, "Fixame_Summary.txt"), index=None, sep='\t')

    extra_features = open(os.path.join(output, "Fixame_features.txt"), "w+")
    if features:
        extra_features.write("contig\tfeature\tcount\n")
        for item in features:
            for ctg_name, v in item.items():
                for feature, v_two in v.items():
                    extra_features.write(f"{ctg_name}\t{feature}\t{v_two}\n")
    extra_features.close()
