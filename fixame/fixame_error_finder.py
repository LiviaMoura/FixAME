import argparse
import os
import logging
import sys
import fixame
import pysam as ps
import pandas as pd
from collections import defaultdict
from xopen import xopen
from concurrent.futures import ProcessPoolExecutor
from Bio.Seq import reverse_complement
from xopen import xopen
from Bio import SeqIO
import pandas as pd
import subprocess
import regex

__author__ = "Rohan Sachdeva"
__copyright__ = "Copyright 2020"
__maintainer__ = "Livia Moura, Rohan Sachdeva"
__email__ = "liviam.moura@gmail.com, rohansach@berkeley.edu"
__status__ = "Development"
   
def calculate_reference_lengths(fasta, minimum_assembly_length):
    "Calculate lengths of sequences and return a dictionary"

    length_dict = {}
    for record in SeqIO.parse(xopen(fasta), "fasta"):
        length = len(record.seq)
        if length >= minimum_assembly_length:
            id_ = record.id
            length_dict[id_] = length
    return length_dict


def check_palindrome(sequence):
    "Check sequence that is a full palindrome"

    sequence_rc = reverse_complement(sequence)
    if sequence == sequence_rc:
        length = len(sequence)

        return length


def check_direct_repeat(sequence):
    "Check sequence for concatenated repeats"

    length = len(sequence)
    fraction_length = int(0.2 * length)
    fraction_seq = sequence[:fraction_length]
    substitutions = "2"
    regex_subs = "{s<=" + substitutions + "}"
    rx = f"({fraction_seq})" + regex_subs
    fraction_count = regex.findall(rx, sequence)
    fraction_count = len(fraction_count)

    if fraction_count >= 2:
        return length


def check_circular(sequence, overlap_length):
    "Check for potentially circular sequences"

    seq_length = len(sequence)
    half_seq_length = int(seq_length / 2)
    beg_seq, end_seq = sequence[:50], sequence[-half_seq_length:]
    beg_seq_in_end_seq_index = end_seq.rfind(beg_seq)

    if beg_seq_in_end_seq_index != -1:
        end_match = end_seq[beg_seq_in_end_seq_index:]
        len_end_match = len(end_match)
        beg_match = sequence[:len_end_match]

        if beg_match == end_match:
            return len_end_match


def check_direct_features(id_sequence):
    "Check sequences for palindromes, direct repeats, and potentially circular sequences"

    id_, sequence = id_sequence

    palindrome_length = check_palindrome(sequence)
    if palindrome_length:
        id_feature = (id_, "palindrome", palindrome_length)

        return id_feature

    else:
        direct_repeat_length = check_direct_repeat(sequence)

        if direct_repeat_length:
            id_feature = (id_, "direct_repeat", direct_repeat_length)

            return id_feature

        else:

            potential_circular_length = check_circular(sequence, 50)
            if potential_circular_length:
                id_feature = (id_, "potential_circular", potential_circular_length)
                return id_feature


def check_direct_features_parallel(fasta, threads):
    "Run direct_features in parallel"

    direct_features_dict = defaultdict(dict)
    reference_tuples = []

    for record in SeqIO.parse(xopen(fasta), "fasta"):
        id_sequence = (record.id, str(record.seq))
        reference_tuples.append(id_sequence)

    with ProcessPoolExecutor(threads) as executor:
        execute_result = executor.map(check_direct_features, reference_tuples)

        for id_feature in execute_result:
            if id_feature:
                id_ = id_feature[0]
                feature, feature_length = id_feature[1], id_feature[2]

                direct_features_dict[id_][feature] = feature_length

    return direct_features_dict


def parse_map(bam_file, num_mm, threads, minimum_assembly_length, reference_to_length):
    "Parse bam file to primarily return information on read locations and mismatches. Also, calculates info on reads and pairing."

    template_lengths, read_lengths = [], []
    bam_dict, reference_read_lengths = defaultdict(list), defaultdict(list)
    bam_parsed = ps.AlignmentFile(bam_file, "rb", threads=threads)
    
    for read in bam_parsed:
        # if read.is_paired and reference == read.next_reference_name and read.is_proper_pair and read.get_tag('NM') <= num_mm:
        # if read.is_paired and get_tag('NM') <= num_mm:
        reference = read.reference_name
        reference = reference.split()[0]
        reference_length = reference_to_length[reference]

        if reference_length >= minimum_assembly_length and not read.is_unmapped:
            mismatches, query_length = read.get_tag("NM"), read.query_alignment_length
            fraction_id = 1 - mismatches / query_length

            if fraction_id >= 0.95:
                read_name = read.query_name
                # mapq = read.mapping_quality
                template_length = read.template_length
                reference_positions = read.get_reference_positions(full_length=False)
                start_pos = reference_positions[0] + 1
                end_pos = reference_positions[-1] + 1
                proper_pair = read.is_proper_pair
                info = (
                    read_name,
                    reference,
                    start_pos,
                    end_pos,
                    proper_pair,
                    mismatches,
                )
                bam_dict[reference].append(info)

                if mismatches <= num_mm:
                    template_lengths.append(abs(template_length))
                    reference_read_lengths[reference].append(query_length)
                    read_lengths.append(query_length)

    average_template_length = pd.Series(template_lengths).mean()
    #average_read_length = pd.Series(read_lengths).mean()

    ## few bam metrics
    cmd = '''samtools stats {} | grep "^SN" | cut -f 2- '''.format(bam_file)

    output = subprocess.check_output(cmd, shell=True, universal_newlines=True)

    for line in output.splitlines():
        array_line = line.split(":")
        if (array_line[0] == 'average length'):
            average_read_length = int(float(array_line[1].strip()))
        if (array_line[0] == 'insert size average'):
            average_gap_length = int(float(array_line[1].strip()))
        if (array_line[0] == 'insert size standard deviation'):
            average_gap_std = int(float(array_line[1].strip()))
            
    #print(mean_gap, mean_gap_dev)
    #average_template_length = int(round(average_template_length))
    #average_read_length = int(round(average_read_length))

    var_template_length = pd.Series(template_lengths).std()
    #var_read_length = pd.Series(read_lengths).std()

    template_length_min = average_template_length - var_template_length
    template_length_max = average_template_length + var_template_length

    template_length_min = int(template_length_min)
    template_length_max = int(template_length_max)

    #average_gap_length = average_template_length - average_read_length * 2
    #average_gap_length = int(average_gap_length)

    #print(average_template_length,'TEMPLETA_LENTH',average_gap_length,'GAP_LENGH')

    return (
        bam_dict,
        reference_read_lengths,
        average_template_length,
        average_read_length,
        average_gap_length,
        template_length_min,
        template_length_max,
        average_gap_std,
    )


def check_local_assembly_errors(reference):
    "Check for local assembly errors and identify regions with high variability"
    look_len = 5 #standard value
    read_lengths = reference_read_lengths[reference]
    reference_length = reference_to_length[reference]
    reference_total_length = sum(read_lengths)
    reference_coverage = reference_total_length / reference_length
    read_to_jump_references_list = []

    if reference_coverage >= fasta_cov:
        infos = bam_dict[reference]
        good_positions = set()
        reference_position_to_reads, read_to_positions = (
            defaultdict(set),
            defaultdict(set),
        )

        perfect_matches_counts_dict = defaultdict(int)
        mismatch_positions_counts_dict = defaultdict(int)

        read_to_mismatch_positions = {}

        for info in infos:
            (
                read,
                reference,
                start_position,
                end_position,
                proper_pair,
                mismatches,
            ) = info

            if mismatches <= num_mm:
                if (
                    proper_pair
                    or end_position <= template_length_max
                    or end_position >= reference_length - template_length_max
                ):
                    for position in range(start_position, end_position + 1):
                        reference_position_to_reads[position].add(read)
                        read_to_positions[read].add(position)
                # elif end_position <= template_length_max:
                #     parse_read_positions(start_position, end_position)
                #
                # elif end_position >= reference_length - template_length_max:
                #     parse_read_positions(start_position, end_position)

                if mismatches == 0:
                    for position in range(start_position, end_position + 1):
                        perfect_matches_counts_dict[position] += 1
                else:
                    for position in range(start_position, end_position + 1):
                        mismatch_positions_counts_dict[position] += 1
            else:
                for position in range(start_position, end_position + 1):
                    mismatch_positions_counts_dict[position] += 1

        high_mismatch_positions = set()

        reference_positions = set(range(1, reference_length + 1))

        for reference_position in sorted(reference_positions):
            mismatch_count = mismatch_positions_counts_dict[reference_position]
            perfect_count = perfect_matches_counts_dict[reference_position]

            total_count = perfect_count + mismatch_count

            if total_count > 0:
                mismatch_fraction = mismatch_count / total_count

                if mismatch_fraction >= 0.5:
                    high_mismatch_positions.add(reference_position)

            pass_count = 0
            reads = reference_position_to_reads[reference_position]
            back_reference_position = reference_position - look_len
            ahead_reference_position = reference_position + look_len

            if back_reference_position < 1:
                back_reference_position = 1
            if ahead_reference_position > reference_length:
                ahead_reference_position = reference_length
            for read in reads:
                if pass_count > 0:
                    break
                read_positions = read_to_positions[read]
                if (
                    back_reference_position in read_positions
                    and ahead_reference_position in read_positions
                ):
                    good_positions.add(reference_position)
                    pass_count += 1
                    break

        error_positions = set(reference_positions) - good_positions
        error_positions = sorted(error_positions)

        high_mismatch_positions = sorted(high_mismatch_positions)

        return reference, reference_coverage, error_positions, high_mismatch_positions


def check_local_assembly_errors_parallel(references, threads, rrl, rtl, fc, bd, nm, tlm): #, ref_read_len, ref_to_len, fast_cov, bam_dic):
    "Run find_regions in parallel"
    global reference_read_lengths, reference_to_length, fasta_cov, bam_dict, num_mm, template_length_max
    reference_read_lengths, reference_to_length, fasta_cov, bam_dict, num_mm, template_length_max= rrl, rtl, fc, bd, nm, tlm
    reference_to_error_regions, coverage_dict, reference_to_high_mismatch_positions = (
        {},
        {},
        {},
    )

    with ProcessPoolExecutor(threads) as executor:
        execute_result = executor.map(check_local_assembly_errors, references)

        for (
            reference,
            reference_coverage,
            error_positions,
            high_mismatch_positions,
        ) in execute_result:
            reference_to_error_regions[reference] = error_positions

            coverage_dict[reference] = reference_coverage

            reference_to_high_mismatch_positions[reference] = high_mismatch_positions

    return (
        reference_to_error_regions,
        coverage_dict,
        reference_to_high_mismatch_positions,
    )
