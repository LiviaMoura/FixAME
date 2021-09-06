import os
import logging
import sys
from datetime import datetime 
import pysam as ps
import pandas as pd
from collections import defaultdict
from xopen import xopen
import traceback
from concurrent.futures import ProcessPoolExecutor
from Bio.Seq import reverse_complement
from xopen import xopen
from Bio import SeqIO
import subprocess
import regex
import shutil


from fixame.fixame_common import *
from fixame.fixame_aligner import aligner

__author__ = "Rohan Sachdeva"
__copyright__ = "Copyright 2021"
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
    look_len = 5  # standard value
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

    else:
        error_positions = set()
        high_mismatch_positions = set()

    return reference, reference_coverage, error_positions, high_mismatch_positions


def check_local_assembly_errors_parallel(references, threads, rrl, rtl, fc, bd, nm, tlm): #, ref_read_len, ref_to_len, fast_cov, bam_dic):
    "Run find_regions in parallel"
    global reference_read_lengths, reference_to_length, fasta_cov, bam_dict, num_mm, template_length_max
    reference_read_lengths, reference_to_length, fasta_cov, bam_dict, num_mm, template_length_max = rrl, rtl, fc, bd, nm, tlm
    reference_to_error_regions, coverage_dict, reference_to_high_mismatch_positions = (
        {},
        {},
        {},
    )

    with ProcessPoolExecutor(threads) as executor:
        execute_result = executor.map(check_local_assembly_errors, references)
        #print(list(execute_result), 'EXECUTE RESULT\n')
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


def organizing_found_errors(av_readlen, dict_replace_0):
    dict_replace = {key: value for key, value in dict_replace_0.items() if value}
    dict_replace_chunked = defaultdict(list)     
    for key, value in dict_replace.items():        
        check_new = False
        final_list = defaultdict(list)        
        ###tracking all pos btw the first and last value from a gap (individuals N as well)
        for index, pos in enumerate(dict_replace[key]):
            if index == 0:
                final_list[key].append(pos)
                pos_temp = pos                
                continue                
            else:                
                if (pos == pos_temp + 1):
                    check_new = False
                else:
                    final_list[key].append(pos_temp)
                    if not check_new:
                        final_list[key].append(pos)
                        check_new = True
                    else:                       
                        final_list[key].append(pos)
                        check_new = False
                pos_temp = pos
        
        final_list[key].append(pos)

        ### merging regions if they are too close (80% of readlength)
        temp_list = list()
        for i in range(1,len(final_list[key])-1,2):
            if final_list[key][i+1] - final_list[key][i] < (av_readlen*0.8):
                temp_list.extend([i,i+1])

        for item in sorted(temp_list, reverse=True):
            del final_list[key][item]
        
        for i in range(0, len(final_list[key]), 2):
            dict_replace_chunked[key].append(tuple(final_list[key][i:i+2]))
    return dict_replace_chunked


def main(**kwargs):
    method = common_validate(**kwargs)

    if (kwargs.get('min_ctg_len') < 800):
        logging.info('Checking minimum contig length')
        logging.error("--min_ctg_len must be >= 800")
        sys.exit()
    else:
        minimum_assembly_length = kwargs.get('min_ctg_len')

    if kwargs.get('r12'):
            read12_in = os.path.realpath(os.path.expanduser(kwargs.get('r12')))
            read1_in =''
            read2_in =''
    else:    
        read1_in = os.path.realpath(os.path.expanduser(kwargs.get('r1')))
        read2_in = os.path.realpath(os.path.expanduser(kwargs.get('r2')))
        read12_in = ''

    # av_readlen = temp_average_read(kwargs.get('r1'))
    output_dir = kwargs.get('output_dir')

    try:
        mydir = os.path.join(output_dir,'fixame_'+datetime.now().strftime('%Y-%b-%d_%H-%M-%S'))
        os.mkdir(mydir)
        os.mkdir(os.path.join(mydir,'tmp'))
        logging_config(mydir)
        logging.info("Fixame output folder created - {}".format(mydir))
    except:
        logging.error("It wasn't possible to create fixame output folder")
    
    if method == 0:
        fasta_in = os.path.realpath(os.path.expanduser(kwargs.get('fasta')))
        name_fasta = os.path.splitext(os.path.basename(fasta_in))[0]

        logging.info("\n --- Analysing the file {} ---\n".format(name_fasta))

        try:
            logging.info("Mapping reads against the reference")
            aligner(mydir,kwargs.get('threads'),kwargs.get('minid'), fasta_in, r1=read1_in, r2=read2_in, r12=read12_in, bam_out=name_fasta)                
        except:
            logging.error("Something went wrong")
            print(traceback.format_exc())
            sys.exit()
        
        fasta_cov, num_mm = kwargs.get('fasta_cov'), kwargs.get('num_mismatch')
        
        try:
            logging.info("Generating some metrics to keep running")
            reference_to_length = calculate_reference_lengths(fasta_in, minimum_assembly_length)

            (   bam_dict, 
                reference_read_lengths,
                average_template_length,
                average_read_length,
                average_gap_length,
                template_length_min,
                template_length_max,
                average_gap_std
            ) = parse_map(mydir+'/tmp/'+name_fasta+'_sorted.bam', num_mm, kwargs.get('threads'), minimum_assembly_length, reference_to_length)
        except:
            logging.error("Something went wrong")
            sys.exit()
        
        try:
            logging.info("Trying to find regions with local assembly errors")
            (
                reference_to_error_regions,
                coverage_dict,
                reference_to_high_mismatch_positions,
            ) = check_local_assembly_errors_parallel(reference_to_length.keys(), kwargs.get('threads'), reference_read_lengths, reference_to_length, fasta_cov, bam_dict, num_mm,template_length_max)                
        except Exception: 
            logging.error("Something went wrong")  
            traceback.print_exc()
            sys.exit()
        
        organized_errors = organizing_found_errors(average_read_length, reference_to_error_regions)

        with open(os.path.join(output_dir, 'Fixame_AssemblyErrors_report.txt'), 'w+') as error_loc:
            error_loc.write("contig_name\terror_start\terror_end\tn_affected_bases\ttype_of_error\n")
            counter_err, counter_contigs = 0,0
            for key,value in organized_errors.items():
                logging.debug("{} possibly local errors at contig {}".format(len(value), key))
                counter_err += len(value)
                counter_contigs += 1
                for item in value:
                    error_loc.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(key,item[0],item[1],abs(item[1]-item[0]),'Local_error'))
            error_loc.close()
        logging.warning("\n\nFixame could detect a total of {} errors in {} contig(s)\n".format(counter_err,counter_contigs))
        logging.info("The file containing the detected errors {} was created".format(output_dir+'/Fixame_Errors_report.txt'))
        if (kwargs.get('keep') == False):
            try:
                logging.info("Removing temporary files")
                shutil.rmtree(os.path.join(mydir,'tmp'))
            except:
                logging.info("It wasn't possible to remove the /tmp folder")
        shutil.rmtree(os.path.join(mydir,'new_fasta'))
        logging.info("\n\nFixame error_finder proccess done!\n")
    
    else: # BINS MODE
        fasta_array = []
        name_sample = 'bins'
        contigs_bins = open(os.path.join(mydir, 'bin_contigs.txt'), 'w+')
        merged = open(os.path.join(mydir, 'tmp', 'bins.fasta'), 'w+')
        for sample in os.listdir(os.path.realpath(os.path.expanduser(kwargs.get('bins')))):
            name = sample.split('.')[0]
            if sample.lower().endswith(".fasta") or sample.lower().endswith(".fa") or sample.lower().endswith(".fna"):
                fasta_array.append(sample)
                for seq_record in SeqIO.parse(os.path.join(kwargs.get('bins'), sample),'fasta'):
                    contigs_bins.write("{}\t{}\n".format(name, seq_record.id))
                with open(os.path.join(kwargs.get('bins'), sample), 'r') as readfile:
                    shutil.copyfileobj(readfile, merged)
        contigs_bins.close()
        merged.close()

        fasta_in = os.path.join(mydir,'tmp','bins.fasta')
        name_sample = 'bins'

        logging.info("\n --- Analysing a metagenome sample with {} bins ---\n".format(len(fasta_array)))

        try:
            logging.info("Mapping reads against the reference")
            aligner(mydir,kwargs.get('threads'),kwargs.get('minid'), fasta_in, r1=read1_in, r2=read2_in, r12=read12_in, bam_out=name_sample)                
        except:
            logging.error("Something went wrong")
            print(traceback.format_exc())
            sys.exit()
        
        fasta_cov, num_mm = kwargs.get('fasta_cov'), kwargs.get('num_mismatch')
        
        try:
            logging.info("Generating some metrics to keep running")
            reference_to_length = calculate_reference_lengths(fasta_in, minimum_assembly_length)

            (   bam_dict, 
                reference_read_lengths,
                average_template_length,
                average_read_length,
                average_gap_length,
                template_length_min,
                template_length_max,
                average_gap_std
            ) = parse_map(mydir+'/tmp/'+name_sample+'_sorted.bam', num_mm, kwargs.get('threads'), minimum_assembly_length, reference_to_length)
        except:
            logging.error("Something went wrong")
            sys.exit()
        
        try:
            logging.info("Trying to find regions with local assembly errors")
            (
                reference_to_error_regions,
                coverage_dict,
                reference_to_high_mismatch_positions,
            ) = check_local_assembly_errors_parallel(reference_to_length.keys(), kwargs.get('threads'), reference_read_lengths, reference_to_length, fasta_cov, bam_dict, num_mm,template_length_max)                
        except Exception: 
            logging.error("Something went wrong")  
            traceback.print_exc()
            sys.exit()
        
        organized_errors = organizing_found_errors(average_read_length, reference_to_error_regions)

        with open(os.path.join(output_dir, 'Fixame_AssemblyErrors_report.txt'), 'w+') as error_loc:
            error_loc.write("contig_name\terror_start\terror_end\tn_affected_bases\ttype_of_error\n")
            counter_err, counter_contigs = 0,0
            for key,value in organized_errors.items():
                logging.debug("{} possibly local errors at contig {}".format(len(value), key))
                counter_err += len(value)
                counter_contigs += 1
                for item in value:
                    error_loc.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(key,item[0],item[1],abs(item[1]-item[0]),'Local_error'))
            error_loc.close()
        logging.warning("\n\nFixame could detect a total of {} errors in {} contig(s)\n".format(counter_err,counter_contigs))
        logging.info("The file containing the detected errors {} was created".format(output_dir+'/Fixame_Errors_report.txt'))
        
        if (kwargs.get('keep') == False):
            try:
                logging.info("Removing temporary files")
                shutil.rmtree(os.path.join(mydir,'tmp'))
            except:
                logging.info("It wasn't possible to remove the /tmp folder")
        shutil.rmtree(os.path.join(mydir,'new_fasta'))
        
        logging.info("\n\nFixame error_finder proccess done!\n")

if __name__ == "__main__":
    main()
