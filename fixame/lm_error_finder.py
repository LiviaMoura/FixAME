import argparse
import os
import logging
import sys
import fixame
import pysam as ps
import pandas as pd
from collections import defaultdict
from xopen import xopen
from Bio import SeqIO
from xopen import xopen

__author__ = "Livia Moura"
__copyright__ = "Copyright 2019"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"

class Error_finder():
    
    def reference_info(fasta):
        
        length_dict, gc_dict = {}, {}

        for record in SeqIO.parse(xopen(fasta), 'fasta'):
            id_, length, gc = record.id, len(record.seq), GC(record.seq)
            length_dict[id_] = length
            gc_dict[id_] = gc

            return length_dict, gc_dict

    def parse_bam(bam_file, num_mm):
        #Parse bam file and return a dict of start,end of each read on reference'
        template_lengths = []
        bam_dict,reference_read_lengths = defaultdict(list),defaultdict(list)
        bam_parsed = ps.AlignmentFile(bam_file, 'rb')
    
        for read in bam_parsed:
            reference = read.reference_name
            if read.is_paired and reference == read.next_reference_name and read.is_proper_pair and read.get_tag('NM') <= num_mm:
                reference = reference.split(' ')[0]
                read_name = read.query_name
                query_length,template_length = read.query_length,read.template_length
                reference_positions = read.get_reference_positions(full_length=False)
                start_pos = reference_positions[0] + 1
                end_pos = start_pos + query_length
                info = (read_name,reference,start_pos,end_pos)
                bam_dict[reference].append(info)
                template_lengths.append(abs(template_length))
                reference_read_lengths[reference].append(query_length) 
        mean_template_length = pd.Series(template_lengths).mean()
        return bam_dict,reference_read_lengths,mean_template_length       
        

    def find_regions(fasta_cov,bam_dict,look_len,reference_lengths,reference_read_lengths,mean_template_length):
        reference_to_bad_regions = defaultdict(set)
        all_read_lengths,coverage_dict = [],{}
        for reference,read_lengths in reference_read_lengths.items():
            all_read_lengths.extend(read_lengths)
            reference_length = reference_lengths[reference]
            reference_total_bp = sum(read_lengths)
            reference_coverage = reference_total_bp / reference_length
            coverage_dict[reference] = reference_coverage
        
            if reference_coverage >= fasta_cov:
                infos = bam_dict[reference]
                all_positions,good_positions = set(),set()
                reference_position_to_reads,read_to_positions = defaultdict(set),defaultdict(set)
                for info in infos:
                    read,reference,start_position,end_position = info
                    for position in range(start_position,end_position + 1):
                        reference_position_to_reads[position].add(read)
                        read_to_positions[read].add(position)
                        all_positions.add(position)
                reference_positions = set(range(1,reference_length + 1))
                for reference_position in reference_positions:
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
                        if back_reference_position in read_positions and ahead_reference_position in read_positions:
                            good_positions.add(reference_position)
                            pass_count += 1
                            break
                bad_positions = set(reference_positions) - good_positions
                bad_positions = sorted(bad_positions)
                if bad_positions:
                    reference_to_bad_regions[reference] = bad_positions
        mean_read_length = int(round(pd.Series(all_read_lengths).mean()))
        mean_gap_length = int(round(mean_template_length - mean_read_length * 2))
        return reference_to_bad_regions,coverage_dict,mean_read_length,mean_gap_length
