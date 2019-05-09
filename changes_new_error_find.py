###########################################################################
#project: Genome curation - Mimicking Banfield's genome curation proccess 
#autor: Livia Maria Silva Moura
#e-mail: liviam.moura@gmail.com
# PDSE Scholarship UC Berkeley - PI's BR/USA: Setubal,J /Banfield, J
##########################################################################

import os,subprocess,re,sys,io,operator,itertools
import argparse
import pysam as ps
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqIO.FastaIO import SimpleFastaParser

def reference_lengths(fasta,ext):
    #Get lengths of all reference sequences and return a dict'
    reference_lengths = {}
    for header,seq in SimpleFastaParser(open(fasta+ext)):
        id_,length = header.split(' ')[0], len(seq)
        reference_lengths[id_] = length
        return reference_lengths

def parse_bam(bam_file):
    #Parse bam file and return a dict of start,end of each read on reference'
    bam_parsed = ps.AlignmentFile(bam_file, 'rb')
    bam_dict = defaultdict(list)
    for read in bam_parsed:
        reference = read.reference_name
        if read.is_proper_pair and read.is_paired and reference == read.next_reference_name:
            read_name,reference = read.query_name,reference.split(' ')[0]
            reference_positions = read.get_reference_positions()
            start_pos,end_pos = reference_positions[0],reference_positions[-1]
            info = (read_name,reference,start_pos,end_pos)
            bam_dict[reference].append(info)
    return bam_dict

def find_regions(bam_dict, look_len, reference_lengths):
    #Find areas in references with no support from at least one full read.
    #Iterating through each base on a reference sequence and checking that there are bases upstream and downstream of the position from at least one read'''
    reference_to_bad_regions = defaultdict(set)
    for reference,infos in bam_dict.items():
        good_positions = set()
        reference_length = reference_lengths[reference]
        reference_position_to_reads,read_to_positions = defaultdict(set),defaultdict(set)
        for info in infos:
            read,reference,start_position,end_position = info
            for position in range(start_position,end_position):
                reference_position_to_reads[position].add(read)
                read_to_positions[read].add(position)
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
        for k,g in itertools.groupby(enumerate(bad_positions), lambda x: x[1] - x[0]):
            bad_region = tuple(map(operator.itemgetter(1), g))
            bad_start,bad_end = bad_region[0],bad_region[-1]
            bad_start_end = (bad_start,bad_end)
            reference_to_bad_regions[reference].add(bad_start_end)
    return reference_to_bad_regions

def build_N(ext,num_mm_filt,fasta,thread,bam_filt=""):
    
    ref_bps = reference_lengths(fasta,ext)
    bam_dict = parse_bam(bam_filt+'_sorted.bam')
    dict_replace_chunked = find_regions(bam_dict,5,ref_bps)    

    print(dict_replace_chunked)

path = "/data7/ThiocyanateBioreactor/2018/curator/Livia/genomas/sample_Christ/SR2-18-Sp2"
os.chdir(path)
fasta_in = "SR2-18-Sp2_coassembly_Potentially_Complete_Berkelbacteria_51_5.contigs.fa"
fasta,ext = os.path.splitext(fasta_in)
read1_in = "SR2-18-Sp2_trim_clean.PE.1.fastq"
read2_in = "SR2-18-Sp2_trim_clean.PE.2.fastq"
thread = 6
num_mm = 1
num_mm_filt = 4
minid_in = 0.99
x_times = 6
#alele_freq = 0.6
dp_cov=1
ext_multifasta=True

av_readlen = 150
#aligner(thread,ext,fasta,read1=read1_in,read2=read2_in,bam_out=fasta,minid=minid_in,rescue_mm=round(av_readlen/3))
#filtering_bam(thread,ext,num_mm,num_mm_filt,fasta,read1=read1_in,read2=read2_in)
target = build_N(ext,num_mm_filt,fasta,thread,bam_filt=fasta+"_filtered")