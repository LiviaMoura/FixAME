from collections import defaultdict
from Bio.SeqUtils import GC
from xopen import xopen
from Bio import SeqIO
import pandas as pd
import pysam as ps

allowed_mismatches = 0

def partial_stats(bam_file, num_mm, threads):
        
    bam_parsed = ps.AlignmentFile(bam_file, 'rb', threads = threads)
    
    template_lengths = []

    read_count = 0
    
    for read in bam_parsed:
        if read_count == 100000:
            break
        
        if read.is_proper_pair and read.get_tag('NM') <= num_mm:
            
                reference = read.reference_name
                mate_reference = read.next_reference_name

                template_length = read.template_length
                template_length = abs(template_length)
                template_lengths.append(template_length)

                read_count += 1 

    template_lengths = pd.Series(template_lengths)

    average_template_length = template_lengths.mean()
    var_template_length = template_lengths.var()
    
    min_template_length = average_template_length - var_template_length
    max_template_length = average_template_length + var_template_length
    
    bam_parsed.close()
    
    return average_template_length, var_template_length, min_template_length, max_template_length
     
average_template_length, var_template_length, min_template_length, max_template_length = partial_stats(bam_file, allowed_mismatches, threads)
 
def alternate_references(bam_file, num_mm, threads, max_template_length):
    
    bam_parsed = ps.AlignmentFile(bam_file, 'rb', threads = threads)    
    
    max_template_length = int(max_template_length)
    
    pair_dict = defaultdict(lambda: defaultdict(list))
    
    count = 0
    
    for read in bam_parsed:
        
        reference = read.reference_name
        mate_reference = read.next_reference_name
        
#         if count == 100000:
#             break
        
        if reference != mate_reference and read.get_tag('NM') == 0: #and not read.mate_is_unmapped:
            reference = reference.split()[0]
            mate_reference = mate_reference.split()[0]

            read_name = read.query_name
            read_base_name = read_name.split()[0]

            reference_length = length_dict[reference]

            length_from_end = reference_length - max_template_length
            length_from_end = int(length_from_end)
            
            query_length,template_length = read.query_alignment_length, read.template_length
            reference_positions = read.get_reference_positions(full_length=False)
            start_pos = reference_positions[0] + 1
            end_pos = start_pos + query_length

            if end_pos <= max_template_length or end_pos >= length_from_end:
            #if end_pos:
                read_reverse = read.is_reverse
                
                info = (reference, start_pos, end_pos, read_reverse)
                
                if read.is_read1:
                    pair_dict[read_base_name][1].append((info))
                
                elif read.is_read2:
                    pair_dict[read_base_name][2].append((info)) 

                count += 1
                    
    multi_reference_to_pair_dict = {}
    
    for read_base_name,pair_info_dict in pair_dict.items():
        
        read_1_info = pair_info_dict[1]
        read_2_info = pair_info_dict[2]
        
        read_1_info_length = len(read_1_info)
        read_2_info_length = len(read_2_info)
        
        references = set()
        
        references_1 = set()
        references_2 = set()
                        
        pair_info_dict = defaultdict(set)

        if read_1_info_length > 0 and read_2_info_length > 0:
            
            for info in read_1_info:
                reference, start_pos, end_pos, read_reverse = info
                print(end_pos)
                references.add(reference)
                references_1.add(reference)
                out_info = (reference, read_reverse)
                pair_info_dict[1].add(out_info)
            
            for info in read_2_info:
                reference, start_pos, end_pos, read_reverse = info
                print(end_pos)
                references.add(reference)
                references_2.add(reference)
                out_info = (reference, read_reverse)
                pair_info_dict[2].add(out_info)
                
        references = tuple(sorted(references))
                
        if len(references) >= 2:
            multi_reference_to_pair_dict[references] = pair_info_dict
            
    bam_parsed.close()
        
    return multi_reference_to_pair_dict
    
multi_reference_to_pair_dict = alternate_references(bam_file, allowed_mismatches, threads, max_template_length)

 
