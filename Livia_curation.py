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

def count_times(x_times):
    n = 1
    while n <= x_times:        
        yield n
        n += 1

def remove_N(fasta_semifinal):
    fasta_temp = io.StringIO(str(fasta_semifinal))
    new_fasta =""
    
    for i, seq_record in enumerate(SeqIO.parse(fasta_temp,'fasta')):
        seq_mutable = seq_record.seq.tomutable()
        len_seq = len(seq_mutable)
        for number in range (len_seq-(av_readlen*3),len_seq):
            if seq_mutable[number] == "N":
                remov_end= number
                break
            
        for number in range ((av_readlen*3)-1,0,-1):
            if seq_mutable[number] == "N":
                remov_start= number+1 ### +1?
                break
                
        if i==0:
            new_fasta += '>'+seq_record.id
            new_fasta += '\n'+seq_mutable[remov_start:remov_end]
        else:
            new_fasta += '\n>'+seq_record.id
            new_fasta += '\n'+seq_mutable[remov_start:remov_end] 
                 
    fasta_temp.close()
    return new_fasta

def check_npos(one_fasta):
    ### Personal remember - genome pos (1..N), list pos (0..N)
    ###taking all N index 
    indexes = [index for index,base in enumerate(one_fasta) if base == "N"]
    
    if indexes:
        ###tracking all pos btw the first and last value from a gap 
        check_new=False
        temp_pos= 0
        temp_list=list()
        for index,pos in enumerate(indexes):
            if index == 0:
                pos_temp = pos
                continue                
            else:
                if (pos == pos_temp + 1) :                    
                    temp_list.append(pos)
                    check_new=False
                else:           
                    #print(pos_temp)
                    if not check_new:
                        if index ==1:
                            temp_list.append(pos_temp) 
                        else:
                            check_new=True
                           
                            del temp_list[-1]
                    else:                       
                        temp_list.append(pos_temp)
                        check_new=False
                pos_temp = pos
        del temp_list[-1]
        
        #print(indexes, "Antes")
       
        ###removing from the main N index, the pos in btw
        for item in temp_list:
             indexes.remove(item)
        
        #print(indexes, " DEPOIS")
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#         ### merging regions if they are too close (80% of readlength)
        temp_list=list()
        for i in range(1,len(indexes)-1,2):
            if indexes[i+1] - indexes[i] < (av_readlen*0.8):
                temp_list.extend([i,i+1])
        
        for item in sorted(temp_list, reverse=True):
            del indexes[item]           
        
#         #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        ###Final N_pos - list with startN, end N, Number of N in total
        N_pos = list()
        N_target = list()
        
#         print("FINAL")
#         print(indexes)
        for i in range(0, len(indexes), 2):
            N_pos.append(tuple((indexes[i],indexes[i+1],indexes[i+1]-indexes[i]+1))) ## return the python index mode
            N_target.append(list([indexes[i]+1,indexes[i+1]+1]))  ### return the position for calculate new target      
        return N_pos,N_target
    else:
        return list(),list()

def aligner(thread,ext,fasta="",bam_out="",read1="",read2="",minid=0.99,rescue_mm="",semi=False,pair=False):
    subprocess.call(['bbmap.sh', 'ref='+str(fasta)+str(ext),'overwrite=t'], shell=False)
    if not semi:
        if pair:
            print("pair)")
            subprocess.call(['bbmap.sh', 'ref='+str(fasta)+str(ext), 'in1='+read1, 'in2='+read2, 
                  'outm='+bam_out+'.bam', 'threads='+str(thread),
                  'minid='+str(minid), 'rescuemismatches='+str(rescue_mm),'ambiguous=random','pairedonly=t','overwrite=t'], shell=False)
        else:
            print("notpair")
            subprocess.call(['bbmap.sh', 'ref='+str(fasta)+str(ext), 'in1='+read1, 'in2='+read2, 
                  'outm='+bam_out+'.bam', 'threads='+str(thread),
                  'minid='+str(minid), 'rescuemismatches='+str(rescue_mm),'ambiguous=random','overwrite=t'], shell=False)
    else:
        print("semi")
        subprocess.call(['bbmap.sh', 'ref='+str(fasta)+str(ext), 'in1='+read1, 'in2='+read2, 
                  'outm='+bam_out+'.bam', 'threads='+str(thread),'ambiguous=random','semiperfectmode=t','overwrite=t'], shell=False)
    subprocess.call(['samtools', 'sort', bam_out+'.bam', '-o', bam_out+'_sorted.bam','--threads',str(thread)], shell=False)
    subprocess.call(['samtools', 'index', fasta+'_sorted.bam','-@',str(thread)], shell=False)
    os.remove(bam_out+'.bam')

def filtering_bam(thread,ext,num_mm,num_mm_filt,fasta="",read1="",read2=""):
    ### Filtering reads with less than 0 or [1] mismatch (default to be decided) 
    samfile = ps.AlignmentFile(fasta+'_sorted.bam',"rb")
    match_reads = list()
    handler=""
    for read in samfile.fetch():
        if not read.has_tag('NM'):
            continue
        if not (read.get_tag('NM') > num_mm): ## if TAG NM is not greater than num_mm (default = 1) -> read_list to be filtered
            match_reads.append(read.qname)
    samfile.close()

    ### creating a file with unique reads names that fulfill the criteria above 
    match_reads = list(set(match_reads)) 
    with open ('matched_reads', 'w') as outfile:
        outfile.write("\n".join(str(item) for item in match_reads))
    outfile.close()
    
    ### Filtering original fastq with the high stringency reads 
    subprocess.call(['filterbyname.sh', 'in='+read1, 'in2='+read2, 'out=res_R1.fastq', 'out2=res_R2.fastq', 'names=matched_reads', 'include=t','overwrite=t'])
    
    aligner(thread,ext,fasta,read1="res_R1.fastq",read2="res_R2.fastq",bam_out=fasta+"_filtered",minid=0.8,rescue_mm=round(av_readlen/2))
    #aligner(thread,ext,fasta,read1="res_R1.fastq",read2="res_R2.fastq",bam_out=fasta+"_semiperf",semi=True)
    
def reference_lengths(fasta,ext):
    length_dict = {}
    for header,seq in SimpleFastaParser(open(fasta+ext)):
        id_,length = header.split(' ')[0], len(seq)
        length_dict[id_] = length
    return length_dict

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
    reference_to_bad_regions = defaultdict(list)
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
        bad_positions = sorted(bad_positions)
        print(bad_positions)
        print("BAD_pos")
        #for k,g in itertools.groupby(enumerate(bad_positions), lambda x: x[1] - x[0]):
        #    bad_region = tuple(map(operator.itemgetter(1), g))
        #    bad_start,bad_end = bad_region[0],bad_region[-1]
        #    bad_start_end = (bad_start,bad_end)
        #    reference_to_bad_regions[reference].add(bad_start_end)
        reference_to_bad_regions[reference] = bad_positions
    return reference_to_bad_regions

def build_N(ext,num_mm_filt,fasta,thread,bam_filt=""):
    if os.path.exists("target"): ### remove target file if it exists
        os.remove("target") 
        
    ref_bps = reference_lengths(fasta,ext)
    bam_dict = parse_bam(bam_filt+'_sorted.bam')
    dict_replace = find_regions(bam_dict,5,ref_bps)    
#     print(dict_replace)
#     print("OPAAA")
    
    for key, value in dict_replace.items():        
        if len(value) < 3:
            continue
        check_new=False
        temp_pos= 0
        temp_list=list()
        for index,pos in enumerate(dict_replace[key]):            
            if index == 0:
                pos_temp = pos
                continue                
            else:                              
                if (pos == pos_temp + 1) :                    
                    temp_list.append(pos)
                    check_new=False
                else:
                    if not check_new:
                        if index == 1:
                            temp_list.append(pos_temp) 
                        else:
                            check_new=True
                            del temp_list[-1]                        
                    else:                       
                        temp_list.append(pos_temp)
                        check_new=False               
            pos_temp = pos        
        del temp_list[-1]        
                
        for item in temp_list:
            dict_replace[key].remove(item)
    
#     print(dict_replace)
#     print("final")
    ### merging regions if they are too close (60% of readlength)
    temp_list=list()
    for i in range(1,len(dict_replace[key])-1,2):
        if dict_replace[key][i+1] - dict_replace[key][i] < (av_readlen*0.6):
            temp_list.extend([i,i+1])

    for item in sorted(temp_list, reverse=True):
        del dict_replace[key][item]           

    dict_replace_chunked = defaultdict(list)
    for i in range(0, len(dict_replace[key]), 2):
        dict_replace_chunked[key].append(tuple(dict_replace[key][i:i+2]))
            
    #print(dict_replace_chunked)  
    #print("Chunked")
     # #----------------Replace region for N---------------------------

    new_fasta =''
    dict_error_pos = defaultdict(list)
    dict_only_errors = defaultdict(list)
    
    for i,seq_record in enumerate(SeqIO.parse(fasta+ext,'fasta')):
        #N_pos_init,N_target_init = check_npos(seq_record.seq) ###Check N pos before insert new N's
        if seq_record.id in dict_replace_chunked.keys():
            
            seq_mutable = seq_record.seq.tomutable()
            count=0
            
            for j,(start,end) in enumerate(dict_replace_chunked[seq_record.id]):   ## Add N after the end position
                seq_mutable[start-1+count:end+count] = (end-start+1)*"N"
                dict_error_pos[seq_record.id].append(list([start+(av_readlen*3)*(j+1),end+(av_readlen*3)*(j+2)]))
                for i in range(av_readlen*3):
                    seq_mutable.insert(end+count,"N")
                count+=av_readlen*3
                
                if seq_record.seq[round( (start+end)/2 )] != "N":
                    dict_only_errors[seq_record.id].append(tuple([start,end]))    
            
            ### This is part of extension script - It adds av_readlen*4 "N"s in the beginning and end from a file
            for i in range(av_readlen*3):
                seq_mutable.insert(0,"N")
                seq_mutable.append("N")
                
            to_check_n = seq_mutable    
            if i==0:                
                new_fasta += '>'+seq_record.id
                new_fasta += '\n'+seq_mutable
            else:                
                new_fasta += '\n>'+seq_record.id
                new_fasta += '\n'+seq_mutable  
            
        else:
            ### This is part of extension script - It adds av_readlen*3 "N"s in the beginning and ending from a file
            seq_mutable = seq_record.seq.tomutable()
            for i in range(av_readlen*3):
                seq_mutable.insert(0,"N")
                seq_mutable.append("N")
            
            to_check_n = seq_mutable
            if i==0:                
                new_fasta += '>'+seq_record.id
                new_fasta += '\n'+seq_mutable
                
            else:                
                new_fasta += '\n>'+seq_record.id
                new_fasta += '\n'+seq_mutable
                
        N_pos_final,N_target = check_npos(to_check_n)         
        dict_error_pos[seq_record.id].append(list([1,av_readlen*3]))
        dict_error_pos[seq_record.id].append(list([len(seq_mutable)+1 - (av_readlen*3),len(seq_mutable)]))
        dict_error_pos[seq_record.id] = sorted(dict_error_pos[seq_record.id], key=lambda i: i[0])
        
        with open('target', 'a') as target:
            for start,end in N_target:
                target.write(seq_record.id+"\t"+str(start)+"\t"+str(end)+"\n")              
    
    target.close()
    
    with open('v_0'+ext, 'w') as fasta_N:
        fasta_N.write(str(new_fasta))
        fasta_N.close()
    
    
    with open('Error_loc_orig','w') as error_loc:
        for key,value in dict_only_errors.items():
            error_loc.write(str(key)+''+str(value)+"\n")
        error_loc.close()                    
    
    #print(dict_only_errors)
    return dict_error_pos

def check_overlap(new_fasta):
    fasta_temp = io.StringIO(str(new_fasta))
    new_fasta =""
    
    if os.path.exists("target"): ### remove target file if it exists
        os.remove("target")
    
    for i, seq_record in enumerate(SeqIO.parse(fasta_temp,'fasta')):
        new_target=list()
        N_pos,N_target_temp = check_npos(seq_record.seq)               
        
        print(N_pos, seq_record.id, "NPOS and ID")
        N_pos = N_pos[1:-1] ### Remove the edges N (they dont need to have overlap checked)
        new_subfasta = ""
        temp_dif_pos=0
        if not N_pos:
            new_subfasta = seq_record.seq
            
        control_index =  list()    
        for j,(start,end,number) in enumerate(N_pos):
            #print(start,end,"START|END")
            if j==0: ### first time doesnt have index modification
                temp_left = seq_record[0:start].seq  ### start = N pos; but here would be Last base pos
                temp_right = seq_record[end+1:-1].seq ### end = N pos; but here would be first base pos
                comp_left = temp_left[-(av_readlen*2):]
                comp_right = temp_right[0:12]
                
                ## Self remember alignment list 0 to N
                alignments = pairwise2.align.localms(comp_left,comp_right,2,-1,-.5,-.1)
                                
                if alignments[0][2] == 24:                    
                    print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######DEU MATCH")                    
                    variat=(av_readlen*2)-alignments[0][3]
                    #print(number,variat)
                    new_subfasta = temp_left[:-variat]+temp_right
                    temp_dif_pos+=number+variat
                    control_index.append(tuple([j+1,N_pos[j][2]+variat]))
                    
####### MAYBE WE'LL HAVE TO CREATE A REVERSE CHECK BECAUSE OF THE SIZE OF READS THAT SOMETIMES LINK TO THE EDGES #################
                else:                    
                    comp_left = temp_left[-12:]
                    comp_right = temp_right[:(av_readlen*2)]
                    alignments = pairwise2.align.localms(comp_right,comp_left,2,-1,-.5,-.1)
                    
                    if alignments[0][2] == 24:
                        print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######ESPECIAL1")                        
                        new_subfasta = temp_left+temp_right[alignments[0][4]:]
                        temp_dif_pos+=number+alignments[0][4]
                        control_index.append(tuple([j+1,N_pos[j][2]+alignments[0][4]]))
                    else:                   
                        new_subfasta = seq_record.seq
###################################################################################################################
            else:                
                temp_left = new_subfasta[0:start-temp_dif_pos]  ### start = N pos; but here would be Last base pos
                temp_right = new_subfasta[end+1-temp_dif_pos:-1] ### end = N pos; but here would be first base pos
                comp_left = temp_left[-(av_readlen*2):]
                comp_right = temp_right[0:12]
               
                ## Self remember alignment list 0-N
                alignments = pairwise2.align.localms(comp_left,comp_right,2,-1,-.5,-.1)
                                
                if alignments[0][2] >= 24:
                    print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######DEU MATCH")
                    variat=(av_readlen*2)-alignments[0][3]  
                    #print(number,variat)
                    new_subfasta = temp_left[:-variat]+temp_right
                    temp_dif_pos+=number+variat
                    control_index.append(tuple([j+1,N_pos[j][2]+variat]))
                ####### MAYBE WE'LL HAVE TO CREATE A REVERSE CHECK BECAUSE OF THE SIZE OF READS THAT SOMETIMES LINK TO THE EDGES #################
                else:                    
                    comp_left = temp_left[-12:]
                    comp_right = temp_right[:(av_readlen*2)]
                    alignments = pairwise2.align.localms(comp_right,comp_left,2,-1,-.5,-.1)
                                        
                    if alignments[0][2] == 24:                        
                        print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######ESPECIAL2")
                        new_subfasta = temp_left+temp_right[alignments[0][4]:]
                        temp_dif_pos+=number+alignments[0][4]
                        control_index.append(tuple([j+1,N_pos[j][2]+alignments[0][4]]))                    
###################################################################################################################
        
        #### Creating new target based on control and N_target_temp ######
#         print("Control_index --- O que houve aqui, deOs mio?")
#         print(control_index)
        
#         print(N_target_temp, "N TARGET TEMP")
        
        if control_index:
            for index,number in control_index:
                for i in range(0,index):
                    new_target.append(list([N_target_temp[i][0],N_target_temp[i][1]]))
                for i in range(index,len(N_target_temp)):
                    new_target.append(list([N_target_temp[i][0]-number,N_target_temp[i][1]-number]))
        else:
            new_target = N_target_temp
            
        with open('target', 'a') as file_target:   
            for start,end in new_target:
                file_target.write(seq_record.id+"\t"+str(start)+"\t"+str(end)+"\n")
                   
#         print(new_target, " New_target to be stocked")
#         print(" ")
        if i==0:
            new_fasta += '>'+seq_record.id
            new_fasta += '\n'+new_subfasta
        else:
            new_fasta += '\n>'+seq_record.id
            new_fasta += '\n'+new_subfasta  
        
        N_pos_final,N_target_temp_final = check_npos(new_subfasta)        
    file_target.close()
    fasta_temp.close()
    return new_fasta

def var_cal_fix(count,ext,thread,x_times,dp_cov,new_target):    
    aligner(thread,ext,fasta="v_"+str(count-1),bam_out="v_"+str(count-1),read1="res_R1.fastq",read2="res_R2.fastq",semi=True)

###################################### Version with alele freq ###################
#     varfind1 = subprocess.Popen("bcftools mpileup -R target -A -B -f v_{}{} v_{}_sorted.bam".format(count-1,ext,count-1), 
#                         shell=True, stdout=subprocess.PIPE,universal_newlines=True)
#     varfind2 = subprocess.Popen("bcftools call --ploidy 1 -Mcv", 
#                         shell=True, stdin=varfind1.stdout,stdout=subprocess.PIPE,universal_newlines=True)
    
#     varfind1.stdout.close()
#     output = varfind2.communicate()[0]
#     varfind2.stdout.close()
###################################################################################

    varfind1 = subprocess.Popen("bcftools mpileup -R target -A -B -f v_{}{} v_{}_sorted.bam".format(count-1,ext,count-1), 
                        shell=True, stdout=subprocess.PIPE,universal_newlines=True)
    output = varfind1.communicate()[0]
    varfind1.stdout.close()

    
    
    ### Select the variants with af >= alele_freq and store it on a list
    var_list= dict()
    for line in output.splitlines():
        if re.match('^#',line):
            continue
        line_list= (line.split('\t')) #0=header;1=pos;3=ref;4=var("," slipted);7=flags(;/= splited)
            #print(line_list)
        if line_list[0] in var_list.keys():
            if int(re.search('DP=(.+?);',line_list[7]).group(1)) >= dp_cov:
                if ((line_list[4].split(",")[0] == 'A') or
                (line_list[4].split(",")[0] == 'C') or
                (line_list[4].split(",")[0] == 'T') or 
                (line_list[4].split(",")[0] == 'G')):
                    var_list[line_list[0]].append((int(line_list[1])-1,len(line_list[3]),len(line_list[4].split(",")[0]),line_list[4].split(",")[0]))
        else:
            if int(re.search('DP=(.+?);',line_list[7]).group(1)) >= dp_cov:
                if ((line_list[4].split(",")[0] == 'A') or
                (line_list[4].split(",")[0] == 'C') or
                (line_list[4].split(",")[0] == 'T') or 
                (line_list[4].split(",")[0] == 'G')):
                    var_list[line_list[0]] = [(int(line_list[1])-1,len(line_list[3]),len(line_list[4].split(",")[0]),line_list[4].split(",")[0])]

###################################### Version with alele freq ###################                    
#         line_list= (line.split('\t')) #0=header;1=pos;3=ref;4=var("," slipted);7=flags(;/= splited)
#         #print(line_list)
#         if line_list[0] in var_list.keys():
#             if float(re.search('AF1=(.+?);',line_list[7]).group(1)) >= alele_freq:
#                 var_list[line_list[0]].append((int(line_list[1])-1,len(line_list[3]),len(line_list[4].split(",")[0]),line_list[4].split(",")[0]))
#         else:
#             if float(re.search('AF1=(.+?);',line_list[7]).group(1)) >= alele_freq:
#                 var_list[line_list[0]] = [(int(line_list[1])-1,len(line_list[3]),len(line_list[4].split(",")[0]),line_list[4].split(",")[0])]
###################################################################################    
  
    new_fasta =''
    for i,seq_record in enumerate(SeqIO.parse("v_"+str(count-1)+ext,'fasta')):
        pos_temp=0
        if seq_record.id in var_list.keys():
            seq_mutable = seq_record.seq.tomutable()
            for value in var_list[seq_record.id]:
                pos,len_ref,len_var,var = value
                seq_mutable[pos+pos_temp:pos+pos_temp+len_ref] = var
                pos_temp+=(len_var-len_ref)
            if i==0:
                new_fasta += '>'+seq_record.id
                new_fasta += '\n'+seq_mutable
            else:
                new_fasta += '\n>'+seq_record.id
                new_fasta += '\n'+seq_mutable
        else:
            if i==0:
                new_fasta += '>'+seq_record.id
                new_fasta += '\n'+seq_record.seq
            else:
                new_fasta += '\n>'+seq_record.id
                new_fasta += '\n'+seq_record.seq

    with open('teste_pre_v_'+str(count)+ext, 'w') as handle:
        handle.write(str(new_fasta))
        handle.close()       
        
    new_fasta2 = check_overlap(new_fasta)
    if count != x_times:
        with open('v_'+str(count)+ext, 'w') as handle:
            handle.write(str(new_fasta2))
            handle.close()        
    else:
        with open('curated_LMSM'+ext, 'w') as handle:
            handle.write(str(new_fasta2))
            handle.close()

%%time 
path = "/data7/ThiocyanateBioreactor/2018/curator/Livia/genomas/QS_7_Complete_huge_phage_40_25"
os.chdir(path)
fasta_in = "QS_7_Complete_huge_phage_40_25.fa"
fasta,ext = os.path.splitext(fasta_in)
read1_in = "qs_7.PE.1.fastq.gz"
read2_in = "qs_7.PE.2.fastq.gz"
thread = 6
num_mm = 2
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
# print(target,"FIRST")
# print(" "),
# number_times = count_times(x_times)

# for count in number_times :
#     print(count)
#     var_cal_fix(count,ext,thread,x_times,dp_cov,target)