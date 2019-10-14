import argparse
import os
import logging
import sys
import re
import subprocess
import fixame
import pysam as ps
from datetime import datetime
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqIO.FastaIO import SimpleFastaParser

from fixame.aligner import aligner
from fixame.lm_common import *
from fixame.lm_error_finder import Error_finder

from xopen import xopen

__author__ = "Livia Moura"
__copyright__ = "Copyright 2019"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"

def lm_curation_validate(**kwargs):
    '''
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
        ESCREVER ALGO AQUI.
    '''

    method = common_validate(**kwargs)
   
    if (kwargs.get('minid') < 0.76) or (kwargs.get('minid') > 1.00):
        logging.info('Checking minimum identify for the first alignment')
        logging.error("Please, verify if -minid is >= 0.76 and <= 1.00 ")
        sys.exit()

    if (kwargs.get('num_mismatch') < 0) or (kwargs.get('num_mismatch') > 5):
        logging.info('Checking number of mismatches allowed for the filtering reads')
        logging.error("-num_mismatch 0>=x>=5")
        sys.exit()
 
    output_dir = os.path.abspath(kwargs.get('output_dir'))
    
    av_readlen = average_read(kwargs.get('r1'))
     
    
    if method == 0:
        fasta_in = os.path.realpath(os.path.expanduser(kwargs.get('fasta')))
        if kwargs.get('r12'):
            read12_in = os.path.realpath(os.path.expanduser(kwargs.get('r12')))
            read1_in =''
            read2_in = ''
        else:    
            read1_in = os.path.realpath(os.path.expanduser(kwargs.get('r1')))
            read2_in = os.path.realpath(os.path.expanduser(kwargs.get('r2')))
            read12_in = ''
        name_fasta = os.path.splitext(os.path.basename(fasta_in))[0]
        
        #Creating folder fixame
        mydir = os.path.join(output_dir,'fixame_'+datetime.now().strftime('%Y-%b-%d_%H-%M-%S'))
        os.mkdir(mydir)
        
        #Creating folder for the new fasta without no reason Ns
        os.mkdir(os.path.join(mydir,'new_fastas'))
        check_overlap(mydir,fasta_in,av_readlen,True)
        
        #Creating folder for temp files - alignments related
        os.mkdir(os.path.join(mydir,'tmp'))
        aligner(mydir,kwargs.get('threads'),kwargs.get('minid'),mydir+'/new_fastas/'+name_fasta+'_renewed.fasta',r1=read1_in,r2=read2_in,r12=read12_in, bam_out=name_fasta+'_renewed')
        
        #Filtering the fastq - Make curation process faster
        filtering_bam(mydir,kwargs.get('threads'),kwargs.get('num_mismatch'),mydir+'/tmp/'+name_fasta+'_renewed',read1_in,read2_in,read12_in)
        
        orig_target,fasta_len,coverage_dict,av_readlen,mean_gap_length = build_N(mydir,kwargs.get('threads'),mydir+'/new_fastas/'+name_fasta+'_renewed.fasta',mydir+'/tmp/'+name_fasta+'_renewed',kwargs.get('fasta_cov'))
        
        for count,r in enumerate(range(kwargs.get('xtimes')),1):
            print(count)
            var_cal_fix(mydir,count,kwargs.get('threads'),kwargs.get('xtimes'),kwargs.get('dp_cov'),orig_target,fasta_len,av_readlen)
        
        os.mkdir(os.path.join(mydir,'curated_seqs'))
        remove_N(mydir,name_fasta,os.path.join(mydir,'tmp','v_'+str(kwargs.get('xtimes'))+'.fasta'),orig_target,fasta_len,av_readlen,mean_gap_length)

    #else:
    #    fasta_array = []
    #    for sample in os.listdir(os.path.realpath(os.path.expanduser(kwargs.get('output_dir')))):	
    #        if sample.lower().endswith(".fasta") or sample.lower().endswith(".fa") or sample.lower().endswith(".fna"):
    #	        fasta_array.append(sample)
    #    for fasta_in in fasta_array:
    #        check_overlap(output_dir,fasta_in,av_readlen,True)
         
           
#    os.path.realpath(os.path.expanduser)

def average_read(r1_fastq):
    '''Initial fast avg read length - It'll be recalculated with more precision later'''
    sum_read_len=0
    for num,record in enumerate(SeqIO.parse(xopen(r1_fastq), 'fastq')):
        sum_read_len += len(record.seq)
        if num == 999:
            av_read_len = int(sum_read_len/1000)
            break
    return av_read_len

def check_overlap(output_dir,fasta,av_readlen,user_file=False,count=''):
    '''Check overlap from the border of N's regions'''
    if user_file:
        name_fasta = os.path.splitext(os.path.basename(fasta))[0]
        new_fasta_temp = open(os.path.join(output_dir,'new_fastas',name_fasta+'_renewed.fasta'),'w')
    else:
        if os.path.exists(os.path.join(output_dir,'tmp','target')): ### remove target file if it exists
            os.remove(os.path.join(output_dir,'tmp','target')) 
        new_fasta_temp = open(os.path.join(output_dir,'tmp','v_'+str(count)+'.fasta'),'w')
        
    for seq_record in SeqIO.parse(fasta,'fasta'):
        new_target=list()
        N_pos,N_target_temp = check_npos(seq_record.seq,av_readlen)               
        if not user_file:
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
                if comp_right in comp_left:
                    alignments = pairwise2.align.localms(comp_left,comp_right,2,-1,-.5,-.1)
                    print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######DEU MATCH")                    
                    variat=(av_readlen*2)-alignments[0][3]                    
                    new_subfasta = temp_left[:-variat]+temp_right
                    temp_dif_pos+=number+variat
                    control_index.append(tuple([j+1,N_pos[j][2]+variat]))
                    
            ####### REVERSE CHECK #################
                else:                    
                    comp_left = temp_left[-12:]
                    comp_right = temp_right[:(av_readlen*2)]
                    
                    if comp_left in comp_right:
                        alignments = pairwise2.align.localms(comp_right,comp_left,2,-1,-.5,-.1)
                        print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######ESPECIAL1")                        
                        new_subfasta = temp_left+temp_right[alignments[0][4]:]
                        temp_dif_pos+=number+alignments[0][4]
                        control_index.append(tuple([j+1,N_pos[j][2]+alignments[0][4]]))
                    else:                   
                        new_subfasta = seq_record.seq
            else:    
                temp_left = new_subfasta[0:start-temp_dif_pos]  ### start = N pos; but here would be Last base pos
                temp_right = new_subfasta[end+1-temp_dif_pos:-1] ### end = N pos; but here would be first base pos
                comp_left = temp_left[-(av_readlen*2):]
                comp_right = temp_right[0:12]
               
                ## Self remember alignment list 0-N
                if comp_right in comp_left:
                    alignments = pairwise2.align.localms(comp_left,comp_right,2,-1,-.5,-.1)
                    print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######DEU MATCH")
                    variat=(av_readlen*2)-alignments[0][3]                      
                    new_subfasta = temp_left[:-variat]+temp_right
                    temp_dif_pos+=number+variat
                    control_index.append(tuple([j+1,N_pos[j][2]+variat]))
                #######  REVERSE CHECK  #################
                else:                    
                    comp_left = temp_left[-12:]
                    comp_right = temp_right[:(av_readlen*2)]
                    
                    if comp_left in comp_right:
                        alignments = pairwise2.align.localms(comp_right,comp_left,2,-1,-.5,-.1)
                        print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######ESPECIAL2")
                        new_subfasta = temp_left+temp_right[alignments[0][4]:]
                        temp_dif_pos+=number+alignments[0][4]
                        control_index.append(tuple([j+1,N_pos[j][2]+alignments[0][4]]))                    
               
        #### Creating new target based on control and N_target_temp ######
        if control_index:
            for index,number in control_index:
                for i in range(0,index):
                    new_target.append(list([N_target_temp[i][0],N_target_temp[i][1]]))
                for i in range(index,len(N_target_temp)):
                    new_target.append(list([N_target_temp[i][0]-number,N_target_temp[i][1]-number]))
        else:
            new_target = N_target_temp
            
        if not user_file:
            with open(os.path.join(output_dir,'tmp','target'), 'a') as file_target:   
                for start,end in new_target:
                    file_target.write(seq_record.id+"\t"+str(start)+"\t"+str(end)+"\n")

        new_fasta_temp.write(">{}\n{}\n".format(seq_record.id,new_subfasta))
        
    new_fasta_temp.close()
    if not user_file:
        file_target.close()    

   
def check_npos(one_fasta,av_readlen,error=False):
    '''Check Npos from contigs'''
    ### Personal remember - genome pos (1..N), list pos (0..N)
    ###taking all N index 
    lst = list(one_fasta)
    indexes=[]
    start=0    
    while True:
        try:
            start = lst.index("N", start)
            indexes.append(start)
            start += 1
        except ValueError: 
            break        
    
    check_new=False
    temp_pos= 0
    final_list=list()
    
    ###tracking all pos btw the first and last value from a gap (individuals N as well)
    if indexes:     
        for index,pos in enumerate(indexes):
            if index == 0:
                final_list.append(pos)
                pos_temp = pos                
                continue                
            else:                
                if (pos == pos_temp + 1):
                    check_new=False
                else:
                    final_list.append(pos_temp)
                    if not check_new:
                        final_list.append(pos)
                        check_new=True
                    else:                       
                        final_list.append(pos)
                        check_new=False
                pos_temp = pos
                
        final_list.append(pos)
        ### merging regions if they are too close (80% of readlength)
        temp_list=list()
        for i in range(1,len(final_list)-1,2):
            if final_list[i+1] - final_list[i] < av_readlen*0.8:
                temp_list.extend([i,i+1])
        
        for item in sorted(temp_list, reverse=True):
            del final_list[item]           
        
        ###Final N_pos - list with startN, end N, Number of N in total
        N_pos = list()
        N_target = list()

        for i in range(0, len(final_list), 2):
            N_pos.append(tuple((final_list[i],final_list[i+1],final_list[i+1]-final_list[i]+1))) ## return the python index mode
            N_target.append(list([final_list[i]+1,final_list[i+1]+1]))  ### return the position for calculate new target      
        
        return N_pos, N_target
    else:
        return list(),list()

def filtering_bam(output_dir,thread,num_mm,bam_sorted,r1,r2,r12):
    ''' Filtering in reads with less than [2] mismatch.
        Step used to keep only interesting reads to reduce the fastq's size'''
    samfile = ps.AlignmentFile(bam_sorted+'_sorted.bam',"rb")
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
    with open (output_dir+'/tmp/matched_reads', 'w') as outfile:
        outfile.write("\n".join(str(item) for item in match_reads))
    outfile.close()
    
    ### Filtering original fastq with the high stringency reads 
    varpath=script_path()
    subprocess.call([os.path.join(varpath,'tools','bbmap','filterbyname.sh'), 'in='+r1, 'in2='+r2, 'out='+output_dir+'/tmp/res_R1.fastq', 'out2='+output_dir+'/tmp/res_R2.fastq', 'names='+output_dir+'/tmp/matched_reads', 'include=t','overwrite=t'])

def build_N(output_dir,thread,fasta,bam_file,fasta_cov):
    ''' Replace error for N and extend the contig's edges for N - pre-curation step'''
    
    if os.path.exists(os.path.join(output_dir,'tmp',"target")): ### remove target file if it exists
            os.remove(os.path.join(output_dir,'tmp',"target"))

    ref_bps = Error_finder.reference_lengths(fasta)
    bam_dict,ref_read_lengths,mean_template_length = Error_finder.parse_bam(bam_file+'_sorted.bam',1)
    dict_replace,coverage_dict,av_readlen,mean_gap_length = Error_finder.find_regions(fasta_cov,bam_dict,5,ref_bps,ref_read_lengths,mean_template_length)
    dict_replace_chunked = defaultdict(list)     
    
    for key, value in dict_replace.items():        
        check_new=False
        temp_pos= 0
        final_list=defaultdict(list)        
        ###tracking all pos btw the first and last value from a gap (individuals N as well)
        for index,pos in enumerate(dict_replace[key]):
            if index == 0:
                final_list[key].append(pos)
                pos_temp = pos                
                continue                
            else:                
                if (pos == pos_temp + 1):
                    check_new=False
                else:
                    final_list[key].append(pos_temp)
                    if not check_new:
                        final_list[key].append(pos)
                        check_new=True
                    else:                       
                        final_list[key].append(pos)
                        check_new=False
                pos_temp = pos
        
        final_list[key].append(pos)

        ### merging regions if they are too close (80% of readlength)
        temp_list=list()
        for i in range(1,len(final_list[key])-1,2):
            if final_list[key][i+1] - final_list[key][i] < (av_readlen*0.8):
                temp_list.extend([i,i+1])

        for item in sorted(temp_list, reverse=True):
            del final_list[key][item]
        
        for i in range(0, len(final_list[key]), 2):
            dict_replace_chunked[key].append(tuple(final_list[key][i:i+2]))            
        
     # #---------------- Replace region for N ---------------------------
    #new_fasta =''
    dict_error_pos = defaultdict(list)
    dict_only_errors = defaultdict(list)
    dict_len = defaultdict()           
    ext_size = av_readlen*3    
    set_keys = set(dict_replace_chunked.keys())  
    
    fasta_N = open(os.path.join(output_dir,'tmp','v_0.fasta'),'w')
    target = open(os.path.join(output_dir,'tmp','target'), 'w')
    
    for seq_record in SeqIO.parse(fasta,'fasta'):                
        dict_len[seq_record.id] = len(seq_record)        
        seq_mutable = ""
        
        if seq_record.id in set_keys:
            set_keys.remove(seq_record.id)
            seq_mutable = seq_record.seq.tomutable()            
            count=0
            target.write("{}\t{}\t{}\n".format(seq_record.id,"1",ext_size))
            
            for j,(start,end) in enumerate(dict_replace_chunked[seq_record.id]):   ## Add N after the end position
                seq_mutable[start-1+count:end+count] = (end-start+1)*"N"
                dict_error_pos[seq_record.id].append(list([start+(ext_size)*(j+1),end+(ext_size)*(j+2)]))
                seq_mutable[end+count:end+count] = "N" * ext_size                  
                count+=ext_size
                target.write("{}\t{}\t{}\n".format(seq_record.id,start+count,end+count+ext_size))
                
                if seq_record.seq[round( (start+end)/2 )] != "N":
                    dict_only_errors[seq_record.id].append(tuple([start,end]))    

            ### This is part of extension script - It adds av_readlen*3 "N"s in the beginning and end from a file
            seq_mutable[0:0] = "N" * ext_size
            seq_mutable.extend(("N"* ext_size))
            
            len_mut= len(seq_mutable)    
            target.write("{}\t{}\t{}\n".format(seq_record.id,len_mut+1 - ext_size, len_mut))
            fasta_N.write(">{}\n{}\n".format(seq_record.id,seq_mutable))
        else:
            ### This is part of extension script - It adds av_readlen*3 "N"s in the beginning and ending from a file
            seq_mutable = seq_record.seq.tomutable()
            len_mut= len(seq_mutable) 
            seq_mutable[0:0] = "N" * ext_size
            seq_mutable.extend(("N"* ext_size))
            target.write("{}\t{}\t{}\n".format(seq_record.id,"1",ext_size))            
            target.write("{}\t{}\t{}\n".format(seq_record.id,len_mut+1 - ext_size, len_mut))
            fasta_N.write(">{}\n{}\n".format(seq_record.id,seq_mutable))
        
        dict_error_pos[seq_record.id].append(list([1,ext_size]))
        dict_error_pos[seq_record.id].append(list([len(seq_mutable)+1 - (ext_size),len(seq_mutable)]))
        dict_error_pos[seq_record.id] = sorted(dict_error_pos[seq_record.id], key=lambda i: i[0]) 
        
    fasta_N.close()
    target.close()
    with open(os.path.join(output_dir,'Errors_location.txt'),'w+') as error_loc:
        error_loc.write("Sequence\terror_position\n")
        for key,value in dict_only_errors.items():
            for item in value:
                error_loc.write("{}\t{}\n".format(key,item))
        error_loc.close()                    
    
    return dict_replace_chunked,dict_len,coverage_dict,av_readlen,mean_gap_length

def var_cal_fix(output_dir,count,thread,x_times,dp_cov,orig_target,fasta_len,av_readlen):    
    '''Find changes on Ns position and replace them'''
    
    aligner(output_dir,thread,0,os.path.join(output_dir,'tmp','v_'+str(count-1)+'.fasta'),r1=os.path.join(output_dir,'tmp','res_R1.fastq'),r2=os.path.join(output_dir,'tmp','res_R2.fastq'),r12="", bam_out='v_'+str(count-1),semi=True)
    varpath = script_path()
    
    # Run variant finder - bcftools
    print('CHEGUEI NO MPILEUP')
    comm = "{}/bcftools mpileup -R {}/target -A -B -f {}/v_{}.fasta {}/v_{}_sorted.bam".format(os.path.join(varpath,'tools','bcftools-1.9'),os.path.join(output_dir,'tmp'),os.path.join(output_dir,'tmp'),count-1,os.path.join(output_dir,'tmp'),count-1)
    varfind = subprocess.Popen(comm, shell=True, stdout=subprocess.PIPE,universal_newlines=True)
    output = varfind.communicate()[0]
    varfind.stdout.close()
    
    # Process output based on dp_cov -default[1]
    var_list= dict()
    for line in output.splitlines():
        if re.match('^#',line):
            continue
        line_list= (line.split('\t')) #0=header;1=pos;3=ref;4=var("," slipted);7=flags(;/= splited)
        
        if line_list[0] in var_list.keys():
            if int(re.search('DP=(.+?);',line_list[7]).group(1)) >= dp_cov:
                if ((line_list[4].split(",")[0].upper() == 'A') or
                (line_list[4].split(",")[0].upper() == 'C') or
                (line_list[4].split(",")[0].upper() == 'T') or 
                (line_list[4].split(",")[0].upper() == 'G')):
                    var_list[line_list[0]].append((int(line_list[1])-1,len(line_list[3]),len(line_list[4].split(",")[0]),line_list[4].split(",")[0]))
        else:
            if int(re.search('DP=(.+?);',line_list[7]).group(1)) >= dp_cov:
                if ((line_list[4].split(",")[0].upper() == 'A') or
                (line_list[4].split(",")[0].upper() == 'C') or
                (line_list[4].split(",")[0].upper() == 'T') or 
                (line_list[4].split(",")[0].upper() == 'G')):
                    var_list[line_list[0]] = [(int(line_list[1])-1,len(line_list[3]),len(line_list[4].split(",")[0]),line_list[4].split(",")[0])]

    new_fasta_temp = open(os.path.join(output_dir,'tmp','snp_fasta.fasta'),'w')    
    for seq_record in SeqIO.parse(os.path.join(output_dir,'tmp','v_'+str(count-1)+'.fasta'),'fasta'):
        pos_temp=0
        if seq_record.id in var_list.keys():
            seq_mutable = seq_record.seq.tomutable()
            for value in var_list[seq_record.id]:
                pos,len_ref,len_var,var = value
                seq_mutable[pos+pos_temp:pos+pos_temp+len_ref] = var
                pos_temp+=(len_var-len_ref)
            new_fasta_temp.write(">{}\n{}\n".format(seq_record.id,seq_mutable))
        else:
            new_fasta_temp.write(">{}\n{}\n".format(seq_record.id,seq_record.seq)) ##removi seq_mutable e substitui por seq_record.seq
    new_fasta_temp.close()     
    
    check_overlap(output_dir,os.path.join(output_dir,'tmp','snp_fasta.fasta'),av_readlen,count=count)


def remove_N(output_dir,name_fasta,fasta_semifinal,orig_target,fasta_len,av_readlen,mean_gap):
    fasta_wo_N = ""
    ext_size = av_readlen*3
    final_fasta = open(os.path.join(output_dir,'curated_seqs',name_fasta+'_curated.fasta'), 'w')
    print(fasta_semifinal)

    for seq_record in SeqIO.parse(fasta_semifinal,'fasta'):
        seq_mutable = seq_record.seq.tomutable()
        actual_start, actual_end = 0,0
        if seq_record.id not in orig_target.keys():        
            fasta_wo_N = remove_N_slave(seq_record.id,seq_mutable,actual_start,actual_end)
            final_fasta.write(str(fasta_wo_N))
        else:
            if orig_target[seq_record.id][0][0] == 1:
                actual_start = orig_target[seq_record.id][0][1] + ext_size
            if orig_target[seq_record.id][-1][1] == fasta_len[seq_record.id]:
                actual_end = orig_target[seq_record.id][-1][1] - orig_target[seq_record.id][-1][0] +1 + ext_size
            fasta_wo_N = remove_N_slave(seq_record.id,seq_mutable,actual_start,actual_end,av_readlen)
            final_fasta.write(str(fasta_wo_N))
    final_fasta.close()   

def remove_N_slave(fasta_header,seq_mutable,actual_start,actual_end,av_readlen):
    len_seq = len(seq_mutable)    
    new_fasta = ""
    ext_size = av_readlen*3
    for number in range (len_seq-(ext_size)- actual_end,len_seq):
        if seq_mutable[number] == "N":
            remov_end = number
            break
        else:
            remov_end= len_seq ## this considers if all edges were fullfilled by bases (same as below)
    for number in range ((ext_size)+actual_start-1,0,-1):
        if seq_mutable[number] == "N":
            remov_start = number+1 ### +1?
            break
        else:
            remov_start = 0 
    new_fasta += '>'+fasta_header+'\n'+seq_mutable[remov_start:remov_end]+'\n' 
    
    return new_fasta
