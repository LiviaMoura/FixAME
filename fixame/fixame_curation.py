import argparse
import os
import logging
import sys
import re
import subprocess
import fixame
import shutil
import pysam as ps
from datetime import datetime
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqIO.FastaIO import SimpleFastaParser

from fixame.fixame_aligner import aligner
from fixame.fixame_common import *
from fixame.fixame_error_finder import *

from xopen import xopen

__author__ = "Livia Moura"
__copyright__ = "Copyright 2019"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"

def fixame_curation_validate(**kwargs):
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
    '''

    method = common_validate(**kwargs)
   
    if (kwargs.get('minid') < 0.76) or (kwargs.get('minid') > 1.00):
        logging.info('Checking minimum identify for the first alignment')
        logging.error("Please, verify if -minid is >= 0.76 and <= 1.00 ")
        sys.exit()

    if (kwargs.get('num_mismatch') < 0) or (kwargs.get('num_mismatch') > 5):
        logging.info('Checking number of mismatches allowed for the filtering reads')
        logging.error("-num_mismatch must be 0>=x>=5")
        sys.exit()

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
    
    av_readlen = temp_average_read(kwargs.get('r1'))
    output_dir = kwargs.get('output_dir')

    try:
        mydir = os.path.join(output_dir,'fixame_'+datetime.now().strftime('%Y-%b-%d_%H-%M-%S'))
        os.mkdir(mydir)
        os.mkdir(os.path.join(mydir,'tmp'))
        os.mkdir(os.path.join(mydir,'new_fastas'))
        logging_config(mydir)
        logging.info("Fixame output folder created - {}".format(mydir))
    except:
        logging.error("It wasn't possible to create fixame output folder")

    # Checking the pipeline - genome/metagenoms vs Bins
    if method == 0:
        fasta_in = os.path.realpath(os.path.expanduser(kwargs.get('fasta')))
        name_fasta = os.path.splitext(os.path.basename(fasta_in))[0]

        logging.info("\n --- Analysing the file {} ---\n".format(name_fasta))
        
                
        try:
            logging.info("Checking overlaping at N regions on {} and fix them".format(fasta_in))
            check_overlap(mydir,fasta_in,av_readlen,True)
            logging.info("A new reference fasta {} was created".format(mydir+'/new_fastas/'+name_fasta+'_renewed.fasta'))
        except:
            logging.error("Something went wrong")
            sys.exit()
        
        
        try:
            logging.info("Mapping reads against the new reference")
            aligner(mydir,kwargs.get('threads'),kwargs.get('minid'),mydir+'/new_fastas/'+name_fasta+'_renewed.fasta',r1=read1_in,r2=read2_in,r12=read12_in, bam_out=name_fasta+'_renewed')                
        except:
            logging.error("Something went wrong")
            sys.exit()
        
        fasta_cov, num_mm = kwargs.get('fasta_cov'), kwargs.get('num_mismatch')

        #Filtering the fastq - Make curation process faster
        try:
            logging.info("Filtering out reads that didn't map the new ref")
            filtering_bam(mydir,kwargs.get('threads'),kwargs.get('num_mismatch'),mydir+'/tmp/'+name_fasta+'_renewed',read1_in,read2_in,read12_in)                
        except:
            logging.error("Something went wrong")
            sys.exit()
        
        try:
            logging.info("Generating some metrics to keep running")
            reference_to_length = calculate_reference_lengths(mydir+'/new_fastas/'+name_fasta+'_renewed.fasta', minimum_assembly_length)

            (   bam_dict, 
                reference_read_lengths,
                average_template_length,
                average_read_length,
                average_gap_length,
                template_length_min,
                template_length_max,
                average_gap_std
            ) = parse_map(mydir+'/tmp/'+name_fasta+'_renewed_sorted.bam', num_mm, kwargs.get('threads'), minimum_assembly_length, reference_to_length)
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
        except: 
            logging.error("Something went wrong")  
            sys.exit()        
        
        
        try:
            logging.info("Selecting the errors regions")
            orig_target,fasta_len = build_N(mydir,kwargs.get('threads'),mydir+'/new_fastas/'+name_fasta+'_renewed.fasta',average_read_length,reference_to_error_regions)
            
        except:
            logging.error("Something went wrong")
            sys.exit()

            
        logging.info("Starting to fix sample {}\n".format(name_fasta))      
        os.mkdir(os.path.join(mydir,'fixing_log'))
        logging.info("Folder {} was created".format(os.path.join(mydir,'fixing_log')))      
            
        for count,r in enumerate(range(kwargs.get('xtimes')),1):
            fixed = open(os.path.join(mydir,'fixing_log','fixame_loop_'+str(count)+'.txt'),'w+')
            try:
                logging.info("Loop {} from {}".format(count,kwargs.get('xtimes')))
                var_cal_fix(mydir,count,fixed,kwargs.get('threads'),kwargs.get('xtimes'),kwargs.get('dp_cov'),orig_target,fasta_len,av_readlen)
            except:
                logging.error("Something went wrong")
                sys.exit()
            fixed.close()
    
        os.mkdir(os.path.join(mydir,'fixame_results'))
        #print(orig_target,'             ',average_gap_length, average_read_length, average_gap_std)
        remove_N(mydir,name_fasta,os.path.join(mydir,'tmp','v_'+str(kwargs.get('xtimes'))+'.fasta'),orig_target,fasta_len,av_readlen, average_gap_length, average_gap_std, kwargs.get('threads'))

        if (kwargs.get('keep') == False):
            shutil.rmtree(os.path.join(mydir,'tmp'))
        

    else: # BINS MODE
        fasta_array = []
        name_sample = 'bins'
        contigs_bins = open(os.path.join(mydir,'bin_contigs.txt'), 'w+')
        merged = open(os.path.join(mydir,'tmp','bins.fasta'),'w+')
        for sample in os.listdir(os.path.realpath(os.path.expanduser( kwargs.get('bins') ))):
            name = sample.split('.')[0]
            if sample.lower().endswith(".fasta") or sample.lower().endswith(".fa") or sample.lower().endswith(".fna"):
                fasta_array.append(sample)
                for seq_record in SeqIO.parse(os.path.join(kwargs.get('bins'),sample),'fasta'):
                    contigs_bins.write("{}\t{}\n".format(name, seq_record.id))
                with open(os.path.join(kwargs.get('bins'),sample), 'r') as readfile:
                    shutil.copyfileobj(readfile, merged)
        contigs_bins.close()
        merged.close()

        fasta_in = os.path.join(mydir,'tmp','bins.fasta')

        logging.info("\n --- Analysing a metagenome sample with {} bins ---\n".format(len(fasta_array)))
        try:
            logging.info("Checking overlaping at N regions on {} and fix them".format(fasta_in))
            check_overlap(mydir,fasta_in, av_readlen,True)
            logging.info("A new reference fasta {} was created".format(mydir+'/new_fastas/'+name_sample+'_renewed.fasta'))
        except:
            logging.error("Something went wrong")
            sys.exit()

        #print(mydir,kwargs.get('threads'),kwargs.get('minid'),mydir+'/new_fastas/'+name_sample+'_renewed.fasta',read1_in,read2_in,read12_in, name_sample+'_renewed')

        try:
            logging.info("Mapping reads against the new reference")
            
            aligner(mydir,kwargs.get('threads'),kwargs.get('minid'),mydir+'/new_fastas/'+name_sample+'_renewed.fasta',r1=read1_in,r2=read2_in,r12=read12_in, bam_out=name_sample+'_renewed')                
        except:
            logging.error("Something went wrong")
            sys.exit()
        
        fasta_cov, num_mm = kwargs.get('fasta_cov'), kwargs.get('num_mismatch')

        #Filtering the fastq - Make curation process faster
        try:
            logging.info("Filtering out reads that didn't map the new ref")
            filtering_bam(mydir,kwargs.get('threads'),kwargs.get('num_mismatch'),mydir+'/tmp/'+name_sample+'_renewed',read1_in,read2_in,read12_in)                
        except:
            logging.error("Something went wrong")
            sys.exit()
        
        try:
            logging.info("Generating some metrics to keep running")
            reference_to_length = calculate_reference_lengths(mydir+'/new_fastas/'+name_sample+'_renewed.fasta', minimum_assembly_length)
            print('opa')
            (   bam_dict, 
                reference_read_lengths,
                average_template_length,
                average_read_length,
                average_gap_length,
                template_length_min,
                template_length_max,
                average_gap_std,
            ) = parse_map(mydir+'/tmp/'+name_sample+'_renewed_sorted.bam', num_mm, kwargs.get('threads'), minimum_assembly_length, reference_to_length)
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
        except: 
            logging.error("Something went wrong") 
            sys.exit()         
        
        
        try:
            logging.info("Selecting the errors regions")
            orig_target,fasta_len = build_N(mydir,kwargs.get('threads'),mydir+'/new_fastas/'+name_sample+'_renewed.fasta',average_read_length,reference_to_error_regions)
            
        except:
            logging.error("Something went wrong")
            sys.exit()

       
        logging.info("Starting to fix all bins\n")      
        os.mkdir(os.path.join(mydir,'fixing_log'))
        logging.info("Folder {} was created".format(os.path.join(mydir,'fixing_log')))      
            
        for count,r in enumerate(range(kwargs.get('xtimes')),1):
            fixed = open(os.path.join(mydir,'fixing_log','fixame_loop_'+str(count)+'.txt'),'w+')
            try:
                logging.info("Loop {} from {}".format(count,kwargs.get('xtimes')))
                var_cal_fix(mydir,count,fixed,kwargs.get('threads'),kwargs.get('xtimes'),kwargs.get('dp_cov'),orig_target,fasta_len,av_readlen)
            except:
                logging.error("Something went wrong")
                sys.exit()
            fixed.close()
    
        os.mkdir(os.path.join(mydir,'fixame_results'))
        #print(orig_target,'             ',average_gap_length)
        remove_N(mydir,name_sample,os.path.join(mydir,'tmp','v_'+str(kwargs.get('xtimes'))+'.fasta'),orig_target,fasta_len,av_readlen, average_gap_length, average_gap_std, kwargs.get('threads'))

        if (kwargs.get('keep') == False):
            try:
                logging.info("Removing temporary files")
                shutil.rmtree(os.path.join(mydir,'tmp'))
            except:
                logging.info("It wasn't possible to remove the /tmp folder")                

        ## Spliting the bins
        df = pd.read_table(os.path.join(mydir,'bin_contigs.txt'), header=None)
        df = df.groupby(0).agg({1:lambda x: list(x)}).reset_index()
        bin_ctg_dict = dict(zip(df[0], df[1]))

        for k,v in bin_ctg_dict.items():
            per_bin = open(os.path.join(mydir,'fixame_results', k+'_fixame.fasta'),'w+')
            for contig in v:
                for record in SeqIO.parse(os.path.join(mydir,'fixame_results','bins_fixame.fasta'),'fasta'):
                    if record.id == contig:
                        per_bin.write(">{}\n{}\n".format(contig, record.seq))
            per_bin.close()
            
            per_bin = open(os.path.join(mydir,'new_fastas', k+'_renewed.fasta'),'w+')
            for contig in v:
                for record in SeqIO.parse(os.path.join(mydir,'new_fastas','bins_renewed.fasta'),'fasta'):
                    if record.id == contig:
                        per_bin.write(">{}\n{}\n".format(contig, record.seq))
            per_bin.close()
            
        os.remove(os.path.join(mydir,'fixame_results', k+'_fixame.fasta'))
        os.remove(os.path.join(mydir,'new_fastas','bins_renewed.fasta'))


         
           
#    os.path.realpath(os.path.expanduser)

def logging_config(workdir):
    '''Configure the logging'''
    log_file = os.path.join(workdir, "fixame.log")
    log_formatter = logging.Formatter("[%(levelname)s] %(asctime)s %(message)s")
    root_logger = logging.getLogger()
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)
    root_logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(log_formatter)
    root_logger.addHandler(ch)
    root_logger.propagate = False


def temp_average_read(r1_fastq):
    '''Initial fast avg read length - It'll be recalculated with more precision later'''
    sum_read_len=0
    for num,record in enumerate(SeqIO.parse(xopen(r1_fastq), 'fastq')):
        sum_read_len += len(record.seq)
        if num == 9999:
            av_read_len = int(sum_read_len/10000)
            break
    return av_read_len

def check_overlap(output_dir,fasta,av_readlen,user_file=False,fixed='',count=''):
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
            
            if j==0: ### first time doesnt have index modification
                temp_left = seq_record[0:start].seq  ### start = N pos; but here would be Last base pos
                temp_right = seq_record[end+1:-1].seq ### end = N pos; but here would be first base pos
                comp_left = temp_left[-(av_readlen*2):]
                comp_right = temp_right[0:12]
                
                ## Self remember alignment list 0 to N
                if comp_right in comp_left:
                    alignments = pairwise2.align.localms(comp_left,comp_right,2,-1,-.5,-.1)
                    if fixed != '':
                        fixed.write("{}\t{}\t{}\tfixed\n".format(seq_record.id,start,end))                   
                    #print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######DEU MATCH")                 #############################   
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
                        if fixed != '':
                            fixed.write("{}\t{}\t{}\tfixed\n".format(seq_record.id,start,end))                   
                        #print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######ESPECIAL1")                    ###########################    
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
                    if fixed != '':
                        fixed.write("{}\t{}\t{}\tfixed\n".format(seq_record.id,start,end)) 
                    #print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######DEU MATCH")                     #####################
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
                        if fixed != '':
                            fixed.write("{}\t{}\t{}\tfixed\n".format(seq_record.id,start,end)) 
                        #print(seq_record.id,start,end,number,"start,  end,   NumberofN, #######ESPECIAL2")              ###################
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
    if fixed != '':
        return fixed       

   
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
    #varpath=script_path()
    subprocess.call([os.path.join('filterbyname.sh'), 'in='+r1, 'in2='+r2, 'out='+output_dir+'/tmp/res_R1.fastq', 'out2='+output_dir+'/tmp/res_R2.fastq', 'names='+output_dir+'/tmp/matched_reads', 'include=t','overwrite=t'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,)

def build_N(output_dir,threads,fasta,av_readlen,dict_replace_0):
    ''' Replace error for N and extend the contig's edges for N - pre-curation step'''
    
    if os.path.exists(os.path.join(output_dir,'tmp',"target")): ### remove target file if it exists
            os.remove(os.path.join(output_dir,'tmp',"target"))    
    
    dict_replace = {key: value for key, value in dict_replace_0.items() if value}
    
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

    #print(dict_only_errors,'dict_only_errors')
    with open(os.path.join(output_dir,'Fixame_Errors_report.txt'),'w+') as error_loc:
        error_loc.write("contig_name\terror_start\terror_end\tn_affected_bases\ttype_of_error\n")
        counter_err, counter_contigs = 0,0
        for key,value in dict_only_errors.items():
            logging.debug("{} possibly local errors at contig {}".format(len(value), key))
            counter_err += len(value)
            counter_contigs += 1
            for item in value:
                error_loc.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(key,item[0],item[1],abs(item[1]-item[0]),'Local_error'))
        error_loc.close()
    logging.warning("\nFixame could detect a total of {} errors in {} contig(s)\n".format(counter_err,counter_contigs))
    logging.info("The file {} was created".format(output_dir+'/Fixame_Errors_report.txt'))                    

    return dict_replace_chunked,dict_len


def var_cal_fix(output_dir,count,fixed,thread,x_times,dp_cov,orig_target,fasta_len,av_readlen):    
    '''Find changes on Ns position and replace them'''
    
    aligner(output_dir,thread,0,os.path.join(output_dir,'tmp','v_'+str(count-1)+'.fasta'),r1=os.path.join(output_dir,'tmp','res_R1.fastq'),r2=os.path.join(output_dir,'tmp','res_R2.fastq'),r12="", bam_out='v_'+str(count-1),semi=True)
    varpath = script_path()
    
    # Run variant finder - bcftools
    
    comm = "bcftools mpileup -R {}/target -A -B -f {}/v_{}.fasta {}/v_{}_sorted.bam".format(os.path.join(output_dir,'tmp'),os.path.join(output_dir,'tmp'),count-1,os.path.join(output_dir,'tmp'),count-1)
    varfind = subprocess.Popen(comm, shell=True, stdout=subprocess.PIPE,universal_newlines=True, stderr=subprocess.DEVNULL,)
    output = varfind.communicate()[0]
    #print(f'AQUI TEMOS O OUTPUT {output}')
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
    fixed = check_overlap(output_dir,os.path.join(output_dir,'tmp','snp_fasta.fasta'),av_readlen,fixed=fixed,count=count)
    return fixed


def remove_N(output_dir,name_fasta,fasta_semifinal,orig_target,fasta_len,av_readlen,mean_gap, mean_gap_std, thread):
    fasta_wo_N = ""
    ext_size = av_readlen*3
    final_fasta = open(os.path.join(output_dir,'fixame_results',name_fasta+'_fixame.fasta'), 'w')
    
    # Alignment needed to selected reads #fasta_semifinal
    aligner(output_dir, thread, 1,os.path.join(output_dir,'tmp','v_0.fasta'),r1=os.path.join(output_dir,'tmp','res_R1.fastq'),r2=os.path.join(output_dir,'tmp','res_R2.fastq'),r12="", bam_out='check_read',semi=True)
    
    for seq_record in SeqIO.parse(fasta_semifinal,'fasta'):
        seq_mutable = seq_record.seq.tomutable()
        
        actual_start, actual_end = 0,0
        if seq_record.id not in orig_target.keys():        
            fasta_wo_N = remove_N_slave(seq_record.id,seq_mutable,actual_start,actual_end,av_readlen)
            final_fasta.write(str(fasta_wo_N))
        else:
            N_pos,N_target_temp = check_npos(seq_record.seq,av_readlen)
            N_pos = N_pos[1:-1]
            #print(N_pos, 'NPOSSSS')

            #Middle
            if N_pos:
                seq_mutable = check_reads_N_edges(output_dir, seq_record.id, seq_mutable, seq_record.id, av_readlen,mean_gap, mean_gap_std, N_pos)

            #Edges
            if orig_target[seq_record.id][0][0] == 1:
                actual_start = orig_target[seq_record.id][0][1] + ext_size
            if orig_target[seq_record.id][-1][1] == fasta_len[seq_record.id]:
                actual_end = orig_target[seq_record.id][-1][1] - orig_target[seq_record.id][-1][0] +1 + ext_size
            fasta_wo_N = remove_N_slave(seq_record.id,seq_mutable,actual_start,actual_end,av_readlen)
            final_fasta.write(str(fasta_wo_N))
    final_fasta.close()   

def check_reads_N_edges(output_dir, contig_name, seq_mutable, seq_name, av_readlen,mean_gap, mean_gap_std, N_pos):
    r_left = os.path.join(output_dir,'tmp', "r_left")
    r_right = os.path.join(output_dir,'tmp', "r_right")
    left_right = os.path.join(output_dir,'tmp', "left_right")
    count = 0 
    
    mean_frag_len = 2*av_readlen + mean_gap
    no_support = open(os.path.join(output_dir,'fixing_log','fixame_without_read_support.txt'),'a')
    for start, end, space in N_pos:
        #print(start,'START', end,'END', space,'SPACE')
        cmd = '''samtools view {}/check_read_sorted.bam {}:{}-{} | \
                grep -v "*" |cut -f 1 | sort | uniq -u > {}'''.format(os.path.join(output_dir,'tmp'), seq_name,(start-mean_gap),start,r_left)
        subprocess.run(cmd, shell=True,)
        cmd = '''samtools view {}/check_read_sorted.bam {}:{}-{} | \
                grep -v "*" |cut -f 1 | sort | uniq -u > {}'''.format(os.path.join(output_dir,'tmp'), seq_name,end,(end+2*mean_gap), r_right)
        subprocess.run(cmd, shell=True,)
        cmd = '''cat {} {}| sort | uniq -d > {}'''.format(r_left, r_right, left_right)
        subprocess.run(cmd, shell=True,)

        if os.stat(left_right).st_size == 0:
            no_support.write("No read support around:\t{}:{}-{}\n".format(contig_name,start,end,))
            continue
        else:
            cmd = '''samtools view {}/check_read_sorted.bam {}:{}-{} | \
                    grep -F -f {}| awk -F $'\t' '$9 > 0 {{ sum += $9; n++ }} END {{print int(sum/n)}}' '''.format(os.path.join(output_dir,'tmp'), seq_name,(start-2*av_readlen), (end+2*av_readlen), left_right)
            distance = int(subprocess.check_output(cmd,universal_newlines=True, shell = True).split()[0])
            
            if distance < (2*av_readlen + (mean_gap + mean_gap_std) ): #rare but possible
                seq_mutable[start+count:start+count] = (mean_frag_len - distance)*"N"                  
                count+= mean_frag_len - distance
            elif distance > (2*av_readlen + (mean_gap - mean_gap_std) ):
                seq_mutable[start+count:start+count+(distance - mean_gap_std)] = ''                
                count-= (distance - mean_gap_std)
            else:
                continue
    no_support.close()
    return seq_mutable

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