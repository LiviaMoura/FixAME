###########################################################################
#project: Genome curation - Mimicking Banfield's genome curation proccess 
#autor: Livia Maria Silva Moura
#e-mail: liviam.moura@gmail.com
# PDSE Scholarship UC Berkeley - PI's BR/USA: Setubal,J /Banfield, J
##########################################################################

import os,subprocess,re,sys,io
import argparse
import pysam as ps
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
import argparse

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
    indexes = [index for index,base in enumerate(one_fasta) if base == 'N']
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

        ###removing from the main N index, the pos in btw
        for item in temp_list:
             indexes.remove(item)

        ### merging regions if they are too close (80% of readlength)
        temp_list=list()
        for i in range(1,len(indexes)-1,2):
            if indexes[i+1] - indexes[i] < (av_readlen*0.8):
                temp_list.extend([i,i+1])
        
        for item in sorted(temp_list, reverse=True):
            del indexes[item]  

        ###Final N_pos - list with startN, end N, Number of N in total
        N_pos = list()
        N_target = list()
        
        for i in range(0, len(indexes), 2):
            N_pos.append(tuple((indexes[i],indexes[i+1],indexes[i+1]-indexes[i]+1))) ## return the python index mode
            N_target.append(list([indexes[i]+1,indexes[i+1]+1]))  ### return the position for calculate new target      
        return N_pos,N_target
    else:
        return list(),list()

def aligner(thread,ext,fasta="",bam_out="",read1="",read2="",minid=0.99,rescue_mm="",semi=False):
    subprocess.call(['bbmap.sh', 'ref='+str(fasta)+str(ext),'overwrite=t'], shell=False)
    if not semi:
        subprocess.call(['bbmap.sh', 'ref='+str(fasta)+str(ext), 'in1='+read1, 'in2='+read2, 
                  'outm='+bam_out+'.bam', 'threads='+str(thread),
                  'minid='+str(minid), 'rescuemismatches='+str(rescue_mm),'overwrite=t'], shell=False)
    else:
        subprocess.call(['bbmap.sh', 'ref='+str(fasta)+str(ext), 'in1='+read1, 'in2='+read2, 
                  'outm='+bam_out+'.bam', 'threads='+str(thread),'semiperfectmode=t','overwrite=t'], shell=False)
    subprocess.call(['samtools', 'sort', bam_out+'.bam', '-o', bam_out+'_sorted.bam'], shell=False)
    subprocess.call(['samtools', 'index', fasta+'_sorted.bam'], shell=False)
    os.remove(bam_out+'.bam')
    
def filtering_bam(thread,ext,num_mm,num_mm_filt,fasta="",read1="",read2=""):
    ### Filtering reads with less than 0 or [1] mismatch (default to be decided) 
    samfile = ps.AlignmentFile(fasta+'_sorted.bam',"rb")
    match_reads = list()
    for read in samfile.fetch():
        if not read.has_tag('NM'):
            continue
        if not (read.get_tag('NM') > num_mm): ## if TAG NM is not greater than num_mm (default = 1) -> read_list to be filtered
            match_reads.append(read.qname.split(" ")[0])
    samfile.close()

    ### creating a file with unique reads names that fulfill the criteria above 
    match_reads = list(set(match_reads)) 
    with open ('matched_reads', 'w') as outfile:
        outfile.write("\n".join(str(item) for item in match_reads))
    outfile.close()
    
    ### Filtering original fastq with the high stringency reads 
    subprocess.call(['filterbyname.sh', 'in='+read1, 'in2='+read2, 'out=res_R1.fastq', 'out2=res_R2.fastq', 'names=matched_reads', 'include=t'])
    
    aligner(thread,ext,fasta,read1="res_R1.fastq",read2="res_R2.fastq",bam_out=fasta+"_filtered",minid=0.6,rescue_mm=round(av_readlen/2))

    

def build_N(ext,num_mm_filt,fasta,bam_filt=""):
    subprocess.call(['samtools', 'index', bam_filt+'_sorted.bam'], shell=False)
    samfile = ps.AlignmentFile(bam_filt+'_sorted.bam',"rb")
    mm_bam = ps.AlignmentFile("MM.bam", "wb", template=samfile) ## creating the bam with mismatch for dp purpose
    
    ### Filtering the 2nd alignment (60% minid) to get reads with MM higher/equal than default=4
    for read in samfile.fetch():
        if not read.has_tag('NM'):
            continue
        if not (read.get_tag('NM') < num_mm_filt): ## if TAG NM not lesser the num_mm_filt (default = 4 mm) -> save a bam with mm reads
            mm_bam.write(read)
    mm_bam.close()
    samfile.close()

    ### Depth coverage comparation
    
    depth_cov_mm = subprocess.Popen("samtools depth MM.bam", shell=True, stdout=subprocess.PIPE,universal_newlines=True).communicate()[0]
    depth_cov_opt = subprocess.Popen("samtools depth {}_sorted.bam".format(bam_filt), shell=True, stdout=subprocess.PIPE,universal_newlines=True).communicate()[0]
    
    df_mm = pd.read_table(io.StringIO(depth_cov_mm), sep='\t', names=['contig','POS','qnt'])
    df_mm['contig'] = df_mm['contig'].str.split(" ", n = 1, expand = True)[0]
    
    #df_mm = df_mm.loc[df_mm['qnt'] >= num_mm_filt] - remove this line?
    
    df_opt = pd.read_table(io.StringIO(depth_cov_opt), sep='\t', names=['contig','POS','qnt'])
    df_opt['contig'] = df_opt['contig'].str.split(" ", n = 1, expand = True)[0]
    df_final = pd.merge(df_mm, df_opt, on=['contig','POS'], how='inner', suffixes=('_mm','_opt'))
    
    ### Dif btw 2nd alignment and MM/bam should be > 50%
    df_final['dif'] = df_final['qnt_mm']/df_final['qnt_opt']
    df_final = df_final.loc[df_final['dif'] > 0.50]
    dict_replace = df_final.groupby('contig')['POS'] .apply(list).to_dict()
    
    # Cleaning the dict_replace... it keep only first and last position and it'll creat a list of list to make it easier to work later
    dict_replace_chunked = defaultdict(list)
    
    if os.path.exists("target"): ### remove target file if it exists
        os.remove("target")
            
    
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
        
        ### merging regions if they are too close (80% of readlength)
        temp_list=list()
        for i in range(1,len(dict_replace[key])-1,2):
            if dict_replace[key][i+1] - dict_replace[key][i] < (av_readlen*0.8):
                temp_list.extend([i,i+1])        
        
        for item in sorted(temp_list, reverse=True):
            del dict_replace[key][item]           
        
        for i in range(0, len(dict_replace[key]), 2):
            dict_replace_chunked[key].append(tuple(dict_replace[key][i:i+2]))
        
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
            
            ### This is part of extension script - It adds av_readlen*3 "N"s in the beginning and end from a file
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
    
    return dict_error_pos


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


def check_overlap(new_fasta):
    fasta_temp = io.StringIO(str(new_fasta))
    new_fasta =""
    
    if os.path.exists("target"): ### remove target file if it exists
        os.remove("target")
    
    for i, seq_record in enumerate(SeqIO.parse(fasta_temp,'fasta')):
        new_target=list()
        N_pos,N_target_temp = check_npos(seq_record.seq)
        print(N_target_temp,"N_target_temp")
        
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
                alignments = pairwise2.align.localms(comp_left,comp_right,2,-1,-.5,-.1)
                
                if alignments[0][2] == 24:
                    print(start,end,number,"start,  end,   NumberofN, #######DEU MATCH")
                    variat=(av_readlen*2)-alignments[0][3]
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
               
        #### Creating new target based on control and N_target_temp ######
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

path = "/data7/ThiocyanateBioreactor/2018/curator/Livia/genomas/sample_Paula/curated"
os.chdir(path)
fasta_in = "ERMZT840_2_Bacteria_67_13_curated.contigs.fa"
fasta,ext = os.path.splitext(fasta_in)
read1_in = "10968.7.186300.TAGTGAC-GGTCACT_trim_clean.PE.1.fastq.gz"
read2_in = "10968.7.186300.TAGTGAC-GGTCACT_trim_clean.PE.2.fastq.gz"
thread = 6
num_mm = 1
num_mm_filt = 4
minid_in = 0.99
x_times = 8
#alele_freq = 0.6
dp_cov=1

av_readlen = 150
#aligner(thread,ext,fasta,read1=read1_in,read2=read2_in,bam_out=fasta,minid=minid_in,rescue_mm=round(av_readlen/3))
#filtering_bam(thread,ext,num_mm,num_mm_filt,fasta,read1=read1_in,read2=read2_in)
target = build_N(ext,num_mm_filt,fasta,bam_filt=fasta+"_filtered")
print(target,"FIRST")
print(" "),
number_times = count_times(x_times)

for count in number_times :
    print(count)
    var_cal_fix(count,ext,thread,x_times,dp_cov,target)


    

