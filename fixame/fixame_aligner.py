import argparse
import os
import logging
import sys
import subprocess
import fixame
from fixame.fixame_common import script_path 

__author__ = "Livia Moura"
__copyright__ = "Copyright 2019"
__maintainer__ = "Livia Moura"
__email__ = "liviam.moura@gmail.com"
__status__ = "Development"

def aligner(output_dir,thread,minid,fasta_path,r1,r2,r12,bam_out,semi=False,last=False):
    varpath = script_path()
    subprocess.call([os.path.join('bbmap.sh'), 'ref='+str(fasta_path),'path='+output_dir+'/new_fastas','overwrite=t'], shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,)
    if not semi:
        if not r12:
            subprocess.call([os.path.join('bbmap.sh'), 'ref='+str(fasta_path),'path='+output_dir+'/new_fastas', 'in1='+r1, 'in2='+r2, 
                  'outm='+output_dir+'/tmp/'+bam_out+'.bam', 'threads='+str(thread),
                  'minid='+str(minid),'ambiguous=random','overwrite=t'], shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, )
        else:
            subprocess.call([os.path.join('bbmap.sh'), 'ref='+str(fasta_path),'path='+output_dir+'/new_fastas', 'in='+r12,
                  'outm='+output_dir+'/tmp/'+bam_out+'.bam', 'threads='+str(thread),
                  'minid='+str(minid),'ambiguous=random','overwrite=t'], shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,)
    else:
        subprocess.call([os.path.join('bbmap.sh'), 'ref='+str(fasta_path),'path='+output_dir+'/new_fastas', 'in1='+r1, 'in2='+r2, 
                  'outm='+output_dir+'/tmp/'+bam_out+'.bam', 'threads='+str(thread),'ambiguous=random','killbadpairs=t','semiperfectmode=t','overwrite=t'], shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,)
        
    subprocess.call(['samtools', 'sort', output_dir+'/tmp/'+bam_out+'.bam', '-o', output_dir+'/tmp/'+bam_out+'_sorted.bam','--threads',str(thread)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,)
    subprocess.call(['samtools', 'index', output_dir+'/tmp/'+bam_out+'_sorted.bam'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,)
    os.remove(output_dir+'/tmp/'+bam_out+'.bam')
