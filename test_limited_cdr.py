#!/opt/annaconda/python3 python

import re
import numpy as np
import os,sys
import json
from subprocess import Popen, PIPE

#Usegae: python test_limited_cdr.py test /usr/local/rsrch1/epi/jchen15/filter_variants/real_data/ov_vaast/ov_combined_cdr/tail_noSNV.cdr
#f=sys.argv[1]
#cdr_file=sys.argv[2]
#cdr_bed=sys.argv[3]
def get_limited_bed(f, cdr_file):
    tmp=open("temp_001.bed","w")
    for j in cdr_file:
        j=j.strip()
        j_list=j.split()
        #ref=j_list[5]
        if re.match("#",j_list[0]):
            pass
        elif re.match("chr",j) and (j_list[3]=="SNV" or j_list[3]=="complex_substitution"):
           #tmp.write("\t".join(j_list[0:3])+"\n")
	    tmp.write(j_list[0][3:]+"\t"+"\t".join(j_list[1:3])+"\n")
        elif re.match("chr",j) and (j_list[3]=="deletion"):# I will have to modeify this line because deletion the start should minus the number of delected basepairs.
            j_list[1]=str(int(j_list[1])-1)
            #tmp.write("\t".join(j_list[0:3])+"\n")
            tmp.write(j_list[0][3:]+"\t"+"\t".join(j_list[1:3])+"\n")
        elif re.match("chr",j) and (j_list[3]=="insertion"): #The insertion in cdr file is the same as vcf file.
            #tmp.write("\t".join(j_list[0:3])+"\n")
            tmp.write(j_list[0][3:]+"\t"+"\t".join(j_list[1:3])+"\n")
        else:pass
    tmp.close()
    log=open("test_output2.af","w")
    log.flush()
    #p=Popen(['bedtools intersect -a /usr/local/rsrch1/epi/jchen15/database/gnomAD/script/test_chr2.af -b temp_001.bed'], shell=True, stderr=PIPE, stdout=log )
##    p=Popen(['bedtools intersect ','-a','/usr/local/rsrch1/epi/jchen15/database/gnomAD/script/test_chr2.af' ,'-b', 'temp_001.bed'], shell=True, stderr=PIPE, stdout=log )#remove_this
    
#    p=Popen(['bedtools intersect -a test_chr2.af -b temp_001.bed >test_output2.af'], shell=False, stderr=PIPE, stdout=PIPE)

#    bedtools intersect -a test_chr2.af -b test.bed >test_output2.af
    #p = Popen(['tail','-1',f],shell=False, stderr=PIPE, stdout=PIPE)
    #number,err = p.communicate()
    #if err:
    #    print (err.decode())
    #else:
    #    # Use split to get the sample number in the cdr
    #    number = int(number.decode().split()[2].strip())
    #return(number)

if __name__ == "__main__":
    f=sys.argv[0]
    cdr_file=open(sys.argv[2])
    get_limited_bed(f,cdr_file)
    #f.close()
    #cdr_file.close()
