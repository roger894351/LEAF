#!/opt/anacnoda/python

## This script read in the Allele frequency and Allele counts from the external data and generate a json dictionary.


import sys
import re
import time
import json


start = time.time()


def check_sequence(x):# This is a function that returns how many str are same and different
    same_str=0
    diff_str=0
    for i in range(len(line_list[3])):
        if line_list[3][i]==line_list[4][i]:
            same_str+=1
        else:diff_str+=1    
    return(same_str,diff_str)


with open(sys.argv[1],"r") as f:
    a={}# this is a dictionary for variant_id
    for i in range(1,23):
        a[str(i)]={}
    for i in "X","Y","M":
        a[str(i)]={}
    for line in f:
        if re.match("#",line):
            pass
        else:
            line.strip()
            line_list=line.split()
            if len(line_list[3])==len(line_list[4]):# this is the snv or even number or variant changes
                if len(line_list[3])==1:
                    pass
                else:
                    if line_list[3][0]==line_list[4][0]:# the first genotype is same
                        same_num, diff_num=check_sequence(line_list)
                        line_list[1]=str(int(line_list[1])+same_num)
                        line_list[2]=str(int(line_list[2])+same_num)
                        line_list[3]=line_list[3][same_num:]
                        line_list[4]=line_list[4][same_num:]
                    elif line_list[3][-1]==line_list[4][-1]:#the last genotype is same
                        same_num, diff_num=check_sequence(line_list)
                        line_list[2]=str(int(line_list[2])-same_num)
                        line_list[3]=line_list[3][0:diff_num]
                        line_list[4]=line_list[4][0:diff_num]
                    else: print(line_list[3],"program can not handle this SNV")
            elif len(line_list[3])>len(line_list[4]):# This is the deletion: This has been take care when I convert to the cdr location
                pass

            elif len(line_list[3])<len(line_list[4]):# This is the insertion
                pass
            key, value=str("_".join(line_list[0:2])+"_"+"_".join(line_list[3:5])),",".join(line_list[6:9])
            a[str(line_list[0])][key]=value


with open('data_counts_split_chr.json', 'w') as fp:
    json.dump(a,fp)


f.close()
fp.close()
