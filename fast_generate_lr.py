#!/usr/local/python
# pwd:/usr/local/rsrch1/epi/jchen15/database/gnomAD/script
#usage:python call_bedtools.py test_head.cdr

import subprocess
import os,sys,re
import numpy as np
import time
import json
from subprocess import Popen, PIPE

def gene_vaast_score(input_allele, log_casm= None, score_only_minor = True,
                     score_only_alternative = True, site_penalty=2.0):

    case_allele1,case_allele0,control_allele1,control_allele0=input_allele[0:4]
    #case_allele1=input_allele[0],case_allele0=input_allele[1]
    ##roger added part
    #print("case_allele1"+str(case_allele1))
    #print("control_allele1"+str(control_allele1))
    # calculate site lrt score (array)
    log_lh = lrt_score(case_allele1, control_allele1,
            case_allele0, control_allele0)
    return(log_lh)
def lrt_score(case_allele1, control_allele1, case_allele0, control_allele0):
    ''' Calculate the unpenalized log likelihood of a gene
        INPUT: allele counts in cases and controls
        OUTPUT: site lrt scores
    '''
    alt_control_freq = 1.0 * control_allele1 / (control_allele0 + control_allele1)
    alt_case_freq = 1.0 * case_allele1 / (case_allele0 + case_allele1)

    null_freq = 1.0 * (case_allele1 + control_allele1) / (control_allele0 +
                control_allele1 + case_allele0 + case_allele1)

    alt_log_lh = log_likelihood(alt_control_freq, control_allele0, control_allele1
                    ) + log_likelihood(alt_case_freq, case_allele0, case_allele1)
    null_log_lh= log_likelihood(null_freq, control_allele0, control_allele1
                    ) + log_likelihood(null_freq, case_allele0, case_allele1)
    #print("the lrt score="+str(alt_log_lh - null_log_lh))
    return (alt_log_lh - null_log_lh)

def log_likelihood(freq, allele0, allele1):
    ''' Given specified frequencies, allele0 count, allele1 count, calculate
        log likelihood
        INPUT: specified frequencies, allele0 count, allele1 count
        OUTPUT: log likelihood from binomial distribution
    '''

    # bound frequency estimate to [1e-9< 1.0-1e-9] to avoid numerical issue
    freq = np.clip(freq, 1e-9, 1.0-1e-9)

    # calculate log likelihood
    loglik = allele1 * np.log(freq) + allele0 * np.log(1.0-freq)

    return loglik


def continous_ids(cdr_ids):
        cdr_list=[]
        for i in cdr_ids:
                if "-" in i:
                        a=i.split("-")
                        b=range(int(a[0]),int(a[1])+1)
                        for num in b:
                                cdr_list.append(str(num))
                else:
                        cdr_list.append(str(i))
        return(cdr_list)
####json1_data = json.loads(json1_str) ##create a json_dictionary
#json2=sys.argv[3]
#casm_file=open("/usr/local/rsrch1/epi/jchen15/apps/VAAST-CASM/casm_chr.json")
#casm_file=open(json2)
#casm_str=casm_file.read()
#casm_data=json.loads(casm_str)

#f="insertion.cdr"
#geno = np.array([[1,0]])
f=sys.argv[1]
def get_cdr_n(f): #give a cdr and return the last sample number is the cdr.
    #f = 'combined_case_control.cdr'
    # Get the last line from the file
    p = Popen(['tail','-1',f],shell=False, stderr=PIPE, stdout=PIPE)
    number,err = p.communicate()
    if err:
        print (err.decode())
    else:
        # Use split to get the sample number in the cdr
        number = int(number.decode().split()[2].strip())
    return(number)

def exac_counts(x): #x =snp_id
    chr_n=x.split("_")[0]
    if x in json1_data[chr_n].keys():
        alt,ref=int(json1_data[chr_n][x].split(",")[0]),int(json1_data[chr_n][x].split(",")[1])
        return(alt,ref)
    else: return(0,0)
def generate_score_input(j_list):
    ref_gt=j_list[5].split("|")[0]# whatever will be true
    n_case=get_cdr_n(f)# delete no call individuals
    case_allele1=0# set 0 first; #case_allele0=2*int(n_case);# get these numbers out:control_allele1=0;control_allele0=0
    gene_ann_index=1#if len(j_list[-1].split("|"))==3:
    snp_id=j_list[0][3:]+"_"+j_list[1]+"_"+ref_gt+"_"+j_list[6].split("|")[gene_ann_index].split(":")[1]
    snp_id_s=j_list[0][3:]+"_"+j_list[1]
    for k in range(6,len(j_list)):
        #print(k)# split the genotype:
        #snp_id="0"# this is snp_id reset
        #snp_id_s="0"        #print(j_list[k].split("|")[gene_ann_index].split(":")[0],j_list[k].split("|")[gene_ann_index].split(":")[1])
        #print(j_list[0:6],j_list[k].split("|")[gene_ann_index].split(":"))
        if j_list[k].split("|")[gene_ann_index]=="^:^":
            miss_count=len(continous_ids(j_list[k].split("|")[0].split(",")))##I change here
            n_case-=miss_count
            if snp_id!="1":
                snp_id=j_list[0][3:]+"_"+"_".join(j_list[1:3])+"_"+ref_gt+"_"+"^"#+j_list[k].split("|")[gene_ann_index].split(":")[0]: if a variant has a genotype other than no call it should has a snp_id other than 1
                snp_id_s=j_list[0][3:]+"_"+j_list[1]
        elif j_list[k].split("|")[gene_ann_index].split(":")[0]==j_list[k].split("|")[gene_ann_index].split(":")[1] and j_list[k].split("|")[gene_ann_index].split(":")[0]!=ref_gt:# this is a homozygus hit
            #snp_id=j_list[0][3:]+"_"+j_list[1]+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[0]
            #snp_id_s=j_list[0][3:]+"_"+j_list[1]#homozygotes
            snp_id=j_list[0][3:]+"_"+"_".join(j_list[1:3])+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[1]
            snp_id_s=j_list[0][3:]+"_"+j_list[1]
            case_allele1+=2*int(len(continous_ids(j_list[k].split("|")[0].split(","))))
        elif j_list[k].split("|")[gene_ann_index].split(":")[0]!=j_list[k].split("|")[gene_ann_index].split(":")[1] and j_list[k].split("|")[gene_ann_index].split(":")[0]==ref_gt:# this is a heterozygotes.
            case_allele1+=len(continous_ids(j_list[k].split("|")[0].split(",")))
            #snp_id=j_list[0][3:]+"_"+j_list[1]+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[1]
            #snp_s=j_list[0][3:]+"_"+j_list[1]
            snp_id=j_list[0][3:]+"_"+"_".join(j_list[1:3])+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[1]
            snp_id_s=j_list[0][3:]+"_"+j_list[1]
        else:
            pass#print("#more than 2 genotypes found",j)
                #print(j_list[k].split("|")[gene_ann_index].split(":"))
    #print(snp_id,snp_id_s,n_case,case_allele1,case_allele0)
    return(snp_id,snp_id_s,n_case,case_allele1,case_allele0)
##example:
#p1 = Popen(["dmesg"], stdout=PIPE)
#p2 = Popen(["grep", "hda"], stdin=p1.stdout, stdout=PIPE)
#p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
#output = p2.communicate()[0]


cdr_file=sys.argv[1]
temp=open("temp_outrandom.bed","w")
p1 = subprocess.Popen(["grep","-v","#", cdr_file], stdout=subprocess.PIPE) # standard output pipe to p2
p2 = subprocess.Popen(["cut","-f","1-3"], stdin=p1.stdout, stdout=subprocess.PIPE)
output = p2.communicate()[0]
temp.write(output)
temp.close()
#/usr/local/rsrch1/epi/jchen15/database/gnomAD/script
#command_line bedtools intersect -a AF_exome_chr.bed -b test.cdr.bed
dictionary=open("temp_dictionary.dic","w")
d=subprocess.Popen(["bedtools","intersect", "-a" ,"/usr/local/rsrch1/epi/jchen15/database/gnomAD/script/AF_exome_chr.bed","-b", "temp_outrandom.bed"],stdout = subprocess.PIPE)
d_data,d_error=d.communicate()

dictionary.write(d_data)
dictionary.close()
# read dictionary in a_data
json1_data={}# this is a dictionary for variant_id
for i in range(1,23):
    json1_data[str(i)]={}
for i in "X","Y","M":
    json1_data[str(i)]={}
for line in d_data.split("\n"):
    if re.match("#",line):
        pass
    elif len(line.split())<=4:
        pass
    else:
        line.strip()
        line_list=line.split()
        if re.match("chr",line):
            key, value=line_list[0][3:]+"_"+"_".join(line_list[1:3])+"_"+"_".join(line_list[4:6]),",".join(line_list[7:10])
            json1_data[str(line_list[0][3:])][key]=value
        else:
            key, value=line_list[0]+"_"+"_".join(line_list[1:3])+"_"+"_".join(line_list[4:6]),",".join(line_list[7:10])
            json1_data[str(line_list[0])][key]=value
#print(a["1"].keys())
#print(json1_data)



#del_temp="rm temp_outrandom.bed temp_dictionary.dic"
#subprocess.Popen(del_temp.split())
#read cdr file:
f=cdr_file
with open(f,"r") as cdr:
    for j in cdr: # this f is handeling a cdr file noted a vcf file
        j=j.strip()
        j_list=j.split("\t") #changed from "\t" to split() some cdr_line is space delimilated
        if re.match("#",j_list[0]):
            pass
        elif re.match("chr",j) and (j_list[3]=="SNV" or j_list[3]=="complex_substitution"):

            ref_gt=j_list[5].split("|")[0]# whatever will be true
            n_case=get_cdr_n(f)# delete no call individuals
            case_allele1=0# set 0 first; #case_allele0=2*int(n_case);# get these numbers out:control_allele1=0;control_allele0=0
            gene_ann_index=1#if len(j_list[-1].split("|"))==3:
            #snp_id="1" # weight has to output every variants
            for k in range(6,len(j_list)):# split the genotype:
                if j_list[k].split("|")[gene_ann_index]=="^:^":
                    miss_count=len(continous_ids(j_list[k].split("|")[0].split(",")))
                    n_case-=miss_count
                    snp_id=j_list[0][3:]+"_"+"_".join(j_list[1:3])+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[1]
                    snp_id_s=j_list[0][3:]+"_"+j_list[1]
                elif j_list[k].split("|")[gene_ann_index].split(":")[0]==j_list[k].split("|")[gene_ann_index].split(":")[1] and j_list[k].split("|")[gene_ann_index].split(":")[0]!=ref_gt:# this is a homozygus hit
                    #alt_gt=j_list[k].split("|")[1]
        #            print("the alt gt=",alt_gt)
                    snp_id=j_list[0][3:]+"_"+"_".join(j_list[1:3])+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[1]
                    snp_id_s=j_list[0][3:]+"_"+j_list[1]#homozygotes
                    #print(snp_id)
                    case_allele1+=2*int(len(continous_ids(j_list[k].split("|")[0].split(","))))
                elif j_list[k].split("|")[gene_ann_index].split(":")[0]!=j_list[k].split("|")[gene_ann_index].split(":")[1] and j_list[k].split("|")[gene_ann_index].split(":")[0]==ref_gt:# this is a heterozygotes.
                    case_allele1+=len(continous_ids(j_list[k].split("|")[0].split(",")))
                    snp_id=j_list[0][3:]+"_"+"_".join(j_list[1:3])+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[1]
                    snp_s=j_list[0][3:]+"_"+j_list[1]
                else:
                    pass#print("#none geno type found",j)
        elif re.match("chr",j) and (j_list[3]=="deletion"):
            case_allele0=0
            snp_id,snp_id_s,n_case,case_allele1,case_allele0=generate_score_input(j_list)

        elif re.match("chr",j) and (j_list[3]=="insertion"):
            case_allele0=0
            snp_id,snp_id_s,n_case,case_allele1,case_allele0=generate_score_input(j_list)
        control_allele1, control_allele0=exac_counts(snp_id)# get ExAC counts
        case_allele0=2*int(n_case)-case_allele1
        if case_allele1+case_allele0==0 or control_allele0+control_allele1==0:
            #pass#input_allele=[case_allele1, case_allele0, control_allele1,control_allele0]
            s=1 # set an arbitrial 0 log lr=
            snp_list=snp_id.split("_")
            bed_id="chr"+snp_list[0]+"\t"+snp_list[1]+"\t"+str(int(snp_list[1])+len(snp_list[4])-1)+"\t"+str(s)
            print(bed_id)
            
        else:
            control_allele1, control_allele0=exac_counts(snp_id)# get ExAC counts
            case_allele0=2*int(n_case)-case_allele1
            input_allele=[case_allele1, case_allele0, control_allele1,control_allele0]
            #if snp_id_s=="1":
            #    pass
            #else:
            if snp_id_s!="1":
                casm=0.0#float(casm_data[snp_id_s.split("_")[0]][snp_id_s])
                s=gene_vaast_score(input_allele, casm)
                if s>100:
                    s=100
                s=np.exp(s)
                snp_list=snp_id.split("_")
                bed_id="chr"+snp_list[0]+"\t"+snp_list[1]+"\t"+str(int(snp_list[1])+len(snp_list[4])-1)+"\t"+str(s)
                print(bed_id)
#        print(snp_id,snp_s)
cdr.close()

