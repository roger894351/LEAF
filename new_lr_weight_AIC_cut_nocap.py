#!/opt/annaconda/python3 python
# coding: utf-8

from __future__ import division
import re,time
import numpy as np
import os,sys
import json
import casm_reader
import OR_adjust
import emplr
from subprocess import Popen, PIPE
from scipy import stats
#from __future__ import division

#note: 1/3/2020 print all likelihood scores and case_allele1, control_allele1, case control0, control0 counts.
#start = time.time()
#casmlr_cap=10, lr_cap=np.log(67.30135)
#singleton_factor=2.26
#doubleton_factor=1.64
#read_or_adj.read_or_adj(sys.argv[])


def gene_vaast_score(input_allele, casm= None, score_only_minor = True,
                     score_only_alternative = True, af_cut=0.1, site_penalty=2.0,casmlr_cap=8,lr_cap=100, casm_cap=np.log(67.30135), sensitivity_level=sys.argv[7]):#log(67.30135)=4.20918,log(27.6429358646029):casmlr:8 lr_cap5
    case_allele1,case_allele0,control_allele1,control_allele0=input_allele[0:4]
#    log_lh = lrt_score(case_allele1, control_allele1,## old lr calculation
#            case_allele0, control_allele0)
    if case_allele1==0:
        log_lh=1#original=0
    elif case_allele1>case_allele0:# allele frequency >0.5
        lr=lr_cap
        log_lh=1
    elif case_allele1==case_allele1+case_allele0: #solve overflow case_allele=0
        log_lh=1#lr_cap # get rid of log(nan error in emplr)
    elif case_allele1+case_allele0==0:
        log_lh=1#original=0
    elif control_allele1+control_allele0==0:
        log_lh=1#original=0
#chrX    2700157 2700157 SNV     4.18764462327   0.867282265337  4.82846795172   1.06394393618   1086    6364    10132   63149 over 
    elif case_allele1>max_copy:
        lr=np.exp(1)
        log_lh=1
    else:
        try:
            ext_freq=1.0*control_allele1/(control_allele1+control_allele0)
        except ZeroDivisionError:
            ext_freq=1# if external observation is much higher => batch effect
        if case_allele1==0:
            binomial_test_p=0
        else:
            binomial_test_p=stats.binom_test(case_allele1,n=case_allele1+case_allele0, p=ext_freq, alternative='less')
#        print(binomial_test_p)
        if binomial_test_p<0.1:
            lr=np.exp(1)
            log_lh=1
        else:
            lr=emplr.get_lr(case_allele1,control_allele1,case_allele1+case_allele0,control_allele1+control_allele0,2)
#            if lr>lr_cap:
#                lr=lr_cap
#            if lr<np.exp(1):# trunck the protect is not woring well
#                lr=np.exp(1)            
#            else:
#                pass
            log_lh=np.log(lr)
#        print(lr)
#        log_lh=np.log(emplr.get_lr(case_allele1,control_allele1,case_allele1+case_allele0,control_allele1+control_allele0,2))
#    print(log_lh,input_allele)
    if casm==0.0 or casm==0:
#        log_casm=0.0-casm_cap
        log_casm=0.0-casm_cap
    else:
        log_casm=np.log(casm)
  
#    if log_lh=="inf": #remove
#        log_lh=0.0-lr_cap #remove
    vaast_site_scores =  2.0 * (log_lh + log_casm) - site_penalty
    mask = (vaast_site_scores<=0)
    
    # mask sites where major allele is more common in cases
    if score_only_minor:
        # scenario 1: 1 is minor allele and it's more frequent in case
        scenario1= (    (control_allele1 <= control_allele0 )
                    &   (case_allele1 * control_allele0 >= case_allele0 *
                        control_allele1)
                    )
        scenario2= (    (control_allele1 >= control_allele0 )
                    &   (case_allele1 * control_allele0 <= case_allele0 *
                        control_allele1)
                    )
        # mask if previous or new masking condition is satisfied
        if score_only_alternative:
            # if we only want to score alternative allele, then mask
            # unless 1 is minor allele and it is more frequent in case
            mask = (mask | (scenario1==False) )
        else:
            # otherwise we accept both scenario 1 and 2
            mask = (mask | ((scenario1==False) & (scenario2==False)))

#    combined_alt_freq=float(case_allele1+control_allele1)/(control_allele1+control_allele0+case_allele1+case_allele0) external counts will be added
    if case_allele1+case_allele0==0.0:
        combined_alt_freq=0.0
    else:
        combined_alt_freq=float(case_allele1)/(case_allele1+case_allele0)
    if control_allele1==0 or case_allele0==0:
        OR=1.0*(case_allele1+0.5)*(control_allele0+0.5)/((1.0*(control_allele1+0.5))*(case_allele0+0.5))
    else:
        OR=1.0*((case_allele1*control_allele0*1.0)/(control_allele1*case_allele0*1.0))
        #print(OR)
        # 
        #print(case_allele1,control_allele0,control_allele1,case_allele0)
###adjust for singletons and doubletons:
#    if case_allele1==1:
#        OR_factor=singleton_factor#8.24#2.1# 8.24 too strong these value are based on synonymous odds ratio 
#    elif case_allele1==2:
#        OR_factor=doubleton_factor#3.63#1.6# 1.64 too strong
#    else:
#        OR_factor=1.0
    if float(sensitivity_level)>1:
        sensitivity_level=round(float(sensitivity_level)/100,3)
        if sensitivity not in or_adjust.keys():
            print("please provide valid sensitivity cutoffs from the following selection"+"\t".join(or_adjust.keys())+'.') 
#    try:
#        OR_factor=or_adjust[sensitivity_level][case_allele1]
#    except KeyError:
#        OR_factor=1.0
    #print(OR,OR_factor)
    OR_sign=1.0 # release this block
#    if OR<=1.0*OR_factor:
#        OR_sign=-1*OR_sign
    #if log_casm<log_lh:#try to get rid of biased likelihood:
    #    log_lh=0
####    if log_lh==0.0:
#        log_lh_e=0.0
#    elif OR_sign*log_lh>=lr_cap:
#        log_lh=lr_cap
#        log_lh_e=np.exp(log_lh)
#    elif OR_sign*log_lh<lr_cap and OR_sign*log_lh>0:
#    log_lh_e=np.exp(OR_sign*log_lh)
#    else:
    log_lh_e=np.exp(OR_sign*log_lh)
#    total_score = ((OR_sign)*log_lh)+log_casm
    if log_lh==1 and log_casm==0.0-casm_cap:
        total_score=0
    else:
        total_score = log_lh+log_casm
#    print('log_lh,casm,log_casm,total_score',log_lh,casm,log_casm,total_score)
#    total_score=total_score# for AIC correction 
    if total_score==0.0: total_score_e=0.0#print0 for total_score and log_lh
#    elif total_score>=casmlr_cap: total_score_e=np.exp(casmlr_cap)
#    elif total_score<=0: total_score_e=0
#    elif casm==0: 
#        total_score_e=0
    else:
        total_score_e = np.exp(total_score)
    #return([str(total_score_e),str(casm), str(log_lh_e),str(OR)])
    return([str('{:3f}'.format(total_score_e)),str(casm), str('{:3f}'.format(log_lh_e)),str(OR)])
def lrt_score(case_allele1, control_allele1, case_allele0, control_allele0):
    ''' Calculate the unpenalized log likelihood of a gene 
        INPUT: allele counts in cases and controls
        OUTPUT: site lrt scores
    '''
    # clip
        ## adjust for singleton and doubleton biases:
    if case_allele1==1:
        control_allele0=control_allele0-(singleton_factor-1)*control_allele1
        control_allele1=control_allele1*singleton_factor
    elif case_allele1==2:
        control_allele0=control_allele0-(doubleton_factor-1)*control_allele1
        control_allele1=control_allele1*doubleton_factor
    else:pass# for now later we can adjust 1.3

    if control_allele0+control_allele1==0:
#        alt_control_freq=1e-9
        return("inf")
    else:
        alt_control_freq = 1.0 * control_allele1 / (control_allele0 + control_allele1)
    if case_allele0+case_allele1==0:
#        alt_case_freq=1e-9
        return("inf")
    else:
        alt_case_freq = 1.0 * case_allele1 / (case_allele0 + case_allele1)
    if control_allele0+control_allele1+case_allele0+case_allele1==0:
#        null_freq = 1e-9
        return("inf")
    else:
        null_freq = 1.0 * (case_allele1 + control_allele1) / (control_allele0 +
                control_allele1 + case_allele0 + case_allele1)
    alt_log_lh = log_likelihood(alt_control_freq, control_allele0, control_allele1
                    ) + log_likelihood(alt_case_freq, case_allele0, case_allele1)
    null_log_lh= log_likelihood(null_freq, control_allele0+case_allele0, control_allele1+case_allele1)
    k=alt_log_lh - null_log_lh
    return(k)

def log_likelihood(freq, allele0, allele1):
    ''' Given specified frequencies, allele0 count, allele1 count, calculate
        log likelihood
        INPUT: specified frequencies, allele0 count, allele1 count
        OUTPUT: log likelihood from binomial distribution
    '''

    # bound frequency estimate to [1e-9< 1.0-1e-9] to avoid numerical issue
    freq = np.clip(freq, 1e-9, 1.0-1e-9)

    # calculate log likelihood
    # to avoid alternative allel is a major allele
    if allele1==0:
#        loglik = allele0 * np.log(1.0-freq)
       loglik=1e-9 * np.log(1e-9) + allele0 * np.log(1.0-freq)
    else:
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

json1=sys.argv[2]
json1_file =open(json1)
json1_str = json1_file.read()
json1_data = json.loads(json1_str) ##create a json_dictionary

def get_cdr_n(f): #give a cdr and return the last sample number is the cdr.
    #f = 'combined_case_control.cdr'
    # Get the last line from the file
    p = Popen(['tail','-1',f],shell=False, stderr=PIPE, stdout=PIPE)
#    p = Popen(['grep','FILE-INDEX',f],shell=False, stderr=PIPE, stdout=PIPE)
    number,err = p.communicate()
    if err:
        print (err.decode())
    else:
        # Use split to get the sample number in the cdr
        number = int(number.decode().split()[2].strip())+1 
#        number= int(number.decode().split("#")[-1].split('\t')[2])+1#modify for the bootstrap run
    return(number)

def exac_counts(x): #x =snp_id with genotype info fro exac id ex 13	32912223	32912223	T	C =13_32912223_T_A
    if x=="0":
        return(0,0)
    else:
        chr_n=x.split("_")[0]
        try: # anoter way to process KeyError issue but the time is a little bit longer
            value=json1_data[chr_n][x]
            alt,total=int(value.split(",")[0]),int(value.split(",")[1])
            ref=total-alt
            return(alt,ref)
        except KeyError:
########## a block to scan if no matched id was found:
#            ref_match="_".join(x.split("_")[:-1])# If a genotype alternative is not found will check other alternative genotype in Exac as external
#            if json1_data[chr_n].keys():
#                for key in json1_data[chr_n]:
#                    if re.match(ref_match, key):
#                        value=json1_data[chr_n][key]
#                        alt,ref=int(value.split(",")[0]),int(value.split(",")[1])
#                       print(key)
##                      print(json1_data[chr_n][key])
#                        return(alt,ref)
#            #elif len(json1_data[chr_n].keys()==0):
            return(0,ext_average_count)
def casm_value(x):# x=snp_id_s
    spda_check="0"
    if x=="0":
        return(0.0)
    else:
       chr_n=x.split("_")[0]
       try:
           casm_value=float(casm_dict[chr_n].pop(x))
       except KeyError:
#          # print(x)
           test_SPDA="_".join(x.split("_")[0:3])+"_SPDA"
           try:
               spda_check="1"
               casm_value=float(casm_dict[chr_n].pop(test_SPDA))
#               if casm_value:
#                   spda_check="1"
           except KeyError:
               spda_check="0"
               casm_value=0.0
       return(casm_value,spda_check)    
   
def generate_score_input(j_list):
    ref_gt=j_list[5].split("|")[0]# whatever will be true
    #print(ref_gt)
    n_case=n_cdr
    case_allele1=0# set 0 first; #
    case_allele0=2*int(n_case);# get these numbers out:control_allele1=0;control_allele0=0
    gene_ann_index=1#if len(j_list[-1].split("|"))==3: 
    snp_id="1"# this is snp_id reset
    snp_id_s="1"
    for k in range(6,len(j_list)):
        if j_list[k].split("|")[gene_ann_index]=="^:^":## nocall: no genotype
            miss_count=2*len(continous_ids(j_list[k].split("|")[0].split(",")))##I change here
            case_allele0-=miss_count
            if snp_id=="1": # control the flow if meet no_call first
                snp_id=j_list[0][3:]+"_"+j_list[1]+"_"+ref_gt+"_"+"?"#+j_list[k].split("|")[gene_ann_index].split(":")[0] #nocalls if no other genotypes
                snp_id_s=j_list[0][3:]+"_"+"_".join(j_list[1:4])           
        elif j_list[k].split("|")[gene_ann_index].split(":")[0]==j_list[k].split("|")[gene_ann_index].split(":")[1] and j_list[k].split("|")[gene_ann_index].split(":")[0]!=ref_gt:# this is a homozygus hit
            if j_list[k].split("|")[gene_ann_index].split(":")[0]=="!":
                snp_id=j_list[0][3:]+"_"+j_list[1]+"_"+ref_gt+"_"+"-" 
                snp_id_s=j_list[0][3:]+"_"+"_".join(j_list[1:4])
            else:
                snp_id=j_list[0][3:]+"_"+j_list[1]+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[0]
#            snp_id_s=j_list[0][3:]+"_"+j_list[1]#homozygotes
                snp_id_s=j_list[0][3:]+"_"+"_".join(j_list[1:4])
            case_allele1+=2*(len(continous_ids(j_list[k].split("|")[0].split(","))))
        elif j_list[k].split("|")[gene_ann_index].split(":")[0]!=j_list[k].split("|")[gene_ann_index].split(":")[1] and ( j_list[k].split("|")[gene_ann_index].split(":")[0]==ref_gt or j_list[k].split("|")[gene_ann_index].split(":")[1]==ref_gt):# this is a heterozygotes.
            snp_id=j_list[0][3:]+"_"+j_list[1]+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[1]
            snp_id_s=j_list[0][3:]+"_"+"_".join(j_list[1:4])
            case_allele1+=len(continous_ids(j_list[k].split("|")[0].split(",")))
        elif j_list[k].split("|")[gene_ann_index].split(":")[0]!=j_list[k].split("|")[gene_ann_index].split(":")[1] and j_list[k].split("|")[gene_ann_index].split(":")[1]!=ref_gt:# this is a secondary heterozygotes that is different than original one.
            snp_id=j_list[0][3:]+"_"+j_list[1]+"_"+ref_gt+"_"+j_list[k].split("|")[gene_ann_index].split(":")[1]# change to another genotype in the future now combined to the same snp_id
            snp_id_s=j_list[0][3:]+"_"+"_".join(j_list[1:4])
            case_allele1+=len(continous_ids(j_list[k].split("|")[0].split(",")))
            #print("the intronic variants")    
        else:
            pass#
            print("#more than 2 genotypes found",j_list[k].split("|")[gene_ann_index].split(":")[0],j_list[k].split("|")[gene_ann_index].split(":")[1],ref_gt)
#        print(snp_id,snp_id_s,n_case,case_allele1,case_allele0)
    return(snp_id,snp_id_s,n_case,case_allele1,case_allele0)

#--main--
#if __name__ == "__main__":
f=sys.argv[1]# please input a AF file here
cdr_bed=open(sys.argv[4],"w")
n_cdr=get_cdr_n(f) #get number of individual in cdr
casm_file=sys.argv[3]
casm_dict=casm_reader.casm_reader(casm_file)
ext_average_count_data=sys.argv[5]
or_adjust=OR_adjust.read_or_adj(sys.argv[6])
#print(sys.argv[6])
sensitivity_level=sys.argv[7]
emplr.init(sys.argv[8])
ext_average_d={"all_gnomad":96425,"gnomad_controls":5338, "exome_controls": 37016, "exome_non_cancer":89262,"exome_non_cancer_nfe_nwe":37121 }
ext_average_count=ext_average_d[ext_average_count_data]#exome_non_cancer

#set_maximum accuracy to use external frequnecy
# OR_adjust_leve_with_different_cuts
#max_copy_cutoff=0.8 # this is setting a maximum copy for use information.
#
#for i in range(len(or_adjust['max_accuracy'].keys())):
#    i=str(i+1)
#    if float(or_adjust['max_accuracy'][i])<=0.8:
#        pass
#    else:
#       #print(or_adjust['max_accuracy'][i])
#       #print(i)
#       max_copy=int(i)
#       print("max_copy",i)
#       break

max_copy=36
with open(f,"r") as cdr:
    for j in cdr: # this f is handeling a cdr file noted a vcf file
        j=j.strip()
        j_list=j.split("\t")
        if re.match("#",j_list[0]):#most index
            snp_id='0'
            snp_id_s="1"
        elif re.match("chr",j) and (j_list[3]=="SNV" or j_list[3]=="complex_substitution" ):
            ref_gt=j_list[5].split("|")[0]# whatever will be true
            n_case=n_cdr# delete no call individuals
            case_allele1=0# set 0 first; #case_allele0=2*int(n_case);# get these numbers out:control_allele1=0;control_allele0=0
            gene_ann_index=1#if len(j_list[-1].split("|"))==3: 
            snp_id,snp_id_s,n_case,case_allele1,case_allele0=generate_score_input(j_list)
        elif re.match("chr",j) and (j_list[3]=="deletion"):# I will have to modeify this line because deletion the start should minus the number of delected basepairs.
            snp_id,snp_id_s,n_case,case_allele1,case_allele0=generate_score_input(j_list)
            #print(snp_id, snp_id_s, "deletion")
        elif re.match("chr",j) and (j_list[3]=="insertion"): #The insertion in cdr file is the same as vcf file.
            snp_id,snp_id_s,n_case,case_allele1,case_allele0=generate_score_input(j_list)
            #print(snp_id, snp_id_s, "insertion")
        control_allele1, control_allele0=exac_counts(snp_id)# get ExAC counts
        case_allele0=case_allele0-case_allele1
        input_allele=[case_allele1, case_allele0, control_allele1,control_allele0]
        input_allele2=[str(case_allele1), str(case_allele0), str(control_allele1),str(control_allele0)]
        if snp_id_s=="1" and snp_id!="0": # use weighted approach all variant are required a weight.
            casm,spda_check=casm_value(snp_id_s)#0.0#float(casm_data[snp_id_s.split("_")[0]][snp_id_s])
            s=gene_vaast_score(input_allele, casm)
            snp_list=snp_id.split("_")
            if spda_check=="1":
                bed_id="chr"+snp_list[0]+"\t"+snp_list[1]+"\t"+str(int(snp_list[1])+len(snp_list[2])-1)+"\t"+"SPDA"+"\t"+"\t".join(s)+"\t"+"\t".join(input_allele2)
            else:

                bed_id="chr"+snp_list[0]+"\t"+snp_list[1]+"\t"+str(int(snp_list[1])+len(snp_list[2])-1)+"\t"+j_list[3]+"\t"+"\t".join(s)+"\t"+"\t".join(input_allele2)
            cdr_bed.write(bed_id+"\n")
        elif snp_id=='0': # File-Index Line will be here
            pass#print(snp_id,"0")
        else:
            casm,spda_check=casm_value(snp_id_s)#0.0#float(casm_data[snp_id_s.split("_")[0]][snp_id_s])
            s=gene_vaast_score(input_allele, casm)
            snp_list=snp_id.split("_")
            if spda_check=="1":
                bed_id="chr"+snp_list[0]+"\t"+snp_list[1]+"\t"+str(int(snp_list[1])+len(snp_list[2])-1)+"\t"+"SPDA"+"\t"+"\t".join(s)+"\t"+"\t".join(input_allele2)
            else:
                bed_id="chr"+snp_list[0]+"\t"+snp_list[1]+"\t"+str(int(snp_list[1])+len(snp_list[2])-1)+"\t"+j_list[3]+"\t"+"\t".join(s)+"\t"+"\t".join(input_allele2)
            cdr_bed.write(bed_id+"\n")
cdr.close()
