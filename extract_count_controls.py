#!/opt/annaconda3/python

#description:
#usage: this script yield the gnomAD information from 
# python extract_count.py ../vcf/exomes/head.vcf test.txt nfe
#CHROM	POS	ID	REF	ALT	controls_AC_ref	controls_AC_nfe	controls_AN_nfe	controls_AF

import sys
import os
import re

#def counts(x)

#script:python extract_count.py <vcf file> <output file> <ancestry_group> <non_cancer,controls> nfe non_cancer
prefix=sys.argv[4] #prefix
if prefix not in ["non_cancer", "controls"]:
    print(" please provide subgroups as non_cancer or controls")
    exit()
pop=sys.argv[3] #this can be "nfe, controls_AFR, AMR, ASJ, EAS,FIN. OTH nfe_nwe"
if pop not in ['asj_female','oth_female','female','nfe_onf','male','sas_female','amr_male','fin_male','nfe_seu','eas','amr_female','nfe_female,asj' \
,'eas_jpn','afr_male','afr','amr','nfe_male','oth','oth_male','afr_female','nfe_nwe','nfe_bgr','eas_female','asj_male','nfe_est','nfe','nfe_swe','fin_female','fin,eas_oea' \
,'sas_male','raw,eas_male','sas','eas_kor,popmax']:
    print("please use the ancestry information form gnomAD such as nfe, sas...etc")
    exit()
out=open(sys.argv[2],'w')
with open(sys.argv[1],'r') as vcf:

    for line in vcf:
        if re.match("##",line):
            pass
        elif re.match("#CHROM",line):
            out.write("\t".join(['#CHROM','START','END','REF','ALT',prefix+'_AC_ref'+pop,prefix+'_AC_alt'+pop,prefix+'_AN_'+pop,prefix+'_AF_'+pop])+"\n")
        else:
            line=line.strip()
            line_list=line.split("\t")
            if len(line_list[3])==len(line_list[4]) and len(line_list[3])==1:# process variant start and stop location info to cdr format.
                info=line_list[0:2]+line_list[1:2]+line_list[3:5]
            elif len(line_list[3])>len(line_list[4]): # this is a deletion
                info=line_list[0:1]+[str(int(line_list[1])+len(line_list[4])),str(int(line_list[1])+len(line_list[3])-len(line_list[4]))]+[line_list[3][len(line_list[4]):],"-"]
               # print(info)
            elif len(line_list[3])<len(line_list[4]): # this is a insertion
                info=line_list[0:2]+line_list[1:2]+['-',line_list[4][len(line_list[3]):]]
            else:
                print(line_list[0:6],"break")
                break
            if re.search(",",line_list[4]):
                genotypes=line_list[4].split(",") #more than one genotypes
                d_len=len(genotypes)
                d={}
                i_list=line_list[7].split(";")
                for j in i_list:
                    #print(len(i_list))
                    j_list=j.split("=") # dictionary of every features
                  
                    #print(j_list)#if len(j_list)>=2:
                    if len(j_list)==1:
                        pass
                    elif len(j_list)>=2: 
                        value, key=j_list[1],j_list[0]
                        #d[key]=value
                        d[key]=value.split(",")
                    else: pass
        #dd        print("K")
                for k in range(0,d_len):
                    #print(k)
                    if d[prefix+"_AC_"+pop][k]==".":
                        #info_d[prefix+"_AC_"+pop][k]=0;info_d[prefix+"_AN_"+pop][k]=0; info_d(prefix+"_AF_"+pop)[k]=0
                        d[prefix+"_AC_"+pop][k]=0;d[prefix+"_AN_"+pop]=0; d(prefix+"_AF_"+pop)[k]=0
                        #no_data=["0","0","0"]
                        #info=info+nodata
                        info_ref=str(0)
                        
                    else:
                        info_ref=str(int(d[prefix+"_AN_"+pop][0])-int(d[prefix+"_AC_"+pop][k]))#print the controls_AN coun#print the controls_AN countt
#info_d[prefix+"_AC_"+pop]=int(info_d[prefix+"_AN_"+pop]-INT(info_d[prefix+"_AC_"+pop]))
                    tri=line_list[0:4]
                    tri.append(genotypes[k])
                    tri.append(info_ref)
                    tri.append(d[prefix+"_AC_"+pop][k])
                    tri.append(d[prefix+"_AN_"+pop][0])
                    tri.append(d[prefix+"_AF_"+pop][k])
                    out.write("\t".join(tri)+"\n")


            else: # regular two genotypes:
                d={}
                i_list=line_list[7].split(";")
                for j in i_list:
                        #print(j)
                    j_list=j.split("=")
                    if len(j_list)==2:
                        value, key=j_list[1],j_list[0]
                        d[key]=value                
                    else: info_ref="?"
                try:
                
                    if d[prefix+"_AC_"+pop]=="." or d[prefix+"_AF_"+pop]==".":
                        d[prefix+"_AC_"+pop]="0";d[prefix+"_AN_"+pop]="0"; d[prefix+"_AF_"+pop]="0"
                        #no_data=["0","0","0"]
                        info_re="0"
                    else:
                        info_ref=str(int(d[prefix+"_AN_"+pop])-int(d[prefix+"_AC_"+pop]))
                except KeyError:
                    d[prefix+"_AC_"+pop]="0";d[prefix+"_AN_"+pop]="0"; d[prefix+"_AF_"+pop]="0"    
                if info_ref=="?":
                    pass
                else:
                    info.append(info_ref)
                    info.append(d[prefix+"_AC_"+pop])
                    info.append(d[prefix+"_AN_"+pop])
                    info.append(d[prefix+"_AF_"+pop])
                    out.write("\t".join(info)+"\n")
vcf.close()
out.close()

