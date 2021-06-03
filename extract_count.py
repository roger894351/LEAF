#!/opt/annaconda3/python

#description:
#usage: this script yield the gnomAD information from 
# python extract_count.py ../vcf/exomes/head.vcf test.txt NFE
#CHROM	POS	ID	REF	ALT	AC_ref	AC_NFE	AN_NFE	AF

import sys
import os
import re

#def counts(x)

#script:python extract_count.py ../vcf/exomes/head.vcf test.txt NFE

pop=sys.argv[3] #thsi can be "NFE, AFR, AMR, ASJ, EAS,FIN. OTH"
out=open(sys.argv[2],'w')
with open(sys.argv[1],'r') as vcf:

    for line in vcf:
        if re.match("##",line):
            pass
        elif re.match("#CHROM",line):
            out.write("\t".join(['#CHROM','POS','ID','REF','ALT','AC_ref','AC_NFE','AN_NFE','AF'])+"\n")
        else:
            line=line.strip()
            line_list=line.split()
            if re.search(",",line_list[4]):
                #for i in line_list:
                genotypes=line_list[4].split(",")
                d_len=len(genotypes)
                #print(d_len)
                #dict_list={}*d_len
                #for j in dict_list:
                d={}
                i_list=line_list[7].split(";")
#                print(i_list) 
                #j_list=j.split("=")
                for j in i_list:
                    #print(len(i_list))
                    j_list=j.split("=")
                  
                    #print(j_list)#if len(j_list)>=2:
                    if len(j_list)==1:
                        pass
#                    elif re.search(",",j_list[1]) and len(j_list)>2:
#                        value, key=j_list[1],j_list[0]
#                        d[key]=value.split(",")
#                    #    info=line_list[0:4]
                    elif len(j_list)>=2: 
                        value, key=j_list[1],j_list[0]
                        #d[key]=value
                        d[key]=value.split(",")
                    else: pass
        #        print("K")
                    #info=list(line_list[0]+line_list(1)+line_list)
#                print("K")    
                for k in range(0,d_len):
                    #print(k)
                    if d["AC_"+pop][k]==".":
                        #info_d["AC_"+pop][k]=0;info_d["AN_"+pop][k]=0; info_d("AF_"+pop)[k]=0
                        d["AC_"+pop][k]=0;d["AN_"+pop]=0; d("AF_"+pop)[k]=0
                        #no_data=["0","0","0"]
                        #info=info+nodata
                        
                    else:
                    #    print(d["AN_"+pop])
                    #    print(d["AC_"+pop][k])
                        info_ref=str(int(d["AN_"+pop][0])-int(d["AC_"+pop][k]))#print the AN coun#print the AN countt
#info_d["AC_"+pop]=int(info_d["AN_"+pop]-INT(info_d["AC_"+pop]))
                    tri=line_list[0:4]
                    tri.append(genotypes[k])
                    tri.append(info_ref)
                    tri.append(d["AC_"+pop][k])
                    tri.append(d["AN_"+pop][0])
                    tri.append(d["AF_"+pop][k])
                    out.write("\t".join(tri)+"\n")


            else:
              #  for i in line_list:
                d={}
                i_list=line_list[7].split(";")
                for j in i_list:
                        #print(j)
                    j_list=j.split("=")
                    if len(j_list)==2:
                        value, key=j_list[1],j_list[0]
                        d[key]=value                
                        info=line_list[0:5]
                    else: pass#print(j_list)
                #print(d)
                #info=list(line_list[0]+line_list(1)+line_list)
                if d["AC_"+pop]=="." or d["AF_"+pop]==".":
                    d["AC_"+pop]="0";d["AN_"+pop]="0"; d["AF_"+pop]="0"
                    #no_data=["0","0","0"]
                    info_ref="0"
                else:
                    info_ref=str(int(d["AN_"+pop])-int(d["AC_"+pop]))
#info_d["AC_"+pop]=int(info_d["AN_"+pop]-INT(info_d["AC_"+pop]))
                info.append(info_ref)
                info.append(d["AC_"+pop])
                info.append(d["AN_"+pop])
                info.append(d["AF_"+pop])
                out.write("\t".join(info)+"\n")
vcf.close()
out.close()

