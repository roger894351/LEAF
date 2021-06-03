#!/usr/local/python

import re
import sys
import shutil
import numpy as np
import OR_adjust
import bed_reader
from scipy import stats
from optparse import OptionParser
import pandas as pd

#usage: python create_filter_bed.py original_bed LR_file new_bed_file.txt
#python  ~/apps/roger_script/VAAST-CASM/real_data/create_bed_filter.py test_out.txt test.lr0.txt oradjust
#from optparse import OptionParser

#def main():
usage = "usage: %prog [options] arg "

parser = OptionParser(usage)
#parser.add_option("-f", "--file", dest="filename",
#                  help="read data from FILENAME")
parser.add_option("-a", action="store_true", dest="verbose_a", default=False)
#parser.add_option("-s", action="store_true", dest="verbose_s", default=False)
#    parser.add_option("-a", "--accuracy_prior",
#                      action="store_true", dest="verbose",
#                      help="set_accuracy_prior than sensitivity", default=False)
#    parser.add_option("-s", "--sensitivy_prior",
#                      action="store_true", dest="verbose",
#                      help="set sensitivity prior than sensitivity", default=False)
#parser.add_option("-c", "--accuracy_cutoff",
#                  action="store", type="int",dest="num",
#                  help="set the max copy allele of accuracy level", default=0.8)
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose")

(options, args) = parser.parse_args()
    
#    options["verbose"]
#    print(options[1])
#print(options.verbose_a)
    #if len(args) != 1:
    #    parser.error("incorrect number of arguments")
    
if options.verbose:
    print("reading %s..." % options.filename)

#if __name__ == "__main__":
#    main()



#lr_filter=[2,4,6.82]
#or_filter=[1.0,1.2,1.4,2.0,3.0]
#or_filter=[1.0,3.0,4.0,4.5,5.0]# change for different or filter
or_filter=[1.0]
sensitivity_level=sys.argv[4]
remove_nonsense_bed=sys.argv[5]
or_adjust=OR_adjust.read_or_adj(sys.argv[3])# OR_adjust_leve_with_different_cuts

#max_copy_cutoff=0.8 # this is setting a maximum copy for use information.
#
#for i in range(len(or_adjust['max_accuracy'].keys())):
#    i=str(i+1)
#    if float(or_adjust['max_accuracy'][i])<=0.8:
#        pass
#    else:
#       #print(or_adjust['max_accuracy'][i])
#       print(i)
#       max_copy=i
#       break

max_copy=36
k=pd.read_csv(sys.argv[3],header=0, sep="\t", index_col=0)
k=k.loc[1:int(max_copy)]
if options.verbose_a:
    priority_list=[]
    for j in k.columns:
        comparison_column_greater = np.where(k['or_to_maximize_accuracy'] > k[j], 1, 0)
        #print(comparison_column_greater)
        temp_sum_g=np.sum(comparison_column_greater)
        if temp_sum_g<int(max_copy) and temp_sum_g>0:
            #print(j)
            k['low_prior_'+j]=k[['or_to_maximize_accuracy',j]].min(axis=1)
            k['high_prior_'+j]=k[['or_to_maximize_accuracy',j]].max(axis=1)
#            print(k['low_prior_'+j])
#            print(k['high_prior'+j])
            #v=[str(x) for x in range(1,int(max_copy)+1)]
            #print(v)
            #k['low_prior_'+j].reindex(v)
            or_adjust['low_prior_'+j]=k['low_prior_'+j]
            or_adjust['low_prior_'+j].to_dict()
            or_adjust['high_prior_'+j]=k['high_prior_'+j]
            or_adjust['high_prior_'+j].to_dict()
            #print(or_adjust['low_prior_'+j])
            priority_list.append('low_prior_'+j)
            priority_list.append('high_prior_'+j)
            or_adjust[j]=k[j]
            or_adjust[j].to_dict()
#            print(or_adjust.keys())
            #print(or_adjust['low_prior_'+j][1])
        else:
            or_adjust[j]=k[j]
            or_adjust[j].to_dict() 
#            print("else", or_adjust.keys())
            #print(or_adjust[j][1])
    print("available cutoff: "+"\t".join(priority_list))

#            conditions = [k['or_to_maximize_accuracy'] > k[j] ,k['or_to_maximize_accuracy'] == k[j]]#accuracy>sensitivity
#            conditions = [(k['or_to_maximize_accuracy'] >= k[j]) & (k['or_to_maximize_accuracy'] <= k[j]),k['or_to_maximize_accuracy'] < k[j]]#accuracy>sensitivity
#            choices = [k['or_to_maximize_accuracy'], k[j]]
#            choices =k[['or_to_maximize_accuracy',j]].min(axis=1)
#            choices = [k['or_to_maximize_accuracy']]
#            k1 = np.select(conditions, choices, default=k['or_to_maximize_accuracy'])
            #print(k1)
            #print(len(k[j]))
            
#            k2=k[['or_to_maximize_accuracy','test']]
#            print(k2)
#            print(k2.equals(k2))
#	    if k2.equals(k2):
#                print(j)
           # if k[['or_to_maximize_accuracy']]!=k[['test']]:
            #    print(j)

            #if k['test']!=k['or_to_maximize_accuracy']:
            #    print(j)
            #k1=k.select(k[['or_to_maximize_accuracy']]>k[[j]],axis=1)
#            print(k1)
            #print(k['que'])

#else:
#    print(k)

#for i in or_adjust[""]
if remove_nonsense_bed!="n":
    nonsense_bed=bed_reader.exact_bed_reader(sys.argv[5])
#    or1=open(sys.argv[1].split("/")[-1].split(".")[-3]+"_sensitivity_"+str(sensitivity_level)+"ex_nonsense.bed","w")
    or1=open(sys.argv[1].split(".")[-3]+"_sensitivity_"+str(sensitivity_level)+"ex_nonsense.bed","w")
else:
#    or1=open(sys.argv[1].split("/")[-1].split(".")[-3]+"_sensitivity_"+str(sensitivity_level)+".bed","w")
    or1=open(sys.argv[1].split(".")[-3]+"_sensitivity_"+str(sensitivity_level)+".bed","w")
LR_file=open(sys.argv[2],"r")

for k in LR_file:
    sign=0
    k=k.strip()
    k_list=k.split()
    or_value=float(k_list[7])
#    print(or_value)
    lr_value=float(k_list[6])
    int_alt=float(k_list[8])#this is the alternative allele count in the internal data
    int_ref=float(k_list[9])
    ext_alt=float(k_list[10])
    ext_ref=float(k_list[11])
    try:
        ext_freq=1.0*ext_alt/(ext_alt+ext_ref)
    except ZeroDivisionError:
        ext_freq=1# if external observation is much higher => batch effect
# use binomial
    if int_alt==0:
#        binomial_test_p=0
        fisher_test_p=0
    else:
#        binomial_test_p=stats.binom_test(int_alt,n=int_alt+int_ref, p=ext_freq, alternative='less')
    #    print(binomial_test_p)
        fisher_test_p=stats.fisher_exact([[int_alt, int_ref], [ext_alt, ext_ref]], alternative='less')[1]
#    if binomial_test_p<0.01:
#or_adjust['max_accuracy'][int_alt]
#    if fisher_test_p<0.1:
    ###dynamic cutoffs based on accuracy(too loose):
    #try:
    #    accuracy=float(or_adjust['max_accuracy'][int_alt])
    #except KeyError:
    #    accuracy=1.0
    #if fisher_test_p<1.0-accuracy:
    #    pass
#    if int_alt==1 and fisher_test_p<0.3:
#        pass
#    elif int_alt==2 and fisher_test_p<0.2:
#        pass
#    elif int_alt>=3 and int_alt<int(max_copy) and  fisher_test_p<0.1:
#        pass
#    elif int(int_alt)<int(max_copy):
    if int_alt>=1 and int(int_alt)<int(max_copy) and fisher_test_p<0.1:
        pass
    elif int_alt>=1 and int(int_alt)<int(max_copy) and fisher_test_p>=0.1:
#    else:
#    if int(int_alt)<int(max_copy):
#        print(int_alt)
        int_alt=int(k_list[8])
        if remove_nonsense_bed!="n":
            variant_id=k_list[0]+"-"+"-".join(k_list[1:3])
            if variant_id in nonsense_bed[k_list[0][3:]]:
                #print("here is the nonsense", variant_id)
                pass
            else:
                try:
                    OR_factor=or_adjust[sensitivity_level][int_alt]
                except KeyError:
                    OR_factor=1.0
                if or_value<or_filter[0]*OR_factor and or_value>=0:
                    or1.write("\t".join(k_list)+"\t"+str(or_filter[0]*OR_factor)+"\t"+str(fisher_test_p)+"\n")
        else:
            try:
                OR_factor=or_adjust[sensitivity_level][int_alt]
            except KeyError:
                OR_factor=1.0
            if or_value<or_filter[0]*OR_factor and or_value>=0:
                or1.write("\t".join(k_list)+"\t"+str(or_filter[0]*OR_factor)+"\t"+str(fisher_test_p)+"\n") 
    else:
        pass

or1.close()

