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

usage = "usage: %prog [options] arg "

parser = OptionParser(usage)
parser.add_option("-a", action="store_true", dest="verbose_a", default=False)
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose")

(options, args) = parser.parse_args()
    
    
if options.verbose:
    print("reading %s..." % options.filename)

or_filter=[1.0]
sensitivity_level=sys.argv[4]
remove_nonsense_bed=sys.argv[5]
or_adjust=OR_adjust.read_or_adj(sys.argv[3])# OR_adjust_leve_with_different_cuts

max_copy=36
k=pd.read_csv(sys.argv[3],header=0, sep="\t", index_col=0)
k=k.loc[1:int(max_copy)]
if options.verbose_a:
    priority_list=[]
    for j in k.columns:
        comparison_column_greater = np.where(k['or_to_maximize_accuracy'] > k[j], 1, 0)
        temp_sum_g=np.sum(comparison_column_greater)
        if temp_sum_g<int(max_copy) and temp_sum_g>0:
            k['low_prior_'+j]=k[['or_to_maximize_accuracy',j]].min(axis=1)
            k['high_prior_'+j]=k[['or_to_maximize_accuracy',j]].max(axis=1)
            or_adjust['low_prior_'+j]=k['low_prior_'+j]
            or_adjust['low_prior_'+j].to_dict()
            or_adjust['high_prior_'+j]=k['high_prior_'+j]
            or_adjust['high_prior_'+j].to_dict()
            priority_list.append('low_prior_'+j)
            priority_list.append('high_prior_'+j)
            or_adjust[j]=k[j]
            or_adjust[j].to_dict()
        else:
            or_adjust[j]=k[j]
            or_adjust[j].to_dict() 
    print("available cutoff: "+"\t".join(priority_list))

if remove_nonsense_bed!="n":
    nonsense_bed=bed_reader.exact_bed_reader(sys.argv[5])
    or1=open(sys.argv[1].split(".")[-3]+"_sensitivity_"+str(sensitivity_level)+"ex_nonsense.bed","w")
else:
    or1=open(sys.argv[1].split(".")[-3]+"_sensitivity_"+str(sensitivity_level)+".bed","w")
LR_file=open(sys.argv[2],"r")

for k in LR_file:
    sign=0
    k=k.strip()
    k_list=k.split()
    or_value=float(k_list[7])
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
        fisher_test_p=stats.fisher_exact([[int_alt, int_ref], [ext_alt, ext_ref]], alternative='less')[1]
    if int_alt>=1 and int(int_alt)<int(max_copy) and fisher_test_p<0.1:
        pass
    elif int_alt>=1 and int(int_alt)<int(max_copy) and fisher_test_p>=0.1:
        int_alt=int(k_list[8])
        if remove_nonsense_bed!="n":
            variant_id=k_list[0]+"-"+"-".join(k_list[1:3])
            if variant_id in nonsense_bed[k_list[0][3:]]:
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

