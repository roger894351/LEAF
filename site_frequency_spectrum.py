#!/usr/local/python

import sys
import re
import pandas as pd
from subprocess import Popen, PIPE

a=pd.read_table(sys.argv[1], header=0, dtype={"#CHROM": "string"})
annoted_a=sys.argv[2]

Popen(['sh','prepare_non_synonymous.sh',sys.argv[1], sys.argv[2]],shell=False, stderr=PIPE, stdout=PIPE)

#a=pd.read_table('all_exome_non_cancer_af_nonsynonymous_nononsense_uniq.bed', header=0, dtype={"#CHROM": "string"})
out=sys.argv[1].split(".")[0]+"_out.txt"
print(a.describe())
a.drop_duplicates(keep="first",inplace=True)
b=a.iloc[:,6].value_counts()
c=b.sum()
print("sum=",c)
b = pd.DataFrame(b).reset_index()
b=b.drop([0])
c=b.sum()
print("exclude 0 sum=",c)
b.columns = ['copy', 'count']
b.to_csv(out,index=False, sep='\t')  
