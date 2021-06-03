#!/usr/local/python

import sys
import re
import pandas as pd

#with open(sys.argv[1],'r') as a:
#    sum={}
#    for i in a:
#        if re("#",i):
#            pass
#        else:
#            i=i.strip()
#            i_list=i_list.split()
a=pd.read_table(sys.argv[1], header=0, dtype={"#CHROM": "string"})
#out=open(sys.argv[1].split(".")[0]+"out.txt","w")
out=sys.argv[1].split(".")[0]+"_out.txt"
print(a.describe())
a.drop_duplicates(keep="first",inplace=True)
b=a.iloc[:,6].value_counts()
#b=a["non_cancer_AC_alt"].value_counts()
c=b.sum()
print("sum=",c)
b = pd.DataFrame(b).reset_index()
#print(a.describe())
b=b.drop([0])
c=b.sum()
print("exclude 0 sum=",c)
b.columns = ['copy', 'count']
b.to_csv(out,index=False, sep='\t')  
