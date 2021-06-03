#!/usr/local/python
#grep nonsnse variant: stops, indel
#note input is a cdr file

import sys
import re

cdr=open(sys.argv[1],'r')
for i in cdr:
    i=i.strip()
    if re.match("#",i):
        pass
    else:
        i_list=i.split("\t")
        if i_list[3]=='deletion':
            print("\t".join(i_list[0:3]))
        elif i_list[3]=='insertion':
            print("\t".join(i_list[0:3]))
        elif i_list[3]=='SNV':
            if re.search('stop',i_list[4]):
                print("\t".join(i_list[0:3]))
        elif i_list[3]=='complex_substitution':
            if re.search('stop',i_list[4]):
                print("\t".join(i_list[0:3]))
cdr.close()
