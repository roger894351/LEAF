#!/opt/anacnoda/python

## This script read in the CASM score in cdr coordinate
#usage: python casm_json_reader.py casm_file f2(.json)
import sys
import re
import time
import json

#f1=sys.argv[1]
#f1="casm_test.txt"
#with open(sys.argv[1],"r") as f:
def casm_reader(f1):
    with open(f1,"r") as f:
        a={}# this is a dictionary for variant_id 
        for i in range(1,23):
            a[str(i)]={}
        for i in "X","Y","M":
            a[str(i)]={}
        for line in f:
            if re.match("#",line):
                pass
            elif len(line.split())<=4:
                pass
            else:
                line.strip()
                line_list=line.split()
                key, value=line_list[0][3:]+"_"+"_".join(line_list[1:4]),line_list[4]#create a more accurate dictionary
                try: #if a casm was here use the larger score
                    k=a[str(line_list[0][3:])][key]
                    if float(k)<float(value):
                        a[str(line_list[0][3:])][key]=value
                except KeyError:
                    a[str(line_list[0][3:])][key]=value
    return(a)
#    print(a)
#f2=sys.argv[2]+".json"
    #           print(a)
#with open(f2, 'w') as fp:
#    json.dump(a,fp)
#print(d)
#end = time.time()
#print(end - start)
#print(casm_reader(sys.argv[1]))

#f1.close()
#fp.close()

