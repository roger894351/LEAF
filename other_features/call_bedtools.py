#!/usr/local/python
# pwd:/usr/local/rsrch1/epi/jchen15/database/gnomAD/script
#usage:python call_bedtools.py test_head.cdr

import subprocess
import os,sys

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
a=subprocess.Popen(["bedtools","intersect", "-a" ,"/usr/local/rsrch1/epi/jchen15/database/gnomAD/script/AF_exome_chr.bed","-b", "temp_outrandom.bed"],stdout = subprocess.PIPE)
a_data,a_error=a.communicate()
#print(a.communicate()[0])
dictionary.write(a_data)
dictionary.close()

#del_temp="rm temp_outrandom.bed temp_dictionary.dic"
#subprocess.Popen(del_temp.split())


