#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args[1])
fh=read.table(args[1],sep="\t")
k=table(fh$V7)
#hist(fh$V7, xlim=c(0,100), breaks = 300000)
write.table(k,args[2],sep='\t')

