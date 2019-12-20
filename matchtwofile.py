#!/usr/bin/env python
#usage: python script.py anno_sample.txt diff.txt out.txt
import sys 

file1=open(sys.argv[1],'r')
lines=file1.readlines()
out=open(sys.argv[3],'w')

with open(sys.argv[2],'r') as f:
    for line in f:
        flag='false'
        line=line.strip()
        for row in lines:
            peak=row.strip().split('\t')[0]
            gene=row.strip().split('\t')[8]
            if peak in line:
                out.write(line+'\t'+gene+'\n')
                flag='true'
                break
        if flag=='true':
            continue
        else:
            out.write(line+'\t'+'NA'+'\n')
file1.close()
out.close()
