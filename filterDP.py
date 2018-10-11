#!/usr/bin/python
# -*- conding:UTF-8 -*-
import os
import sys
import argparse

__author__='chenyulong'
__doc__='filter vcf by depth'

def parse_args():
    usage = "python filterDP.py -i path/input [-o path/output] [-d depth]"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-i',dest='input',metavar='input file',required=True,help='The input file')
    parser.add_argument('-o',dest='output',metavar='output file',default='./filter.vcf',help='The output file')
    parser.add_argument('-d',dest='depth',metavar='read depth',type=int,default='4',help="The min depth. Default is 4")
    return parser.parse_args()

def filter(infile,outfile,depth):
    if os.path.isfile(infile):
         f_in=open(infile,'r')
    else:
        print('The file is not exists')
        sys.exit()
    f_out=open(outfile,'w')
    
    for line in f_in:
        line=line.strip()
        if line.startswith("#CHROM"):
            f_out.write(line+"\n")
        else:
            info=line.split("\t")[7]
            DP=info.split(';')[7]
            DPvalue=int(DP.split('=')[1])
            if DPvalue >= depth:
                f_out.write(line+'\n') 
    f_in.close()
    f_out.close()

def main():
    args=parse_args()
    if args.input and args.output and args.depth:
        f1=args.input
        f2=args.output
        dep=args.depth
        filter(f1,f2,dep)

if __name__=="__main__":
    main()
