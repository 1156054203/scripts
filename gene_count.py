#! /usr/bin/env python
# _*_ coding: utf-8 _*_
#python gene_count.py /online/home/chenyl/dongfang/intersection/group1/
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import re

if __name__ == '__main__':
    r = ArgumentParser()
    r.add_argument(action="store",dest="dir",help="the group directory")
    args = r.parse_args()
    os.chdir("%s" % args.dir)
    group = (args.dir.strip().split("/"))[-2]
    result = open("%s_gene_count_dongfang.txt" % group,"w")
    result.write("gene\tbase_num\tsample_num\n")
    dict = {}
    sample_name = ""
    name = ""
    base = 0
    sample = 0
    for file in os.popen('find ./ -name "*.ready.txt"'):
        #print file
        sample_name = (re.split('[/.]', file.strip()))[2]
        print sample_name
        file = open("%s" % file.strip(),"r")
        for line in file:
            if "pos" not in line:
                gene = (line.strip().split("\t"))[-2]
                if gene not in dict:
                    name = sample_name
                    base = 1
                    sample = 1
                    dict[gene] = [base,sample]
                elif gene in dict and name != sample_name:
                    name = sample_name
                    base += 1
                    sample += 1
                    #sample = int((dict[gene])[-1]) + 1
                    dict[gene] = [base,sample]
                else:
                    base += 1
                    #sample = (dict[gene])[-1]
                    dict[gene] = [base, sample]
        file.close()
    for g in dict.keys():
        result.write("%s\t" % g)
        for num in dict[g]:
            result.write("%s\t" % num)
        result.write("\n")
    result.close()
