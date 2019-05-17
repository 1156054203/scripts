#!/usr/bin/python
# -*- coding: UTF-8 -*-
#get gene-path relations form kegg table
#chenyl:20190515
import os
import sys


def get_path(x):
    data={}
    with open(x, "r") as f:
        for row in f:
            row=row.strip('\n')
            if row.startswith('C') and '[' in row:
                path=row.split('[')[1].split(']')[0]
                data[path]=[]
            elif row.startswith('D') and 'AT' in row:
                 gene=row.split('      ')[1].split(' ')[0]
                 data[path].append(gene)
            else:
                continue
    return data


if __name__ == "__main__":
    data=get_path(sys.argv[1])
    for key,value in data.items():
        for i in value:
             sys.stdout.write(i+"\t"+key+"\n")
