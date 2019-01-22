#!/usr/bin/env python
"""This is a python script which:
1. reads annotated *.cnv.tsv from wgs;
2. output bed file by dup & del.

Command like:
    python2.7 get_cnvDupDelBed.py raw.cnv.xls sample.cnv.gain.bed sample.cnv.loss.bed
"""

__author__ = "chenyulong"
__email__ = "chenyulong@gmail.com"

import os
import sys
import glob
import csv

if(len(sys.argv) != 4):
    sys.exit(__doc__)

cnv_file, dup_file, del_file = sys.argv[1:]
allowed_chr = ['chr%s' % x for x in [x+1 for x in range(22)]+['X', 'Y', 'M']]

with open(cnv_file, 'r') as infiler, \
        open(dup_file, 'w') as dupr, \
        open(del_file, 'w') as delr:
    dupr.write('chr\tstart\tend\tPvalue\n')
    delr.write('chr\tstart\tend\tPvalue\n')
    for row in csv.reader(infiler, delimiter='\t'):
        if len(row) != 13 or row[0].startswith('#'):
            continue
        corr = row[4]
        size = row[5]
        pvalue = row[7]
        if float(size) < 300000 or float(pvalue) > 1e-5:
            continue
        colon_index = corr.rfind(':')
        if colon_index == -1:
            continue
        chrname = corr[:colon_index]
        if chrname not in allowed_chr:
            continue
        pos_range = corr[colon_index+1:]
        cnv_type = row[3]
        start, end = pos_range.split('-')

        if cnv_type == 'duplication':
            dupr.write('%s\t%s\t%s\t%s\n' % (
                chrname, str(int(start)), end, pvalue))
        elif cnv_type == 'deletion':
            delr.write('%s\t%s\t%s\t%s\n' % (
                chrname, str(int(start)), end, pvalue))
        else:
            continue
