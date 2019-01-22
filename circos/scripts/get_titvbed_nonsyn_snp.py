#!/usr/bin/env python
"""This is a python script which:
1. reads annotated snp.tsv from wgs;
2. output bed file with Tv & Ti.

Command like:
    python2.7 get_titvbed_nonsyn_snp.py raw.snv.xls sample.ti.bed sample.tv.bed
"""

__author__ = "yulong"
__email__ = "yulong@gmail.com"

import os
import sys
import glob
import csv

if(len(sys.argv) != 4):
    sys.exit(__doc__)

tsv_file, tv_file, ti_file = sys.argv[1:]
tv_sets = [
    set(['A', 'C']), set(['A', 'T']),
    set(['G', 'C']), set(['G', 'T'])
]
ti_sets = [
    set(['A', 'G']), set(['C', 'T'])
]

allowed_chr = ['chr%s' % x for x in [x+1 for x in range(22)]+['X', 'Y', 'M']]

with open(tsv_file, 'r') as infiler, \
        open(tv_file, 'w') as tvr, \
        open(ti_file, 'w') as tir:
    tvr.write('chr\tstart\tend\tref\talt\n')
    tir.write('chr\tstart\tend\tref\talt\n')
    for row in csv.reader(infiler, delimiter='\t'):
        chrname = row[0]
        start = row[1]
        end = row[2]
        ref = row[3]
        alts = row[4]
        func = row[13]
        if chrname not in allowed_chr:
            continue
        for alt in alts:
                if set([ref, alt]) in tv_sets:
                    tvr.write('%s\t%s\t%s\t%s\t%s\n' % (
                        chrname, str(int(start)), end, ref, alt))
                elif set([ref, alt]) in ti_sets:
                    tir.write('%s\t%s\t%s\t%s\t%s\n' % (
                        chrname, str(int(start)), end, ref, alt))
                else:
                    continue
       # if 'nonsynonymous' in func:
       #     for alt in alts:
       #         if set([ref, alt]) in tv_sets:
       #             tvr.write('%s\t%s\t%s\t%s\t%s\n' % (
       #                 chrname, str(int(start)), end, ref, alt))
       #         elif set([ref, alt]) in ti_sets:
       #             tir.write('%s\t%s\t%s\t%s\t%s\n' % (
       #                 chrname, str(int(start)), end, ref, alt))
       #         else:
       #             continue
