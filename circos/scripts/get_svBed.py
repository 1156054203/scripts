#!/usr/bin/env python
"""This is a python script which:
1. reads annotated *_sv.tsv from wgs;
2. output bed file of sv;
3. Pvalue less than 1e-5;
4. Duplication & Deletion not included.

Command like:
    python2.7 get_svBed.py raw.sv.xls sample.sv.bed
"""

import os
import sys
import glob
import csv

if(len(sys.argv) != 3):
    sys.exit(__doc__)

allowed_chr = ['chr%s' % x for x in [x+1 for x in range(22)]+['X', 'Y']]
sv_file, outfile= sys.argv[1:]

with open(sv_file, 'r') as infiler, \
        open(outfile, 'w') as outfiler:
    outfiler.write('chr1\tstart1\tend1\tchr2\tstart2\tend2\tType\tPvalue\n')
    for row in csv.reader(infiler, delimiter='\t'):
        if (row[0]=='chrom1' and row[5]=='chrom2'):
            continue
        chr1, start1, end1 = row[0:3]
        chr2, start2, end2 = row[5:8]
        pvalue = row[10]
        sv_type = row[13]
        if (chr1 not in allowed_chr) or \
           (chr2 not in allowed_chr) or \
           (chr1 == chr2) or \
           float(pvalue) > 1e-5:
            continue
        if ('DUPLICATION' in sv_type) or ('DELETION' in sv_type):
            continue
        else:
            outfiler.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                chr1, str(int(start1)+1), end1,
                chr2, str(int(start2)+1), end2, sv_type, pvalue))
