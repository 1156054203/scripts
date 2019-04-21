#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-04-15 11:03:07
# @Author  : chenyl
# Usage    : convert chromsome length to bed 

import os
import sys

if len(sys.argv) != 2:
   print(".....\n     Usage: python sys.argv[0] infile >outfile")

   sys.exit(1)


with open(sys.argv[1], "r") as f:
    for l in f:
        chrom,length = l.rstrip().split()
        for i in range(0, int(length), 1000):
            start, stop = i, i + 1000
            if stop < int(length):
                sys.stdout.write("%s\t%s\t%s\n" % (chrom, start, stop))
            else:
                sys.stdout.write("%s\t%s\t%s\n" % (chrom, start, int(length)))
