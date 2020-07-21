#!/usr/bin/env python3
import gzip
import itertools
import sys
from io import BufferedReader, TextIOWrapper


def read_fastq(fastq):
    """
    Return fastq records
    """
    if fastq.endswith('gz'):
        fastq_file = TextIOWrapper(BufferedReader(gzip.open(fastq, mode='rb')))
    else:
        fastq_file = open(fastq, mode='rt')

    while True:
        element = ''.join(itertools.islice(fastq_file, 4))
        if element is not '':
            yield element
        else:
            break
    fastq_file.close()

    return element

reads = read_fastq("testdir/fastq/1.rawfastq/SRR11600559_R1.fastq.gz")
for read in reads:
    sys.stdout.write(read)
