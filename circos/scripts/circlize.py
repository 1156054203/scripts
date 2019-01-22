#!/usr/bin/env python2
"""This is a python script which:
1. reads inputs & outputs of circlize;
2. generates bed files and circos graph.
"""

import os
import sys
import glob
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--copy", action="store_true",
                    help="copy variant files to desired directory, optional")
parser.add_argument("-3", "--hg38", action="store_true",
                    help="use hg38 as assembly, optional")
parser.add_argument("-o", "--outdir", action="store", dest="outdir",
                    help="Output directory, optional")
parser.add_argument("-e", "--indel", action="store", dest="indelfile",
                    help="Annotated INDEL tsv file, optional")
parser.add_argument("-n", "--sampleacc", action="store", dest="sampleacc",
                    help="Sample Acc.")
parser.add_argument("-s", "--snpfile", action="store", dest="snpfile",
                    help="Annotated SNP tsv file")
parser.add_argument("-c", "--cnvfile", action="store", dest="cnvfile",
                    help="Annotated CNV tsv file")
parser.add_argument("-v", "--svfile", action="store", dest="svfile",
                    help="Annotated SV tsv file")
args = parser.parse_args()

if args.snpfile is None or args.cnvfile is None or \
        args.svfile is None or args.sampleacc is None:
    parser.print_help()
    sys.exit()

if args.outdir is None:
    outdir = '/offline/data-sharing/vcf/%s/healthcare' % args.sampleacc
else:
    outdir = args.outdir

sampleacc = args.sampleacc
snpfile = args.snpfile
cnvfile = args.cnvfile
svfile = args.svfile
outdir = outdir[:-1] if outdir[-1]=='/' else outdir
python2 = '/online/software/python-2.7.3/bin/python2'
r_script = '/online/software/R-3.1.2/bin/Rscript'
scripts_dir = '/online/scripts'
circos_scripts_dir = '%s/Circos' % scripts_dir
BASE_DIR = '/offline/data-sharing'
circlize_dir = '%s/%s' % (BASE_DIR, 'circlize_resources')

if not os.path.isdir(outdir): os.system('mkdir -p %s' % outdir)

if args.copy:
    annodir = '%s/vcf/%s/anno' % (BASE_DIR, sampleacc)
    if not os.path.isdir(annodir): os.system('mkdir -p %s' % annodir)
    os.system('cp %s %s/%s.snp.tsv' % (snpfile, annodir, sampleacc))
    if args.indelfile is not None:
        os.system('cp %s %s/%s.indel.tsv' % (
            args.indelfile, annodir, sampleacc))
    os.system('cp %s %s/%s.cnv.tsv' % (cnvfile, annodir, sampleacc))
    os.system('cp %s %s/%s.sv.tsv' % (svfile, annodir, sampleacc))

# Generate bed files for SNP
os.system(
    '%s %s/get_TvTiBED_nonsyn_snp_csv.py '
    '%s %s/%s_tv.bed %s/%s_ti.bed' % (
        python2, circos_scripts_dir, snpfile,
        circlize_dir, sampleacc, circlize_dir, sampleacc
    )
)

# Generate bed files for CNV
os.system(
    '%s %s/get_DupDelBED_cnv.py '
    '%s %s/%s_cnv_dup.bed %s/%s_cnv_del.bed' % (
        python2, circos_scripts_dir, cnvfile,
        circlize_dir, sampleacc, circlize_dir, sampleacc)
)
# Generate bed files for SV
os.system('%s %s/get_BED_sv.py %s %s/%s_sv.bed' % (
    python2, circos_scripts_dir, svfile, circlize_dir, sampleacc))

# Create a graph of whole genome profile
if args.hg38:
    os.system('%s %s/draw_circlize_by_sample.R %s hg38' % (
        r_script, circos_scripts_dir, sampleacc))
else:
    os.system('%s %s/draw_circlize_by_sample.R %s' % (
        r_script, circos_scripts_dir, sampleacc))

# Copy created graph to desired directory
os.system('cp %s/%s.circos.eps %s/circos-%s.eps' %(
    circlize_dir, sampleacc, outdir, sampleacc))
