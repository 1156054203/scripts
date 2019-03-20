#!/usr/bin/env python3

"""This script is used to get variation context, Used to replace get_spectra.R, 
   because R script run slowly,which need ten minutes per sample
"""

__author__ = "yulong chen"
__email__ = "chenyl@cloudhealth99.com"

import os
import sys
import time

start=time.clock()

# Check argv and folder existance
if not len(sys.argv) == 2: 
    sys.exit('Usage: python3 %s file' % sys.argv[0])
if not os.path.isfile(sys.argv[1]): 
    sys.exit("%s not exists." % sys.argv[1])

dir = os.path.dirname(sys.argv[1])
os.chdir(dir)
file = os.path.basename(sys.argv[1])
outfile = open(file.replace('.tmp.txt','.ready.txt'),'w')
outfile.write('chr'+'\t'+'pos'+'\t'+'ref'+'\t'+'alt'+'\t'+'flank5seq'+'\t'+'flank3seq'+'\t'+'before' \
              +'\t'+'refcontext'+'\t'+'after'+'\t'+'altcontext'+'\t'+'gene'+'\t'+'strand'+'\n')

##('chr','pos','ref','alt','flank5seq','flank3seq','before','refcontext','after','altcontext','gene','strand')
dict={'A':'T','T':'A','C':'G','G':'C'}

with open(file,'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        chr,pos,ref,alt,flank5seq,flank3seq,gene,strand=line.strip().split('\t')
        before=flank5seq[9]
        after=flank3seq[0]
        if strand == '1':
            refcontext=ref
            altcontext=alt
        else:
            refcontext=dict[ref]
            altcontext=dict[alt]

        outfile.write(chr+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+flank5seq+'\t'+flank3seq+'\t'+before \
                              +'\t'+refcontext+'\t'+after+'\t'+altcontext+'\t'+gene+'\t'+strand+'\n')
outfile.close()
end=time.clock()
print("Time used:",start-end)
