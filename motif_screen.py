#!/bin/env python
##chenyl
##screen motif seq of DNA_AR.fa

#from Bio import motifs
from Bio import SeqIO
import re,sys

class motif_screen:
    def __init__(self,fasta):
        #self.motif=motif
        self.fasta=fasta
 
    def getseq(self):
        record=SeqIO.read(self.fasta,'fasta')
        ntseq=record.seq
        aaseq=record.seq.translate()
        pattern = re.compile(r'[AVP].[ST]{3}')
        allmatch = re.findall(pattern,str(aaseq))
        mitseq=[]
        for i in allmatch:
            start=re.search(i,str(aaseq)).span()[0]
            end=re.search(i,str(aaseq)).span()[1]
            mitseq.append(('start: %s' % start,'end: %s' % end,'motif: %s' % str(i)))
        mitseq.append(('protein: %s' % str(aaseq)))
        return mitseq

if __name__=='__main__':
    if len(sys.argv) < 2:
        print('Usage:\n     python  %s  DNA_AR.fa' % sys.argv[0])
        sys.exit()
    motif=motif_screen(sys.argv[1])            
    print(motif.getseq())
