#!/bin/env python
## Stat nt freq and GC content of orf_trans_all_R64-2-1_20150113.fasta
## chenyl
from __future__ import division
import sys,re
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as mat

class baseCount:
    def __init__(self,fasta):
        self.fasta=fasta
 
    def get_stat(self):
        res={}
        with open(self.fasta) as file:
            for line in file:
                if line.startswith(">"):
                    flag=False
                    #name = line[1:].rstrip().split(' ')[0]
                    chrname=re.search('(?<=,\s)[\S\s]+?(?=\sfrom)',line).group(0)
                    if chrname not in res.keys():res[chrname]=[]
                    if 'reverse complement' in line:flag=True
                    statA,statC,statG,statT,statGC,statSum=0,0,0,0,0,0
                    continue
                if flag:
                    line=self.re_complement(line)
                    statA=line.count('A')
                    statC=line.count('C')
                    statG=line.count('G')
                    statT=line.count('T')
                    res[chrname].append([statA,statC,statG,statT])
                    
                else:
                    statA=line.count('A')
                    statC=line.count('C')
                    statG=line.count('G')
                    statT=line.count('T')
                    res[chrname].append([statA,statC,statG,statT])
        return res

    def re_complement(self,seq):
        seq = seq.upper()
        seq = seq.replace('A', 't')
        seq = seq.replace('T', 'a')
        seq = seq.replace('C', 'g')
        seq = seq.replace('G', 'c')
        return seq.upper()[::-1]
    
    def get_res(self):
        res=self.get_stat()
        gclist,ntlist,outnt=[],[],[]
        for key,value in res.items():
            nta=sum([x[0] for x in value])
            ntc=sum([x[1] for x in value])
            ntg=sum([x[2] for x in value])
            ntt=sum([x[3] for x in value])
            base=nta+ntc+ntg+ntt
            gcp=round(((ntc+ntg)/base),3)
            nap,ncp,ngp,ntp=round((nta/base),3),round((ntc/base),3),round((ntg/base),3),round((ntt/base),3)
            gclist.append([key,gcp])
            ntlist.append([key,'A',nap])
            ntlist.append([key,'C',ncp])
            ntlist.append([key,'G',ngp])
            ntlist.append([key,'T',ntp])
            outnt.append([key,nap,ncp,ngp,ntp,gcp])
        return gclist,ntlist,outnt

    def plot(self):
        gclist,ntlist,outnt=self.get_res()
        outtable=pd.DataFrame(outnt,columns=['Chr','A','C','G','T','GC'])
        gctable=pd.DataFrame(gclist,columns=['Chr','GC'])
        nttable=pd.DataFrame(ntlist,columns=['Chr','base','value'])
        #mat.rc('font',family='SimHei',size="15")
        mat.figure(figsize=(8,10))
        mat.subplot(211)
        ntfig=sn.boxplot(x='base',y='value',data=nttable)
        ntfig.set_xlabel('Nucleotide');ntfig.set_ylabel('Freq')
        mat.subplot(212)
        chr=list(gctable['Chr']);chr.sort();chrorder=chr[1:];chrorder.append(chr[0])
        gcfig=sn.barplot(x='Chr',y='GC',data=gctable,order=chrorder)
        gcfig.set_xlabel('Chromosome');gcfig.set_ylabel('GC content');gcfig.set_xticklabels(gcfig.get_xticklabels(),rotation=30)
        mat.savefig('nt_gc_freq.png') 
        return outtable

if __name__=='__main__':
    if len(sys.argv) < 2:
        print('Usage:\n     python  %s  orf_trans_all_R64-2-1_20150113.fasta' % sys.argv[0])
        sys.exit()
    instance=baseCount(sys.argv[1])
    out=instance.plot()
    out.to_csv('nt_gc_freq.txt',header=True,index=False,sep='\t',na_rep='NA')
    
