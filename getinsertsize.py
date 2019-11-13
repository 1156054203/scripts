#!/bin/env python
# stat insert size and Standard Deviation from *.flagstat.txt
# *.flagstat.txt file -- first column:number second column:read length
from __future__ import division
import sys,os,math

def getinsertsize(infile,flag):
   file=open(infile,'r')
   name=os.path.split(infile)[1].split('.')[0]
   out=open(name+'insertsize.txt','w')
   readlist=[]
   for line in file:
      line=line.strip()
      num=int(line.split()[0])
      length=int(line.split()[1])
      if flag=='true':
         if length <= int(max) and length > 0:
            readlist.extend([length]*num)
            continue
      else:
         if length > 0:
            readlist.extend([length]*num)
            continue

   average = sum(readlist)/len(readlist)
   cusum = 0
   for value in readlist:
      cusum += (value-average)**2

   std = math.sqrt(cusum/len(readlist))
   out.write(str(average)+'\t'+str(std))
   return average,std

if __name__=="__main__":
   if len(sys.argv) == 3:
      flag='true'
      max=sys.argv[2]
   elif len(sys.argv) == 2:
      flag='false'
   else:
      sys.exit()
   getinsertsize(sys.argv[1],flag)
