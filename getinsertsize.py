#!/bin/env python
#stat insert size and Standard Deviation from *fragmentLength.txt
#usage: python fragmentLength.txt [1000] 
#out: fragmentLength.rm1000.txt
#fragmentLength.txt: two columns,frist is num, second is length
#[1000]: how long segments to keep

from __future__ import division
import sys,os,math,random



def getinsertsize(infile,flag):
   file=open(infile,'r')
   name=os.path.split(infile)[1].split('.')[0]
   out=open(name+'.rm'+base+'.insertsize.txt','w')
   prop=random.uniform(0.08,0.12)
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
   propstd=average*prop
   cusum = 0
   for value in readlist:
      cusum += (value-average)**2

   std = round(math.sqrt(cusum/len(readlist)))
   out.write(name+'\t'+str(round(average))+'\t'+str(round(std))+'\t'+str(round(propstd))+'\n')
   return round(average),round(std),round(propstd)

if __name__=="__main__":
   if len(sys.argv) == 3:
      flag='true'
      max=sys.argv[2]
      base=sys.argv[2]
   elif len(sys.argv) == 2:
      flag='false'
      base='0'
   else:
      sys.exit()
   getinsertsize(sys.argv[1],flag)
