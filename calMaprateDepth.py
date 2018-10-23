#!/usr/bin/env python
import os
import sys
import re
import argparse

def parse_args():
    usage = "python mapDepth.py -d pathToDir -o outfile"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-t', dest='type', metavar='protype',default='wgs',choices=['wgs','panel'],help="The project type")
    parser.add_argument('-o', dest='output', metavar='outfile',default='./cal.xls', help="The outpath and outfile")
    parser.add_argument('-p', dest='pro', metavar='project',required=True,help="The project name")
    return parser.parse_args()

def cal_wgsmapDepth(pro,output):
    prefix='/offline/Analysis/WGS'
    filepath=prefix+'/'+pro
    if os.path.isdir(filepath):
        if os.listdir(filepath):
           sample=[]
           sample.extend(list(os.walk(filepath))[0][1])
        else:
            print("The dir is empty")
            sys.exit()
    f_out=open(output,'w')
    f_out.write('Project'+'\t'+'Sample'+'\t'+'Gb'+'\t'+'MeanCov'+'\t'+'MapRate'+'\n')
    for sam in sample:
        tmpdir=os.path.join(filepath,sam+'/out'+'/'+sam+'.'+'ready_stats')
        tarf1=os.path.join(tmpdir,'genome_results.txt')
        if os.path.isfile(tarf1):
            f_in=open(tarf1,'r')
            for line in f_in:
                if 'number of mapped reads' in line:
                    maptmp=line.split('(')[1]
                    map=maptmp[0:len(maptmp)-2]
                elif 'number of sequenced bases' in line:
                    numtmp=line.split('=')[1]
                    num=round(int(numtmp[1:len(numtmp)-3].replace(',',''))*10**-9,3)
                elif 'mean coverageData' in line:
                    cover=line.split('=')[1].strip(' ').split('.')[0]
            f_out.write(pro+'\t'+sam+'\t'+str(num)+'\t'+cover+'\t'+map+'\n')
            f_in.close()
        else:
            print('%s is not existed!' % tarf1)
    f_out.close()

def cal_panelmapDepth(pro,output):
    prefix='/offline/Analysis/Panel'
    filepath=prefix+'/'+pro+'/'+'V170001-413snp'
    if os.path.isdir(filepath):
        if os.listdir(filepath):
           sample=[]
           sample.extend(list(os.walk(filepath))[0][1])
        else:
            print("The dir is empty")
            sys.exit()
    f_out=open(output,'w')
    f_out.write('Project'+'\t'+'Sample'+'\t'+'Mb'+'\t'+'MeanCov'+'\t'+'MapRate'+'\n')
    for sam in sample:
        tmpdir=os.path.join(filepath,sam+'/'+'qualimap_'+sam)
        tarf1=os.path.join(tmpdir,'genome_results.txt')
        if os.path.isfile(tarf1):
            f_in=open(tarf1,'r')

            for line in f_in:
                if 'number of mapped reads' in line:
                    maptmp=line.split('(')[1]
                    map=maptmp[0:len(maptmp)-2]
                elif 'number of sequenced bases' in line:
                    numtmp=line.split('=')[1]
                    num=round(int(numtmp[1:len(numtmp)-3].replace(',',''))*10**-6,3)
                elif 'mean coverageData' in line:
                    cover=line.split('=')[1].strip(' ').split('.')[0].replace(',','')
            f_out.write(pro+'\t'+sam+'\t'+str(num)+'\t'+cover+'\t'+map+'\n')
            f_in.close()
        else:
            print('%s is not existed!' % tarf1)
    f_out.close()

def main():
    args = parse_args()
    if args.pro and args.output and args.type=='panel':
        pro=args.pro
        output= args.output
        cal_panelmapDepth(pro,output)
    elif args.pro and args.output:
        pro=args.pro
        output= args.output
        cal_wgsmapDepth(pro,output)
   
if __name__ == "__main__":
    main()
