#!/usr/bin/env python
import os
import sys
import glob
import argparse
import subprocess

def parse_args():
    usage = "python mapDepth.py -d pathToDir -o outfile"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-t', dest='type', metavar='protype',default='wgs',choices=['wgs','wes','panel'],help="The project type")
    parser.add_argument('-o', dest='output', metavar='outfile',default='./cal.xls', help="The outpath and outfile")
    parser.add_argument('-p', dest='pro', metavar='project',required=True,help="The project name")
    return parser.parse_args()

def cal_wgsmapDepth(pro,output):
    type='WGS'
    prefix='/online/Analysis/WGS'
    filepath=prefix+'/'+pro
    if os.path.isdir(filepath):
        if os.listdir(filepath):
           sample=[]
           sample.extend(list(os.walk(filepath))[0][1])
        else:
            print("The dir is empty")
            sys.exit()
    f_out=open(output,'w')
    f_out.write('Project'+'\t'+'Sample'+'\t'+'RawData(Gb)'+'\t'+'CleanData(Gb)'+'\t'+'MapBase(Gb)'+'\t'+'MeanCov'+'\t'+'MapRate'+'\t'+'cov(Y/X)'+'\t'+'Tpye'+'\t'+'Path'+'\t'+'Items'+'\n')
    for sam in sample:
        dir1=os.path.join(filepath,sam+'/out'+'/'+sam+'.'+'ready_stats')
        f1=os.path.join(dir1,'genome_results.txt')
        if os.path.isfile(f1):
            f_in=open(f1,'r')
            for line in f_in:
                if 'number of mapped reads' in line:
                    maptmp=line.split('(')[1]
                    map=maptmp[0:len(maptmp)-2]
                elif 'number of mapped bases' in line:
                    numtmp=line.split('=')[1]
                    num=round(int(numtmp[1:len(numtmp)-3].replace(',',''))*10**-9,3)
                elif 'mean coverageData' in line:
                    cover=line.split('=')[1].strip(' ').split('.')[0]
                elif 'chrX\t' in line:
                    Xcov=float(line.split('\t')[4])
                elif 'chrY\t' in line:
                    Ycov=float(line.split('\t')[4])
            YX=round(Ycov/Xcov,5)
            f_in.close() 
        else:
            print('%s is not existed!' % f1)
        
        dir2=os.path.join(filepath,sam+'/out')
        os.chdir(dir2)
        if glob.glob(r'./*.json'):
            rawdata=subprocess.Popen(["ls |grep json|while read var;do cat $var|sed -n '5p'|sed 's/^[ \t]*//g;s/[,\"]*//g'|awk -F: '{print$2}';done|awk '{sum+=$1}END{print sum}'"],shell=True,stdout=subprocess.PIPE)
            tmp1=rawdata.communicate()[0].decode("utf-8")
            cleandata=subprocess.Popen(["ls |grep json|while read var;do cat $var|sed -n '16p'|sed 's/^[ \t]*//g;s/[,\"]*//g'|awk -F: '{print$2}';done|awk '{sum+=$1}END{print sum}'"],shell=True,stdout=subprocess.PIPE)
            tmp2=cleandata.communicate()[0].decode("utf-8")
            raw=round(int(tmp1)/10**9,2)
            clean=round(int(tmp2)/10**9,2)
            ratio=round(int(tmp2)/int(tmp1),2) 
            f_out.write(pro+'\t'+sam+'\t'+str(raw)+'\t'+str(clean)+'\t'+str(num)+'\t'+cover+'\t'+map+'\t'+str(YX)+'\t'+type+'\t'+dir2+'\n')
        else:
            print('The fastq files is not in %s' %dir2)
    f_out.close()

def cal_wesmapDepth(pro,output):
    type='WES'
    prefix='/online/Analysis/WES'
    filepath=prefix+'/'+pro
    if os.path.isdir(filepath):
        if os.listdir(filepath):
           sample=[]
           sample.extend(list(os.walk(filepath))[0][1])
        else:
            print("The dir is empty")
            sys.exit()
    f_out=open(output,'w')
    f_out.write('Project'+'\t'+'Sample'+'\t'+'RawData(Gb)'+'\t'+'CleanData(Gb)'+'\t'+'MapBase(Gb)'+'\t'+'MeanCov'+'\t'+'MapRate'+'\t'+'cov(Y/X)'+'\t'+'Tpye'+'\t'+'Path'+'\t'+'Items'+'\n')
    for sam in sample:
        tmpdir=os.path.join(filepath,sam+'/out'+'/'+'qualimap_'+sam)
        f1=os.path.join(tmpdir,'genome_results.txt')
        if os.path.isfile(f1):
            f_in=open(f1,'r')
            for line in f_in:
                if 'number of mapped reads' in line:
                    maptmp=line.split('(')[1]
                    map=maptmp[0:len(maptmp)-2]
                elif 'number of mapped bases' in line:
                    numtmp=line.split('=')[1]
                    num=round(int(numtmp[1:len(numtmp)-3].replace(',',''))*10**-9,3)
                elif 'mean coverageData' in line:
                    cover=line.split('=')[1].strip(' ').split('.')[0]
                elif 'chrX\t' in line:
                    Xcov=float(line.split('\t')[4])
                elif 'chrY\t' in line:
                    Ycov=float(line.split('\t')[4])
            YX=round(Ycov/Xcov,5)
            f_in.close()
        else:
            print('%s is not existed!' % f1)

        dir2=os.path.join(filepath,sam+'/out')
        os.chdir(dir2)
        if glob.glob(r'./*.json'):
            rawdata=subprocess.Popen(["ls |grep json|while read var;do cat $var|sed -n '5p'|sed 's/^[ \t]*//g;s/[,\"]*//g'|awk -F: '{print$2}';done|awk '{sum+=$1}END{print sum}'"],shell=True,stdout=subprocess.PIPE)
            tmp1=rawdata.communicate()[0].decode("utf-8")
            cleandata=subprocess.Popen(["ls |grep json|while read var;do cat $var|sed -n '16p'|sed 's/^[ \t]*//g;s/[,\"]*//g'|awk -F: '{print$2}';done|awk '{sum+=$1}END{print sum}'"],shell=True,stdout=subprocess.PIPE)
            tmp2=cleandata.communicate()[0].decode("utf-8")
            raw=round(int(tmp1)/10**9,2)
            clean=round(int(tmp2)/10**9,2)
            ratio=round(int(tmp2)/int(tmp1),2)
            f_out.write(pro+'\t'+sam+'\t'+str(raw)+'\t'+str(clean)+'\t'+str(num)+'\t'+cover+'\t'+map+'\t'+str(YX)+'\t'+type+'\t'+dir2+'\n')
        else:
            print('The fastq files is not in %s' %dir2)
    f_out.close()

def cal_panelmapDepth(pro,output):
    type='Panel'
    prefix='/online/Analysis/Panel'
    filepath=prefix+'/'+pro
    if os.path.isdir(filepath):
        if os.listdir(filepath):
           sample=[]
           sample.extend(list(os.walk(filepath))[0][1])
        else:
            print("The dir is empty")
            sys.exit()
    f_out=open(output,'w')
    f_out.write('Project'+'\t'+'Sample'+'\t'+'RawData(Mb)'+'\t'+'CleanData(Mb)'+'\t'+'MapBase(Mb)'+'\t'+'MeanCov'+'\t'+'MapRate'+'\t'+'cov(Y/X)'+'\t'+'Tpye'+'\t'+'Path'+'\t'+'Items'+'\n')
    for sam in sample:
        tmpdir=os.path.join(filepath,sam+'/'+'qualimap_'+sam)
        f1=os.path.join(tmpdir,'genome_results.txt')
        if os.path.isfile(f1):
            f_in=open(f1,'r')
            for line in f_in:
                if 'number of mapped reads' in line:
                    maptmp=line.split('(')[1]
                    map=maptmp[0:len(maptmp)-2]
                elif 'number of mapped bases' in line:
                    numtmp=line.split('=')[1]
                    num=round(int(numtmp[1:len(numtmp)-3].replace(',',''))*10**-6,3)
                elif 'mean coverageData' in line:
                    cover=line.split('=')[1].strip(' ').split('.')[0].replace(',','')
                elif 'chrX\t' in line:
                    Xcov=float(line.split('\t')[4])
                elif 'chrY\t' in line:
                    Ycov=float(line.split('\t')[4])
            YX=round(Ycov/Xcov,5)
            f_in.close()
        else:
            print('%s is not existed!' % f1)

        dir2=os.path.join(filepath,sam)
        os.chdir(dir2)
        if glob.glob(r'./*.json'):
            rawdata=subprocess.Popen(["ls |grep json|while read var;do cat $var|sed -n '5p'|sed 's/^[ \t]*//g;s/[,\"]*//g'|awk -F: '{print$2}';done|awk '{sum+=$1}END{print sum}'"],shell=True,stdout=subprocess.PIPE)
            tmp1=rawdata.communicate()[0].decode("utf-8")
            cleandata=subprocess.Popen(["ls |grep json|while read var;do cat $var|sed -n '16p'|sed 's/^[ \t]*//g;s/[,\"]*//g'|awk -F: '{print$2}';done|awk '{sum+=$1}END{print sum}'"],shell=True,stdout=subprocess.PIPE)
            tmp2=cleandata.communicate()[0].decode("utf-8")
            raw=round(int(tmp1)/10**6,2)
            clean=round(int(tmp2)/10**6,2)
            ratio=round(int(tmp2)/int(tmp1),2)
            f_out.write(pro+'\t'+sam+'\t'+str(raw)+'\t'+str(clean)+'\t'+str(num)+'\t'+cover+'\t'+map+'\t'+str(YX)+'\t'+type+'\t'+dir2+'\n')
        else:
            print('The fastq files is not in %s' %dir2)
    f_out.close()

def main():
    args = parse_args()
    if args.pro and args.output and args.type=='panel':
        pro=args.pro
        output= args.output
        cal_panelmapDepth(pro,output)
    elif args.pro and args.output and args.type=='wes':
        pro=args.pro
        output= args.output
        cal_wesmapDepth(pro,output)
    elif args.pro and args.output:
        pro=args.pro
        output= args.output
        cal_wgsmapDepth(pro,output)
   
if __name__ == "__main__":
    main()
