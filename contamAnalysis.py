#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Modified on Tue Oct 22 2018

"""
import os
from collections import OrderedDict
from argparse import ArgumentParser
import random
import time
import gzip
import shlex,subprocess

def parse_args():
    global usage
    usage='contamAnalysis.py -m 100000 -t 8 -r ref -f1 path/read1 -f2 path/read2 -b path/blastn -n path/nt -o outDir'
    parser = ArgumentParser(usage=usage)
    parser.add_argument("-f1",dest="Read1",action="store")
    parser.add_argument("-f2",dest="Read2",action="store",required=False)
    parser.add_argument("-b",dest="blast",action="store")
    parser.add_argument("-n",dest="nt",action="store")
    parser.add_argument("-o",dest="outDir",action="store")
    parser.add_argument("-r",dest="ref",action="store")
    parser.add_argument("-m",dest='num',default=100000,type=int)
    parser.add_argument("-t",dest='threads',default=8,type=int)
    return parser.parse_args()


start_time = time.time()
print("Run start: %s" % time.strftime('%Y-%m-%d %H:%M:%S',time.localtime()))

# 解压随机抽取
def extractFastq(filename,readsNum,tmpName):
    fn_open = gzip.open if filename.endswith(".gz") else open
    cwd = os.getcwd()
    faTmp = open(cwd+os.sep+str(tmpName),"w")
    
    with fn_open(filename) as fh:
        num_lines = sum([1 for line in fh])
        total_records = int(num_lines/4)
        fh.seek(0)        

        random.seed(100)
        output_set = set(random.sample(range(total_records),readsNum))
        
        record_number = 0
        for line1 in fh:
            line2 = next(fh)
            next(fh)
            next(fh)
            if record_number in output_set:
                faTmp.write(">"+line1.decode("utf-8"))
                faTmp.write(line2.decode("utf-8"))
            record_number += 1 
    faTmp.close()

# 取前xx行
def shell_head(filename,readsNum,tmpName):
    fn_open = gzip.open if filename.endswith(".gz") else open
    cwd = os.getcwd()
    faTmp = open(cwd+os.sep+str(tmpName),"w")

    with fn_open(filename) as fh:
        record_num = 0
        for line1 in fh:
            line2 = next(fh)
            next(fh)
            next(fh)
            record_num += 1
            faTmp.write(">"+line1.decode("utf-8"))
            faTmp.write(line2.decode("utf-8"))
            if record_num >= readsNum:
                break
    faTmp.close()

# 排序 
def order(outDir,ref):
    tmp_dict = OrderedDict()
    with open(args.outDir+os.sep+args.ref+".result") as fi:
         for line in fi:
             species = line.split("\t")[8]
             if not species in tmp_dict.keys():
                 tmp_dict[species] = 1
             else:
                 tmp_dict[species] += 1

    tmp_lst = sorted(tmp_dict.items(), key=lambda item: item[1],reverse = True)  
    with open(output+os.sep+args.ref+".20","w") as fii:
        i = 0
        for line in tmp_lst:
            fii.write("\t".join(map(str,line))+"\n")
            i += 1
            if i >= 20:
                break

def main():
    args = parse_args()
    if args.Read1 is not None:
        Read1 = args.Read1
    else:
        print(usage)
        exit("Error: a Read1 is required!")
    if args.Read2 is not None:
        Read2 = args.Read2
   
    if args.blast is not None:
        blast = args.blast
    else:
        print(usage)
        exit("Error: a blast is required!")
    if args.nt is not None:
        nt = args.nt
    else:
        print(usage)
        exit("Error: nt database is required!")
    if args.outDir is not None:
        outDir = args.outDir
        outDir = outDir.rstrip("/")
    else:
        print(usage)
        exit("Error: Output directory is required!")
    if args.ref is not None:
        ref = args.ref
    else:
        print(usage)
        exit("Error: a ref is required!")
    if args.num is not None:
        num = args.num
    else:
        print(usage)
        exit("Error: a num is required!")
    if args.threads is not None:
        threads = args.threads
    else:
        print(usage)
        exit("Error: a threads is required!")

    hell_head(Read1,num,"tmp_R1.fa")
    if args.Read2:
        shell_head(Read2,num,"tmp_R2.fa")
    else:
        pass

    cat_ps = subprocess.Popen(["cat","tmp_R1.fa","tmp_R2.fa"],stdout=subprocess.PIPE)
    faTmp = open("random_extract.fa","w")
    faTmp.write(cat_ps.communicate()[0].decode("utf-8"))
    faTmp.close()
    
    cmd = "{} -query random_extract.fa -db {} -task megablast  -num_threads {} \
    -out {} -evalue 0.01 -outfmt '6 qseqid qlen sseqid slen length pident evalue sscinames \
    stitle staxids sskingdoms' -max_target_seqs 1".format(blast,nt,threads,outDir+os.sep+ref+".result")
    subprocess.call(shlex.split(cmd))
    
    order(outDir,ref)

    subprocess.call("rm -f tmp_R1.fa".split(" "))
    subprocess.call("rm -f tmp_R2.fa".split(" "))
    subprocess.call("rm -f random_extract.fa".split(" "))    

if __name__ == "__main__":
    main()

run_time=int(time.time() - start_time)
print("Run ended: %s" % time.strftime('%Y-%m-%d %H:%M:%S',time.localtime()))
print()
print("Run time: %s" % time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(run_time)))
