#!/usr/bin/env python
#encoding:utf-8
# @Author  : chenyl

import sys

list=[]
fr=sys.argv[1]
fw=sys.argv[2]

def multi2one(x,y):
    fr=open(x, 'r')
    fw=open(y, 'w')
    seq={}
    for line in fr:
        if line.startswith('>'):    #判断字符串是否以‘>开始’
            name=line.split()[0]    #以空格为分隔符，并取序列为0的项。
            seq[name]=''
        else:
            seq[name]+=line.replace('\n','')

    for i in seq.keys():
        if i not in list:
            fw.write(i+'\n')
            fw.write(seq[i]+'\n')
        else:
            continue
    fr.close()
    fw.close()


multi2one(fr,fw)
