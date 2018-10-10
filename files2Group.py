#!/usr/bin/python
# -*- conding:UTF-8 -*-
import os
import sys
import argparse

def parse_args():
    usage = "python files2Group.py -d pathToDir -k fileNameKeyWord"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-k', dest='key', metavar='filename key', required=True,help="The keyword of file name")
    parser.add_argument('-d', dest='path', metavar='pathto dir', default='.', help="The file directory. Default is the current directory")
    return parser.parse_args()


def merge(filePath,keyword):
    filelist=[]
    if os.path.isdir(filePath):
        if os.listdir(filePath):
           filelist.extend(os.listdir(filePath))
        else:
            print("The dir is empty") 
            sys.exit(1)
    else:
        print('The dir is not exist')
        sys.exit(1)
   
    filelist=[x for x in filelist if keyword in x]
    start=0
    end=20
    tmp=[]
    index=1
    tag=1
    while tag:
        if end < len(filelist):
            tmp=filelist[start:end]
            f_write=open(filePath+'/'+'group'+str(index)+'.stat.xls','w')
            head=open(filePath+'/'+tmp[0],'r').readlines()[0].split('\t')
            head=head[0]+'\t'+head[3]+'\t'+head[4]+'\t'+head[5]
            f_write.write(head)
            for file in tmp:
                f_read=open(filePath+'/'+file,'r')
                context=f_read.readlines()[1:]
                for i in context:
                    i=i.split('\t')
                    i=i[0]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]
                    f_write.write(i)
        else:
            tmp=filelist[start:len(filelist)]
            f_write=open(filePath+'/'+'group'+str(index)+'.stat.xls','w')
            head=open(filePath+'/'+tmp[0],'r').readlines()[0].split('\t')
            head=head[0]+'\t'+head[3]+'\t'+head[4]+'\t'+head[5]
            f_write.write(head)
            for file in tmp:
                f_read=open(filePath+'/'+file,'r')
                context=f_read.readlines()[1:]
                for i in context:
                    i=i.split('\t')
                    i=i[0]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]
                    f_write.write(i)
            tag=0
        f_write.close()
        f_read.close()
        start+=20
        end+=20
        index+=1
def main():
    args = parse_args()
    if args.key and args.path:
        path=args.path
        key= args.key
        merge(path,key)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
