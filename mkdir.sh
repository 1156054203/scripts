#!/bin/bash

echo The mkdir.sh script started...

cat list|awk '{print $2}'|while read a;do mkdir -p $a;done

awk '{print "mv *"$1"* "$2}' list|sh

awk '{print " echo "$2"/Sample_"$1"*"}' list|sh|sed 's/\//\\/g'|awk '{print "D:\\Data\\clinical project\\data\\"$1}'

echo The ended...

#cat list |awk '{print "D:\\Data\\clinical project\\data\\WES\\"$2"\\Sample_"$1"LU01"}'|tr "/" "\\"
#awk '{print " echo "$2"/Sample_"$1"*"}' list |sh|sed 's/\//\\/g'|awk '{print "D:\\Data\\clinical project\\data\\WES\\"$1}'|grep -v PRD|sed 's/WES/Clearseq\\HL/g'
#awk '{print " echo "$2"/Sample_"$1"*"}' list |sh|sed 's/\//\\/g'|awk '{print "D:\\Data\\clinical project\\data\\WES\\"$1}'|grep  PRD
