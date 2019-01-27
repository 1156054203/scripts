#!/bin/bash
prefix=/offline/Analysis/Panel/V18000-P003
pro=$1

for sample in `ls -l $prefix/$pro|grep '^d'|awk '{print$NF}'`;do

    Read1=`md5sum -t $prefix/$pro/$sample/*${sample}*_1.fq.gz|cut -d' ' -f1`
    Read2=`md5sum -t $prefix/$pro/$sample/*${sample}*_2.fq.gz|cut -d' ' -f1`
    md1=`grep "$Read1" $prefix/$pro/$sample/*${sample}*.txt|cut -d' ' -f1`
    md2=`grep "$Read2" $prefix/$pro/$sample/*${sample}*.txt|cut -d' ' -f1`

    if [ $Read1 == $md1 ];then
       echo "$sample read1 is OK : $md1"
    else
       echo "$sample read1 is failed"
    fi

    if [ $Read2 = $md2 ];then
       echo -e "$sample read2 is OK : $md2 \n"
    else
       echo -e "$sample read2 is failed \n"
    fi
done

echo 'The script is end'
