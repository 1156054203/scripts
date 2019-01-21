#!/bin/bash

usage="bash $0 sheet"
if [ $# -eq 1 ];then
   sheet=$1
else 
   echo $usage
   exit 1
fi

gn=(`cat $sheet|sed -n '1~2p'|cut -f4|cut -c1`)
tumor=(`cat $sheet|sed -n '1~2p'|cut -f1`)
normal=(`cat $sheet|sed -n '2~2p'|cut -f1`)
prefix=/online/home/chenyl/dongfang
bcftools=/online/home/chenyl/software/bcftools-1.9/bcftools
bgzip=/online/home/chenyl/software/htslib-develop/bgzip
tabix=/online/home/chenyl/software/htslib-develop/tabix 



    

for i in `seq 0 $((${#tumor[@]}-1))`;do
    t=${tumor[$i]}
    n=${normal[$i]}
    pair=${t}-${n}
    group=group${gn[$i]}
    test -d $prefix/intersection/$group||mkdir $prefix/intersection/$group
    outdir=$prefix/intersection/$group
    tmp=`find $prefix/gatk/group* -type d -name "$pair"|cut -d/ -f7`
    $bcftools view $prefix/gatk/$tmp/$pair/${pair}.pass.snv.vcf -o $prefix/gatk/$tmp/$pair/${pair}.pass.snv.vcf.gz -O z
    $bcftools index $prefix/gatk/$tmp/$pair/${pair}.pass.snv.vcf.gz
    wait
    $bcftools view $prefix/strelka/$tmp/$pair/results/${pair}.pass.snv.vcf -o $prefix/strelka/$tmp/$pair/results/${pair}.pass.snv.vcf.gz -O z
    $bcftools index $prefix/strelka/$tmp/$pair/results/${pair}.pass.snv.vcf.gz
    wait
    $bcftools view $prefix/muse/$tmp/$pair/${pair}.filter.vcf -o $prefix/muse/$tmp/$pair/${pair}.filter.vcf.gz -O z
    $bcftools index $prefix/muse/$tmp/$pair/${pair}.filter.vcf.gz
    wait
    $bcftools isec -n +2 $prefix/gatk/$tmp/$pair/${pair}.pass.snv.vcf.gz $prefix/strelka/$tmp/$pair/results/${pair}.pass.snv.vcf.gz $prefix/muse/$tmp/$pair/${pair}.filter.vcf.gz -o $outdir/${pair}.isec.vcf -O v
done
