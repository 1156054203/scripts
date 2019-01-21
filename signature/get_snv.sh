#!/bin/bash
usage="bash $0 pro"

if [ $# -eq 1 ];then
   pro=$1
else
   echo $usage
   exit 1
fi

prefix=/online/home/chenyl/dongfang/$pro
vcftools=/online/software/vcftools-v0.1.14/bin/vcftools
bed=/online/databases/Homo_sapiens/hg_exome/Agilent_V6/hg38_Agilent_V6-col6.bed

for var in `find $prefix/group*/CHG0* -name "*.filter.vcf"`;do
    group=`echo $var|cut -d/ -f7`
    sample=`echo $var|cut -d/ -f8`
    outdir=$prefix/$group/$sample
    $vcftools --vcf $var --bed $bed --remove-indels --recode --recode-INFO-all --stdout >$outdir/${sample}.pass.snv.vcf
    #$vcftools --vcf $var --bed $bed --remove-filtered-all --remove-indels --recode --recode-INFO-all --stdout >$outdir/${sample}.pass.snv.vcf
done
