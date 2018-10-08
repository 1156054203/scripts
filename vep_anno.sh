#!/bin/bash

vcftools=/online/software/vcftools-v0.1.14/bin/vcftools
cache=/online/home/wanghn/data/vep
hg38_fa=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
hg19_fa=/online/software/speedseq/genomes/ucsc_ref_hg19/ucsc.hg19.fasta
rna_fa=/online/databases/Homo_sapiens_20170517/ensemble_GRCh38.p10/genome/hisat2_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa
vep_path=/online/home/chenyl/software/ensembl-vep-release-92

prefix=/online/home/chenyl/temp
user=$1
pro=$2

for var in `find ${prefix}/${pro} -maxdepth 1 -name "*.vcf"`;do
    sample=`echo $var|awk -F'/' '{print$NF}'|cut -d'.' -f1`
    echo
    echo "${sample}...generating pbs file!"
    outdir=`pwd`
    test -d ${outdir}/${sample}||mkdir -p ${outdir}/${sample}

    cat >${outdir}/run_${sample}_anno.pbs <<stop

#PBS -N ${sample}
#PBS -o ${outdir}/${sample}/${sample}.o
#PBS -e ${outdir}/${sample}/${sample}.e
#PBS -l nodes=1:ppn=20
#PBS -r y
#PBS -u $user

export tabix=/online/software/tabix

time $vcftools --vcf ${var} --minDP 4 --recode --recode-INFO-all --out ${outdir}/${sample}/${sample}_filter

time ${vep_path}/vep --cache --dir $cache -i ${outdir}/${sample}/${sample}_filter.recode.vcf --format vcf -o ${outdir}/${sample}/${sample}_anno.vcf --a GRCh38 --offline --merged --fasta ${hg38_fa} --hgvs --sift b --polyphen b --af_1kg --af_esp --af_gnomad --pubmed --fork 4 --force --vcf  ##Specify input formatted as VCF,and by --vcf specify output as vcf

time ${vep_path}/filter_vep -i ${outdir}/${sample}/${sample}_anno.vcf --filter "SIFT < 0.05 or PolyPhen > 0.15 and EAS_AF < 0.05 or gnomAD_EAS_AF < 0.05" -o ${outdir}/${sample}/${sample}_anno.filter.vcf --force

stop

echo "qsub run_${sample}_anno.pbs" >>${outdir}/run_${pro}_anno.sh
echo "${sample}  ...Finished!"
done

sed -i "1i#!/bin/bash" ${outdir}/run_${pro}_anno.sh

chmod 755 ${outdir}/run_${pro}_anno.sh
