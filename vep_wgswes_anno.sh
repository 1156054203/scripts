#!/bin/bash

vcftools=/online/software/vcftools-v0.1.14/bin/vcftools
cache=/online/home/wanghn/data/vep
hg38_fa=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
hg19_fa=/online/software/speedseq/genomes/ucsc_ref_hg19/ucsc.hg19.fasta
rna_fa=/online/home/chenyl/software/databases/Homo_sapiens.GRCh38.dna.primary_assembly.fa
vep_path=/online/home/chenyl/software/ensembl-vep-release-92
bed=/online/databases/Homo_sapiens/hg_exome/Agilent_V6/hg38_Agilent_V6-col6.bed


prefix=/online/home/chenyl/dongfang/wgs_wes_vcf
outdir=/online/home/chenyl/dongfang/wgs_wes_vcf/annotation
user=$1
pro=$2
tp=`echo $pro|cut -d- -f1`
ref=$(eval echo \$${pro##*-}_fa)

if [ $tp = wes ];then
   suffix=snp.vcf
else
   suffix=norm.vcf
fi

for var in `cat $prefix/$pro`;do
    sample=${var%%.*}_${tp}
    echo
    echo "${sample}...generating pbs file!"
    test -d ${outdir}/${sample}||mkdir -p ${outdir}/${sample}

    cat >${outdir}/run_${sample}_anno.pbs <<stop

#PBS -N ${sample}
#PBS -o ${outdir}/${sample}/${sample}.o
#PBS -e ${outdir}/${sample}/${sample}.e
#PBS -l nodes=1:ppn=12
#PBS -r y
#PBS -u $user

export tabix=/online/software/tabix

time $vcftools --vcf $prefix/${sample%%_*}.${suffix} --bed $bed --remove-indels --recode --recode-INFO-all --stdout | \
${vep_path}/vep --cache --dir $cache -i stdin --format vcf -o stdout --no_stats --a GRCh38 --offline --merged --fasta $ref --hgvs --sift b --polyphen b --af_1kg --af_esp --af_gnomad --pubmed --fork 4 --vcf | \
${vep_path}/filter_vep -filter "Consequence is not intergenic_variant" -o $outdir/${sample}/${sample}.anno.vcf  ##Specify input formatted as VCF,and by --vcf specify output as vcf


stop

echo "qsub run_${sample}_anno.pbs" >>${outdir}/run_anno.sh
echo "${sample}  ...Finished!"
done

sed -i "1i#!/bin/bash" ${outdir}/run_anno.sh

chmod 755 ${outdir}/run_anno.sh
