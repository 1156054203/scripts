#!/bin/bash
cache=/online/home/wanghn/data/vep
hg38_fa=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
hg19_fa=/online/software/speedseq/genomes/ucsc_ref_hg19/ucsc.hg19.fasta
vep_path=/online/home/chenyl/software/ensembl-vep-release-92
wd=`pwd`

prefix=/online/home/chenyl/temp
user=$1
pro=$2

for var in `find ${prefix}/${pro} -maxdepth 1 -name "*.vcf"`;do
    sample=`echo $var|awk '{print$NF}'|cut -d'.' -f1`
    echo
    echo "${sample}...generating pbs file!"
    outdir=$wd
    test -d ${outdir}/${sample}||mkdir -p ${outdir}/${sample}

    cat >${outdir}/run_${sample}_anno.pbs <<stop

#PBS -N ${sample}
#PBS -o ${outdir}/${sample}/${sample}.o
#PBS -e ${outdir}/${sample}/${sample}.e
#PBS -l nodes=1:ppn=20
#PBS -r y
#PBS -u $user

export tabix=/online/software/tabix

time ${vep_path}/vep --cache --dir $cache -i ${var} --format vcf -o ${outdir}/${sample}/${sample}_anno.txt --a GRCh38 --offline --merged --fasta ${hg38_fa} --sift b --polyphen b --gene_phenotype --regulatory --domains --allele_number --af --pubmed --vcf_info_field ANN --fork 4 --force --vcf  ##Specify input formatted as VCF,and by --vcf specify output as vcf

time ${vep_path}/filter_vep -i ${outdir}/${sample}/${sample}_anno.txt --filter "SIFT > 0.5 and PolyPhen " -o ${outdir}/${sample}/${sample}_anno.filter.txt --force

stop

echo "qsub run_${sample}_anno.pbs" >>${outdir}/run_${pro}_anno.sh
echo "${sample}  ...Finished!"
done
sed -i "1i#!/bin/bash" ${outdir}/run_${pro}_anno.sh

chmod 755 ${outdir}/run_${pro}_anno.sh