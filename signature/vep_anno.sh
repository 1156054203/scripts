#!/bin/bash

Usage="bash $0 chenyl"

if [ $# -eq 1 ];then
   user=$1
else
   echo $Usage
   exit 1
fi

vcftools=/online/software/vcftools-v0.1.14/bin/vcftools
vep_dir=/online/home/chenyl/software/ensembl-vep-release-92
agilent_v6=/online/databases/Homo_sapiens/hg_exome/Agilent_V6/hg38_Agilent_V6-col6.bed
cache=/online/home/wanghn/data/vep
hg38_fa=/online/databases/Homo_sapiens/hg38/hg38bundle/Homo_sapiens_assembly38.fasta
hg19_fa=/online/software/speedseq/genomes/ucsc_ref_hg19/ucsc.hg19.fasta
rna_fa=/online/home/chenyl/software/databases/Homo_sapiens.GRCh38.dna.primary_assembly.fa

prefix=/online/home/chenyl/dongfang/intersection

for var in `find ${prefix}/group* -maxdepth 1 -name "*.isec.vcf"`;do
    sample=`echo $var|awk -F'/' '{print$NF}'|cut -d'.' -f1`
    echo
    echo "${sample}...generating pbs file!"
    outdir=$(dirname $var)
   
    cat >${prefix}/run_${sample}_anno.pbs <<stop

#PBS -N ${sample}
#PBS -o ${outdir}/${sample}.o
#PBS -e ${outdir}/${sample}.e
#PBS -l nodes=1:ppn=12
#PBS -r y
#PBS -u $user

export tabix=/online/software/tabix

cat $var | awk '{OFS="\t"}{print\$1,\$2,".",\$3,\$4,".",\$5}'| $vep_dir/vep -i stdin --cache --dir $cache --format vcf --a GRCh38 --offline --merged --fasta ${hg38_fa} --sift b --polyphen b --af_1kg --af_esp --af_gnomad  --fork 4 --vcf --no_stats -o stdout | $vep_dir/filter_vep -filter "Consequence is not intergenic_variant" --force -o $outdir/${sample}.isec.anno.vcf
stop

chmod 755 ${prefix}/run_${sample}_anno.pbs
echo "qsub $prefix/run_${sample}_anno.pbs" >>$prefix/run_anno.sh

done

sed -i "1i #!/bin/bash" $prefix/run_anno.sh

chmod 755 $prefix/run_anno.sh
