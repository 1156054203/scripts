#!/bin/bash


usage="Usage: $0 run_dir1 run_dir2 run_dir3 [run_dir4]"

list=()
if [ $# -eq 3 ];then
   run_dir1=$1
   run_dir2=$2
   run_dir3=$3
   for value1 in ${list[@]};do
      for sample1 in `ls $run_dir1|grep $value1`;do
         if [ -d $run_dir2/$sample1 ]; then
            if [ -d $run_dir3/$sample1 ]; then
               echo $sample1
               sample2_R1="$run_dir2/$sample1/${sample1#Sample_}_combined_R1.fastq.gz"
               sample2_R2="$run_dir2/$sample1/${sample1#Sample_}_combined_R2.fastq.gz"
               sample3_R1="$run_dir3/$sample1/${sample1#Sample_}_combined_R1.fastq.gz"
               sample3_R2="$run_dir3/$sample1/${sample1#Sample_}_combined_R2.fastq.gz"
                  cp -r $run_dir1/$sample1/ .
                  cat $sample2_R1 $sample3_R1 >> ./$sample1/${sample1#Sample_}_combined_R1.fastq.gz
                  cat $sample2_R2 $sample3_R2 >> ./$sample1/${sample1#Sample_}_combined_R2.fastq.gz
            else
               echo "The $run_dir3/$sample1 is not nonexistent "
            fi
         else 
            echo "The $run_dir2/$sample1 is not nonexistent "
         fi      
      done
   done
elif [ $# -eq 4 ];then
   run_dir1=$1
   run_dir2=$2
   run_dir3=$3
   run_dir4=$4
   for value1 in ${list[@]};do
      for sample1 in `ls $run_dir1|grep $value1`;do
         if [ -d $run_dir2/$sample1 ]; then
            if [ -d $run_dir3/$sample1 ]; then
               if [ -d $run_dir4/$sample1 ]; then
                  echo $sample1
                  sample2_R1="$run_dir2/$sample1/${sample1#Sample_}_combined_R1.fastq.gz"
                  sample2_R2="$run_dir2/$sample1/${sample1#Sample_}_combined_R2.fastq.gz"
                  sample3_R1="$run_dir3/$sample1/${sample1#Sample_}_combined_R1.fastq.gz"
                  sample3_R2="$run_dir3/$sample1/${sample1#Sample_}_combined_R2.fastq.gz"
                  sample4_R1="$run_dir4/$sample1/${sample1#Sample_}_combined_R1.fastq.gz"
                  sample4_R2="$run_dir4/$sample1/${sample1#Sample_}_combined_R2.fastq.gz"
                     cp -r $run_dir1/$sample1/ .
                     cat $sample2_R1 $sample3_R1 $sample4_R1 >> ./$sample1/${sample1#Sample_}_combined_R1.fastq.gz
                     cat $sample2_R2 $sample3_R2 $sample4_R2 >> ./$sample1/${sample1#Sample_}_combined_R2.fastq.gz
               else
                  echo "The $run_dir4/$sample1 is not nonexistent "
               fi
            else
               echo "The $run_dir3/$sample1 is not nonexistent "
            fi
         else
            echo "The $run_dir2/$sample1 is not nonexistent "
         fi
      done
   done
else
   echo $usage
fi
