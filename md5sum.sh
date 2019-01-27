#!/bin/bash

if [ $# -eq 1 ];then
#
   run_dir=$1 #
fi


for folder in $run_dir;do
        echo -n "$folder" "["`date +"%Y-%m-%d %H:%M:%S"`"]" "start "
        cd $folder
        cur_dir=`pwd`
        for file in `find . -name '*.fastq.gz'`;do
                dir_name=`dirname $file`
                filename=`basename $file`
                cd ${dir_name}
                md5sum ${filename} >${filename}.md5
                cd ${cur_dir}
        done
        echo -e "["`date +"%Y-%m-%d %H:%M:%S"`"]" "done"
        cd ..
done

echo "["`date +"%Y-%m-%d %H:%M:%S"`"]" "All done"
