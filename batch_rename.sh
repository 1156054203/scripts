#!/usr/bin/bash

DIR=`pwd`

## parse md5 file
if [[ -f md5sum.txt ]]
then
	cat md5sum.txt | while read line
	do
		data=$(echo ${line} | cut -d ' ' -f 2)
		id=$(echo ${data} | cut -d '_' -f 1)
		if [[ -f ${data} ]]
		then
			echo "start move ${data} to Sample_${id}"
			if [[ ! -d Sample_${id} ]]
			then
				mkdir Sample_${id}
			fi
			mv ${data} Sample_${id}
			echo "${line}" >> Sample_${id}/${data}.md5
			echo "finished move ${data} to Sample_${id}"
		fi
	done
else
	echo "md5sum.txt done exist!!!"
fi

## md5 check
ls | grep 'Sample' | sed 's/Sample_//g' | while read line
do
	cd ${DIR}/Sample_${line} && md5sum -c *md5
done
