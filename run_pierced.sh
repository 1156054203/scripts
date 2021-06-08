#!/bin/bash

usage="Usage:\n\tbash $0 [localhost]\n"
stfdir=/home/chenyl/software/pierced/linux

if [ $# -eq 1 ];then
   port=$1
   $stfdir/ding -config=$stfdir/ding.cfg -subdomain=chenyulong $port:8080
   exit 1
elif [ $# -eq 0 ];then
   port=localhost
   $stfdir/ding -config=$stfdir/ding.cfg -subdomain=chenyulong $port:8080
else
   echo -e $usage
   exit 1
fi
