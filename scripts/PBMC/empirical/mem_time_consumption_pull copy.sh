#!/bin/bash

DIR=/path/to/output/PBMC/benchmarks
OUT=$DIR/qstat_usage

mkdir -p $OUT

files=`ls $DIR/*qstat.txt`

for file in $files
do
    name=`basename $file`
    echo $name
    grep "usage" $file | sed 's/usage         1:            //g' | sed 's/ //g' > $OUT/usage_$name
done
