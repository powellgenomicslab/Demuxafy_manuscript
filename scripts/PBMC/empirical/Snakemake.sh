#!/bin/bash

SCRIPTS="/path/to/scripts/PBMC" 
SNAKEFILE="$SCRIPTS/Snakefile"
OUT="/path/to/output/PBMC"
LOG="$OUT/logs"
mkdir -p $LOG

cd $SCRIPTS



nohup snakemake --snakefile $SNAKEFILE \
    --rerun-incomplete \
    --jobs 150 \
    --use-singularity \
    --keep-going \
    --cluster "qsub -S /bin/bash -q short.q -r yes -pe smp {threads} -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -e $LOG -o $LOG -j y -V " > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &
