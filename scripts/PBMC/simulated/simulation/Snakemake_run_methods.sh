#!/bin/bash

SCRIPTS="/path/to/scripts/PBMC/simulated/simulation" 
SNAKEFILE=$SCRIPTS/Snakefile_run_methods.smk
YAML=$SCRIPTS/config.yaml
OUT="/path/to/output/PBMC/simulate/simulation"
LOG="$OUT/logs"
mkdir -p $LOG

cd $SCRIPTS


nohup snakemake --snakefile $SNAKEFILE \
    --rerun-incomplete --jobs 200 --use-conda --use-singularity --keep-going --restart-times 0 --configfile $YAML \
    --cluster "qsub -S /bin/bash -q short.q -r yes -pe smp {threads} -l h_vmem={resources.disk_per_thread_gb}G -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -e $LOG -o $LOG -j y -V" \
    > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &

