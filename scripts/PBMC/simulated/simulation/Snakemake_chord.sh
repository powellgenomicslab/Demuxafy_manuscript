#!/bin/bash

SCRIPTS=/path/to/scripts/PBMC/simulated/simulation
SNAKEFILE=$SCRIPTS/Snakefile_chord.smk
LOG=/path/to/output/PBMC/simulate/simulation/chord/logs
mkdir -p $LOG

cd $LOG



nohup snakemake --snakefile $SNAKEFILE \
    --rerun-incomplete --jobs 100 --use-conda --use-singularity --keep-going --restart-times 1 \
    --cluster "qsub -S /bin/bash {params.q} -r yes -pe smp {threads} -l h_vmem={resources.disk_per_thread_gb}G -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -e $LOG -o $LOG -j y -V" \
    > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &

