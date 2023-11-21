#!/usr/local/envs/py36/bin python3

### This script runs chord for comparison to Demuxafy results
### This should be run after completing Snakefile_run_methods.smk

import glob
import os 
import pandas as pd
import re
import fnmatch
import gzip

outdir = "/path/to/output/PBMC/simulation" ### the directory that all the results are written to, this should be the same path used in the Snakefile_simulate.smk
singlet_dir = "path/to/output/PBMC/Round1Overlap" ### out path used in 'PullSingletsDoublets.R' script located at 'Demuxafy_manuscript/scripts/PBMC/empirical/PullSingletsDoublets.R'
tool_scripts = "/path/to/tool_scripts" ### path to the tool_scripts directory containing tool-specific scripts available on github and zenodo


def get_expected_Ndoublet(pool):
    lines_in_file = gzip.open(outdir + "/" + pool + "/matrix_out/barcodes.tsv.gz", 'r').readlines()
    if len(lines_in_file) < 20000:
        number_of_doublets = int(((len(lines_in_file)**2)*0.008)/1000)
    else:
        number_of_doublets = int(0.2*len(lines_in_file))
    return(number_of_doublets)

def get_expected_doubletRate(pool):
    lines_in_file = gzip.open(outdir + "/" + pool + "/matrix_out/barcodes.tsv.gz", 'r').readlines()
    if len(lines_in_file) < 20000:
        doublet_rate = ((len(lines_in_file)*0.008)/1000)
    else:
        doublet_rate = 0.2
    return(doublet_rate)


all_files = []
DoubletDecon_rhops = []
scrublet_rules = []
DoubletDetection_rules = []



files = []
dirnames = [x[0] for x in os.walk(singlet_dir)]
for dirname in dirnames:
    for fname in os.listdir(dirname):
        if re.match(r"[0-9]+_[0-9]+_updated.bam$", fname):
            files.append(fname)


if not os.path.exists(outdir + "/updated_metadata_unique.tsv"):
    metadata = pd.read_csv(outdir + "/updated_metadata.tsv", sep = "\t")
    metadata_updated = metadata.drop_duplicates()
    metadata_updated['DoubletRate'] = 0

    for pool in metadata_updated.Pool:
        ## get the number of droplets to get the doublet rate
        with gzip.open(outdir + "/" + pool + "/matrix_out/barcodes.tsv.gz", 'rb') as f:
            for i, l in enumerate(f):
                pass
        metadata_updated.loc[metadata_updated['Pool'].str.match(pool), 'DoubletRate'] = metadata_updated.loc[metadata_updated['Pool'].str.match(pool), 'Ndoublets']/(i + 1)

    metadata_updated.to_csv(outdir + "/updated_metadata_unique.tsv", sep = "\t")



sample_file = outdir + "/updated_metadata_unique.tsv"
samples = pd.read_csv(sample_file, sep = "\t")




rule all:
    input:
        expand(outdir + "/{pool}/chord/chord_results_doublet.csv", pool=samples.iloc[pd.np.r_[0:3,10:13,20:23,30:33,40:43,50:53,60:63]].Pool),
        expand(outdir + "/SimulatedOverlap/{pool}/chord_metrics.tsv", pool=samples.iloc[pd.np.r_[0:3,10:13,20:23,30:33,40:43,50:53,60:63]].Pool),


rule chord:
    input:
        matrix = outdir + "/{pool}/matrix_out/matrix.mtx.gz"
    output:
        results = outdir + "/{pool}/chord/chord_results_doublet.csv",
        variables = outdir + "/{pool}/chord/chord_variables.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 256,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 256,
    threads: 2
    params:
        q = "-q short.q",
        qstat = outdir + "/benchmarks/{pool}_chord_qstat.txt",
        matrix_dir = outdir + "/{pool}/matrix_out/",
        out = outdir + "/{pool}/chord/",
        script =tool_scripts + "/chord.R",
        dbl_rate = lambda wildcards: samples.DoubletRate[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "benchmarks/{pool}_chord.benchmark"
    conda:
        simulation_scripts_dir + "/../../../conda_environments/scrunchy.yaml"
    shell:
        """
        echo {wildcards.pool} > {output.variables}
        echo {params.out} >> {output.variables}
        echo {params.matrix_dir} >> {output.variables}
        echo {params.dbl_rate} >> {output.variables}
        Rscript {params.script} {output.variables}
        qstat -j $JOB_ID > {params.qstat}
        [[ -s {output.results} ]]
        echo $?
        """


rule chord_metrics:
    input:
        common = ancient(outdir + "/SimulatedOverlap/{pool}/singlet_counts_per_barcode_single_soft.rds"),
        barcodes = ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
        meta = sample_file
    output:
        metrics = outdir + "/SimulatedOverlap/{pool}/chord_metrics.tsv",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 24,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 24,
    threads: 2
    params:
        q = "-q short.q",
        outdir = outdir + "/{pool}/chord/",
        script = tool_scripts + "../simulated/simulation/chord_metrics.R"
    conda:
        simulation_scripts_dir + "/../../../conda_environments/generalR.yaml"
    shell:
        """
        Rscript {params.script} -o {params.outdir} -p {wildcards.pool}
        """