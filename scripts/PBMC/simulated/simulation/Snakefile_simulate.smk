#!/usr/local/envs/py36/bin python3

import glob
import os 
import pandas as pd
import re
import fnmatch
import gzip

datadir = "/path/to/10x/data/dir/" ### the directory that contains all the results from cellranger for the fibroblast pools
outdir = "/path/to/output/PBMC/simulation" ### the directory that all the results are written to
singlet_dir = "path/to/output/PBMC/Round1Overlap" ### out path used in 'PullSingletsDoublets.R' script located at 'Demuxafy_manuscript/scripts/PBMC/empirical/PullSingletsDoublets.R'
gtf = "/path/to/genes/genes.gtf" ### path to gtf file
all_soft_sif = "/path/to/AllSoftwares.sif" ### the path to the singularity image "AllSoftwares.sif", provided on zenodo
demuxafy_sif = "/path/to/Demuxafy.sif" ### the path to the singularity image "Demuxafy.sif", provided on zenodo
simulation_scripts_dir = "/path/to/Deuxafy_manuscript/scripts/simulated/simulation"


all_files = []



files = []
dirnames = [x[0] for x in os.walk(singlet_dir)]
for dirname in dirnames:
    for fname in os.listdir(dirname):
        if re.match(r"[0-9]+_[0-9]+_updated.bam$", fname):
            files.append(fname)


if len(files) == 936:
    stage="simulating"
    print("simulating")
    sample_file = singlet_dir + "/SimulatedPoolsbams_increasingSizes.tsv"
    samples = pd.read_csv(sample_file, sep = "\t")
    unequal_doublets = pd.read_csv(outdir + "/unequalN_doublets.tsv", sep = "\t", header=None, names = ['Pool', 'pctl', 'Ndoublets'])
    samples_hard = samples.iloc[pd.np.r_[0:3,10:13,20:23,30:33,40:43,50:53,60:63]] ### only simulated challenging pools for three of each pool
    all_files.append(expand(outdir + "/{pool}/assigned_sorted.bam", pool=samples.Pool))
    all_files.append(expand(outdir + "/{pool}/matrix_out/matrix.mtx.gz", pool=samples.Pool))
    all_files.append(expand(outdir + "/{pool_hard}_subsampled/pooled.sorted.bam", pool_hard=samples_hard.Pool))
    all_files.append(expand(outdir + "/{pool_hard}_ambient_{pctl}pctl/pooled.sorted.bam", pool_hard=samples_hard.Pool, pctl=config["ambient_pctl"]))
    all_files.append(expand(outdir + "/{pool_hard}_mt_{pctl}pctl/pooled.sorted.bam", pool_hard=samples_hard.Pool, pctl=config["mt_pctl"]))
    all_files.append(expand(outdir + "/{pool_hard}_mt_{pctl}pctl/pooled.sorted.bam", pool_hard=samples_hard.Pool, pctl=config["mt_pctl"]))
    all_files.append(expand(outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/pooled.sorted.bam", pool_hard=samples_hard.Pool, pctl_uneven=config["pctl_uneven"]))


else:
    stage="preprocessing"
    print("preprocessing")
    pool_indiv_meta = pd.read_csv(singlet_dir + "/Pool_individuals_meta.tsv", sep = "\t")
    all_files.append(expand(singlet_dir + "/{pool_original}/{indiv}_updated.bam", zip, pool_original=pool_indiv_meta.Pool, indiv = pool_indiv_meta.Individuals))



wildcard_constraints:
    pool= '|'.join([re.escape(x) for x in samples.Pool]),
    pool_hard= '|'.join([re.escape(x) for x in samples_hard.Pool]),
    pool_original= '|'.join([re.escape(x) for x in samples_hard.Pool]),
    indiv = '|'.join([re.escape(str(x)) for x in pool_indiv_meta.Individuals]),


rule all:
    input:
        all_files,



if stage == "preprocessing":

    rule split_bams_by_individual:
        input:
            bam = datadir + "/{pool_original}_V1/outs/possorted_genome_bam.bam",
            barcodes = singlet_dir + "/{pool_original}/{pool_original}_droplet_barcodes.tsv"
        output:
            done = singlet_dir + "/{pool_original}/subset_bam.done",
            qstat = outdir + "/benchmarks/{pool_original}_split_bams_by_individual.qstat"
        params:
            out = singlet_dir + "/{pool_original}/",
        threads: 12
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 8
        benchmark:
            outdir + "/benchmarks/{pool_original}_split_bams_by_individual.benchmark"
        conda:
            simulation_scripts_dir + "/../../../conda_environments/sinto.yaml"
        shell:
            """
            sinto filterbarcodes -b {input.bam} -c {input.barcodes} --barcodetag CB --outdir {params.out} --nproc {threads}
            qstat -j $JOB_ID > {output.qstat}
            echo "done with splitting bam" > {output.done}
            """

    rule update_bam_barcodes:
        input:
            done = singlet_dir + "/{pool_original}/subset_bam.done"
        output:
            sam = temp(singlet_dir + "/{pool_original}/{indiv}_updated.sam"),
            bam = singlet_dir + "/{pool_original}/{indiv}_updated.bam",
            bai = singlet_dir + "/{pool_original}/{indiv}_updated.bam.bai",
            qstat = outdir + "/benchmarks/{pool_original}_{indiv}_update_bam_barcodes.qstat"
        params:
            sif = demuxafy_sif,
            bam = singlet_dir + "/{pool_original}/{indiv}.bam"
        threads: 1
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 12,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 12
        benchmark:
            outdir + "/benchmarks/{pool_original}_{indiv}_update_bam_barcodes.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools view -@ {threads} -H {params.bam} > {output.sam}
            singularity exec {params.sif} samtools view -@ {threads} {params.bam} | sed "s/-1\tUR/-{wildcards.indiv}\tUR/g" >> {output.sam}
            singularity exec {params.sif} samtools view -@ {threads} -S -b {output.sam} > {output.bam}
            singularity exec {params.sif} samtools index -@ {threads} {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            rm {params.bam}
            """


if stage == "simulating":

    rule merge_bams4simulation:
        input:
            barcodes = outdir + "/{pool}/barcodes4simulation.tsv"
        output:
            bam = temp(outdir + "/{pool}/combined_singlets.bam"),
            bai = temp(outdir + "/{pool}/combined_singlets.bam.bai"),
            qstat = outdir + "/benchmarks/{pool}_merge_bams4simulation.qstat"
        params:
            bams=lambda wildcards: samples.Bams[samples.Pool == wildcards.pool].iloc[0].replace(",", " "),
            sif = demuxafy_sif,
        threads: 8
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 4
        benchmark:
            outdir + "/benchmarks/{pool}_merge_bams4simulation.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools merge -@ {threads} {output.bam} {params.bams}
            singularity exec {params.sif} samtools index -@ {threads} {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            """



    rule simulate_bam:
        input:
            barcodes = outdir + "/{pool}/barcodes4simulation.tsv",
            bam = outdir + "/{pool}/combined_singlets.bam",
            bai = outdir + "/{pool}/combined_singlets.bam.bai"
        output:
            bam = temp(outdir + "/{pool}/simulated_doublets.bam"),
            qstat = outdir + "/benchmarks/{pool}_simulate_bam.qstat"
        params:
            N=lambda wildcards: samples.Ndoublets[samples.Pool == wildcards.pool].iloc[0],
            sif = demuxafy_sif
        threads: 8
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 24,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 24
        benchmark:
            outdir + "/benchmarks/{pool}_simulate_bam.benchmark"
        shell:
            """
            singularity exec {params.sif} java -jar /opt/Drop-seq_tools-2.5.4/jar/dropseq.jar GenerateSyntheticDoublets \
                                                                    --NUMBER_MULTICELL {params.N} \
                                                                    --CELL_BARCODE_TAG "CB" \
                                                                    --OUTPUT {output.bam} \
                                                                    --INPUT {input.bam} \
                                                                    --CELL_BC_FILE {input.barcodes} \
                                                                    --EMIT_SINGLETONS false
            qstat -j $JOB_ID > {output.qstat}
            """


    rule combine_simulated_bams:
        input:
            bam = outdir + "/{pool}/simulated_doublets.bam"
        output:
            bam = temp(outdir + "/{pool}/pooled.bam"),
            bai = temp(outdir + "/{pool}/pooled.bam.bai"),
            qstat = outdir + "/benchmarks/{pool}_combine_simulated_bams.qstat"
        params:
            bams=lambda wildcards: samples.Bams[samples.Pool == wildcards.pool].iloc[0].replace(",", " "),
            sif = demuxafy_sif,
        threads: 8
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 4
        benchmark:
            outdir + "/benchmarks/{pool}_combine_simulated_bams.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools merge -@ {threads} {output.bam} {params.bams} {input.bam}
            singularity exec {params.sif} samtools index {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            """


    ##### Sort and index bam #####
    rule sort_index_count:
        input:
            bam = outdir + "/{pool}/pooled.bam",
            bai = outdir + "/{pool}/pooled.bam.bai"
        threads: 4
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
        output:
            bam = temp(outdir + "/{pool}/pooled.sorted.bam"),
            bai = temp(outdir + "/{pool}/pooled.sorted.bam.bai"),
            qstat = outdir + "/benchmarks/{pool}_sort_index_count.qstat"
        params:
            sif=demuxafy_sif,
            out=outdir + "/{pool}/"
        benchmark:
            outdir + "/benchmarks/{pool}_sort_index_count.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools sort -@ {threads} -T {params.out} -o {output.bam} {input.bam}
            singularity exec {params.sif} samtools index {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            """


    #########################################################
    ########## CREATE MATRICES FOR SIMULATED POOLS ##########
    #########################################################
    ##### Assign genes to bam reads #####
    rule feature_counts:
        input:
            gtf=GTF,
            bam = outdir + "/{pool}/pooled.sorted.bam",
            bai = outdir + "/{pool}/pooled.sorted.bam.bai"
        output:
            count_bam = temp(outdir + "/{pool}/pooled.sorted.bam.featureCounts.bam"),
            gene_assigned = temp(outdir + "/{pool}/gene_assigned"),
            qstat = outdir + "/benchmarks/{pool}_feature_counts.qstat",
        params:
            o=outdir + "/{pool}/gene_assigned",
            sif=demuxafy_sif,
            out = outdir + "/{pool}/"
        threads: 10
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 16,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 16
        benchmark:
            outdir + "/benchmarks/{pool}_feature_counts.benchmark"
        conda:
            simulation_scripts_dir + "/../../../conda_environments/subread.yaml"
        shell:
            """
            featureCounts -a {input.gtf} -o {params.o} -R BAM -T {threads} {input.bam}
            qstat -j $JOB_ID > {output.qstat}
            """



    ##### Sort and index bam #####
    rule sort_index:
        input:
            bam = outdir + "/{pool}/pooled.sorted.bam.featureCounts.bam",
            gene_assigned = outdir + "/{pool}/gene_assigned",
        threads: 4
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
        output:
            bam = temp(outdir + "/{pool}/featureCounts.sorted.bam"),
            bai = temp(outdir + "/{pool}/featureCounts.sorted.bam.bai"),
            qstat = outdir + "/benchmarks/{pool}_sort_index.qstat"
        params:
            sif=demuxafy_sif,
            out = outdir + "/{pool}/"
        benchmark:
            outdir + "/benchmarks/{pool}_sort_index.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools sort -@ {threads} -T {params.out} -o {output.bam} {input.bam}
            singularity exec {params.sif} samtools index {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            """


    rule deduplicate_combined_bam:
        input:
            bam = outdir + "/{pool}/featureCounts.sorted.bam",
            bai = outdir + "/{pool}/featureCounts.sorted.bam.bai",
        output:
            bam = temp(outdir + "/{pool}/assigned.bam"),
            qstat = outdir + "/benchmarks/{pool}_deduplicate_combined_bam.qstat"
        params:
            sif = demuxafy_sif,
            log = outdir + "/{pool}/deduplicate_combined_bam.log",
            out = outdir + "/{pool}/"
        threads: 8
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 8
        benchmark:
            outdir + "/benchmarks/{pool}_deduplicate_combined_bam.benchmark"
        conda:
            simulation_scripts_dir + "/../../../conda_environments/umi_tools.yaml"
        shell:
            """
            umi_tools dedup -I {input.bam} --paired --log={params.log} --umi-tag=UB --cell-tag=CB --extract-umi-method=tag --cell-tag-split=None --gene-tag=XT --assigned-status-tag=XS --per-cell --per-gene --method unique --temp-dir={params.out} > {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            """


    ##### Sort and index bam #####
    rule sort_index_assigned:
        input:
            outdir + "/{pool}/assigned.bam",
        threads: 1
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
        output:
            bam = outdir + "/{pool}/assigned_sorted.bam",
            bai = outdir + "/{pool}/assigned_sorted.bam.bai",
            qstat = outdir + "/benchmarks/{pool}_sort_index.qstat"
        params:
            sif=demuxafy_sif,
            out = outdir + "/{pool}/"
        benchmark:
            outdir + "/benchmarks/{pool}_sort_index.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools sort -@ {threads} -T {params.out} -o {output.bam} {input}
            singularity exec {params.sif} samtools index {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            """

    ##### create umi_tools count files from simulated bams for scrublet #####
    rule umi_tools:
        input:
            outdir + "/{pool}/assigned_sorted.bam"
        output:
            counts = outdir + "/{pool}/counts.tsv.gz",
            qstat = outdir + "/benchmarks/{pool}_umi_tools.qstat"
        threads: 5
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 5
        params: 
            sif=demuxafy_sif
        benchmark:
            outdir + "/benchmarks/{pool}_umi_tools.benchmark"
        conda:
            simulation_scripts_dir + "/../../../conda_environments/umi_tools.yaml"
        shell:
            """
            umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --umi-tag=UB --cell-tag=CB --extract-umi-method=tag --cell-tag-split=None -I {input} -S {output.counts}
            qstat -j $JOB_ID > {output.qstat}
            """

    ##### create matrix files from counts files simulated bams for scrublet #####
    rule umi_tools_2_mtx:
        input:
            # ancient(outdir + "/{pool}/counts.tsv.gz")
            outdir + "/{pool}/counts.tsv.gz"
        output:
            matrix = outdir + "/{pool}/matrix_out/matrix.mtx.gz",
            features = outdir + "/{pool}/matrix_out/features.tsv.gz",
            barcodes = outdir + "/{pool}/matrix_out/barcodes.tsv.gz",
            variables = outdir + "/{pool}/matrix_out/variables.txt",
            qstat = outdir + "/benchmarks/{pool}_umi_tools_2_mtx.qstat"
        params:
            pool_dir = outdir + "/{pool}/",
            script = simulation_scripts_dir +"/umi_counts_mtx.R"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 5
        threads: 4
        benchmark:
            outdir + "/benchmarks/{pool}_umi_tools_2_mtx.benchmark"
        conda:
            simulation_scripts_dir + "/../../../conda_environments/scrunchy.yaml"
        shell:
            """
            echo {params.pool_dir} > {output.variables}
            Rscript {params.script} {output.variables}
            qstat -j $JOB_ID > {output.qstat}
            """

    ######################################
    ########## HARD SIMULATIONS ##########
    ######################################
    rule new_meta_sample:
        input:
            # meta_file = ancient(sample_file)
            meta_file = sample_file
        output:
            ref = outdir + "/sample_meta.tsv",
            new_sample_file = outdir + "/updated_metadata.tsv",
        params:
            sif = demuxafy_sif
        threads: 1
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 1
        shell:
            """
            singularity exec {params.sif} cat {input.meta_file} > {output.ref}
            singularity exec {params.sif} cat {input.meta_file} > {output.new_sample_file}
            """

    rule ambient:
        input:
            matrix = outdir + "/{pool_hard}/matrix_out/matrix.mtx.gz",
            features = outdir + "/{pool_hard}/matrix_out/features.tsv.gz",
            barcodes = outdir + "/{pool_hard}/matrix_out/barcodes.tsv.gz",
            bam = outdir + "/{pool_hard}/assigned_sorted.bam",
            sample_file = outdir + "/sample_meta.tsv",
        output:
            bam = outdir + "/{pool_hard}_ambient_{pctl}pctl/pooled.sorted.bam",
            variables = outdir + "/{pool_hard}_ambient_{pctl}pctl/variables.txt",
            qstat = outdir + "/benchmarks/{pool_hard}_{pctl}pctl_ambient.qstat"
        params:
            sif = demuxafy_sif,
            script = simulation_scripts_dir + "/Ambient_simulations.R",
            outdir = outdir + "/{pool_hard}_ambient_{pctl}pctl/",
            indir = outdir + "/{pool_hard}",
            new_sample = outdir + "/updated_metadata.tsv",
            matrix_dir = outdir + "/{pool_hard}/matrix_out/",
        threads: 32
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 32,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 32
        benchmark:
            outdir + "/benchmarks/{pool_hard}_{pctl}pctl_ambient.benchmark"
        conda:
            simulation_scripts_dir + "/../../../conda_environments/simulations.yaml"
        shell:
            """
            echo {wildcards.pool_hard} > {output.variables}
            echo {wildcards.pctl} >> {output.variables}
            echo {params.indir} >> {output.variables}
            echo {params.outdir} >> {output.variables}
            echo {params.matrix_dir} >> {output.variables}
            echo {input.bam} >> {output.variables}
            Rscript {params.script} {output.variables}
            grep {wildcards.pool_hard} {input.sample_file} | sed 's/{wildcards.pool_hard}\t/{wildcards.pool_hard}_ambient_{wildcards.pctl}pctl\t/g' >> {params.new_sample}
            qstat -j $JOB_ID > {output.qstat}
            """


    rule mt_percent:
        input:
            bam = outdir + "/{pool_hard}/assigned_sorted.bam",
            sample_file = outdir + "/sample_meta.tsv",
            matrix = outdir + "/{pool_hard}/matrix_out/matrix.mtx.gz",
        output:
            bam = outdir + "/{pool_hard}_mt_{pctl}pctl/pooled.sorted.bam",
            variables = outdir + "/{pool_hard}_mt_{pctl}pctl/variables.txt",
            qstat = outdir + "/benchmarks/{pool_hard}_{pctl}pctl_mt_percent.qstat"
        params:
            script = simulation_scripts_dir + "/Mt_percent_simulations.R",
            outdir = outdir + "/{pool_hard}_mt_{pctl}pctl/",
            indir = outdir + "/{pool_hard}/",
            matrix_dir = outdir + "/{pool_hard}/matrix_out/",
            new_sample = outdir + "/updated_metadata.tsv"
        threads: 32
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 8
        benchmark:
            outdir + "/benchmarks/{pool_hard}_{pctl}pctl_mt_percent.benchmark"
        conda:
            simulation_scripts_dir + "/../../../conda_environments/simulations.yaml"
        shell:
            """
            echo {wildcards.pool_hard} > {output.variables}
            echo {wildcards.pctl} >> {output.variables}
            echo {params.indir} >> {output.variables}
            echo {params.outdir} >> {output.variables}
            echo {params.matrix_dir} >> {output.variables}
            echo {input.bam} >> {output.variables}
            echo {threads} >> {output.variables}
            Rscript {params.script} {output.variables}
            grep "{wildcards.pool_hard}\t" {input.sample_file} | sed 's/{wildcards.pool_hard}\t/{wildcards.pool_hard}_mt_{wildcards.pctl}pctl\t/g' >> {params.new_sample}
            qstat -j $JOB_ID > {output.qstat}
            """

    rule unevenN:
        input:
            sample_file = sample_file,
            bam = outdir + "/{pool_hard}/assigned_sorted.bam",
            barcodes = outdir + "/{pool_hard}/barcodes4simulation.tsv",
        output:
            barcodes = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/barcodes4simulation.tsv",
            bam = temp(outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/combined_singlets.bam"),
            bai = temp(outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/combined_singlets.bam.bai"),
            variables = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/variables.txt",
            qstat = outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl.qstat"
        params:
            script = simulation_scripts_dir + "Unequal_N_simulations.R",
            outdir = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/",
            outbase = outdir,
            new_sample = outdir + "/updated_metadata.tsv"
        threads: 4
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 256,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 256
        benchmark:
            outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl.benchmark"
        conda:
            simulation_scripts_dir + "/../../../conda_environments/simulations.yaml"
        shell:
            """
            echo {wildcards.pool_hard} > {output.variables}
            echo {wildcards.pctl_uneven} >> {output.variables}
            echo {params.outdir} >> {output.variables}
            echo {threads} >> {output.variables}
            echo {input.sample_file} >> {output.variables}
            echo {input.barcodes} >> {output.variables}
            echo {params.outbase} >> {output.variables}
            Rscript {params.script} {output.variables}
            grep -P '{wildcards.pool_hard}\t' {input.sample_file} | sed 's/{wildcards.pool_hard}\t/{wildcards.pool_hard}_unevenN_{wildcards.pctl_uneven}pctl\t/g' >> {params.new_sample}
            qstat -j $JOB_ID > {output.qstat}
            """


    rule simulate_bam_uneven:
        input:
            barcodes = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/barcodes4simulation.tsv",
            bam = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/combined_singlets.bam",
            bai = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/combined_singlets.bam.bai",
            qstat = outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl.qstat"
        output:
            bam = temp(outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/simulated_doublets.bam"),
            qstat = outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl_simulate_bam.qstat"
        params:
            N = lambda wildcards: unequal_doublets.Ndoublets[(unequal_doublets.Pool == wildcards.pool_hard) & (unequal_doublets.pctl == float(wildcards.pctl_uneven))].iloc[0],
            sif = demuxafy_sif,
        threads: 4
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 256,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 256,
            java_mem = lambda wildcards, attempt: attempt * 256
        benchmark:
            outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl_simulate_bam.benchmark"
        shell:
            """
            echo {params.N}
            export _JAVA_OPTIONS="-Xmx{resources.java_mem}g -Xms{resources.java_mem}g -XX:MaxHeapSize={resources.java_mem}g -XX:+UseCompressedOops -XX:MaxPermSize=1G"

            singularity exec {params.sif} java -jar /opt/Drop-seq_tools-2.5.4/jar/dropseq.jar GenerateSyntheticDoublets \
                                                                                                --NUMBER_MULTICELL {params.N} \
                                                                                                --CELL_BARCODE_TAG "CB" \
                                                                                                --OUTPUT {output.bam} \
                                                                                                --INPUT {input.bam} \
                                                                                                --CELL_BC_FILE {input.barcodes} \
                                                                                                --EMIT_SINGLETONS false \
                                                                                                --TMP_DIR $TMPDIR \
                                                                                                --MAX_RECORDS_IN_RAM 5000000
            qstat -j $JOB_ID > {output.qstat}
            """


    rule combine_simulated_bams_uneven:
        input:
            singlet_bam = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/combined_singlets.bam",
            bam = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/simulated_doublets.bam"
        output:
            bam = temp(outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/pooled.bam"),
            bai = temp(outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/pooled.bam.bai"),
            qstat = outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl_combine_simulated_bams.qstat"
        params:
            sif = demuxafy_sif,
            
        threads: 8
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 4
        benchmark:
            outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl_combine_simulated_bams.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools merge -@ {threads} {output.bam} {input.singlet_bam} {input.bam}
            singularity exec {params.sif} samtools index {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            """


    ##### Sort and index bam #####
    rule sort_index_count_uneven:
        input:
            bam = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/pooled.bam",
            bai = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/pooled.bam.bai"
        threads: 4
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
        output:
            bam = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/pooled.sorted.bam",
            bai = outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/pooled.sorted.bam.bai",
            qstat = outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl_sort_index_count.qstat"
        params:
            sif=demuxafy_sif,
            out=outdir + "/{pool_hard}_unevenN_{pctl_uneven}pctl/"
        benchmark:
            outdir + "/benchmarks/{pool_hard}_unevenN_{pctl_uneven}pctl_sort_index_count.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools sort -@ {threads} -T {params.out} -o {output.bam} {input.bam}
            singularity exec {params.sif} samtools index {output.bam}
            qstat -j $JOB_ID > {output.qstat}
            """


    rule subsample:
        input:
            bam = outdir + "/{pool_hard}/assigned_sorted.bam",
            sample_file = outdir + "/sample_meta.tsv",
        output:
            bam = outdir + "/{pool_hard}_subsampled/pooled.sorted.bam",
            qstat = outdir + "/benchmarks/{pool_hard}_subsample.qstat"
        params:
            sif = demuxafy_sif,
            seed = 17,
            proportion = 667,
            new_sample = outdir + "/updated_metadata.tsv"
        threads: 4
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 8
        benchmark:
            outdir + "/benchmarks/{pool_hard}_subsample.benchmark"
        shell:
            """
            singularity exec {params.sif} samtools view -bs {params.seed}.{params.proportion} {input.bam} > {output.bam}
            singularity exec {params.sif} grep {wildcards.pool_hard} {input.sample_file} | sed 's/{wildcards.pool_hard}\t/{wildcards.pool_hard}_subsampled\t/g' >> {params.new_sample}
            qstat -j $JOB_ID > {output.qstat}
            """
