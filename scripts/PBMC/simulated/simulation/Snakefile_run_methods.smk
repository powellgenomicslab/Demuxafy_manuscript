#!/usr/local/envs/py36/bin python3

### This script runs all the demultiplexing and doublet detecting methods as well as scripts to assess each
### This should be run after completing Snakefile_simulate.smk

import glob
import os 
import pandas as pd
import re
import fnmatch
import gzip


datadir = "/path/to/10x/data/dir/" ### the directory that contains all the results from cellranger for the fibroblast pools
outdir = "/path/to/output/PBMC/simulation" ### the directory that all the results are written to
singlet_dir = "path/to/output/PBMC/Round1Overlap" ### out path used in 'PullSingletsDoublets.R' script located at 'Demuxafy_manuscript/scripts/PBMC/empirical/PullSingletsDoublets.R'
FAI="/path/to/genome.fa.fai" ### GRCh38
FASTA="/path/to/genome.fa" ### GRCh38
gtf = "/path/to/genes/genes.gtf" ### path to gtf file
demuxafy_sif = "/path/to/Demuxafy.sif" ### the path to the singularity image "Demuxafy.sif", provided on zenodo
simulation_scripts_dir = "/path/to/Deuxafy_manuscript/scripts/simulated/simulation" ### path to directory containg the simulation scripts that can be accessed from Zenodo or Github


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
    metadata_updated = metadata_updated[(metadata_updated['Pool']!='size128_SimulatedPool2_unevenN_0.75pctl') &  (metadata_updated['Pool'] !='size8_SimulatedPool3_unevenN_0.5pctl') &  (metadata_updated['Pool'] !='size8_SimulatedPool3_unevenN_0.75pctl') & (metadata_updated['Pool'] != 'size128_SimulatedPool2_unevenN_0.5pctl') & ~(metadata_updated['Pool'].str.contains('unevenN_0.95'))]
    metadata_updated['DoubletRate'] = 0

    for pool in metadata_updated.Pool:
        ## get the number of droplets to get the doublet rate
        with gzip.open(outdir + "/" + pool + "/matrix_out/barcodes.tsv.gz", 'rb') as f:
            for i, l in enumerate(f):
                pass
        metadata_updated.loc[metadata_updated['Pool'].str.match(pool), 'DoubletRate'] = metadata_updated.loc[metadata_updated['Pool'].str.match(pool), 'Ndoublets']/(i + 1)

    metadata_updated.to_csv(outdir + "/updated_metadata_unique.tsv", sep = "\t")


print("demultiplexing softwares")
sample_file = outdir + "/updated_metadata_unique.tsv"
samples = pd.read_csv(sample_file, sep = "\t")


all_files.append(expand(outdir + "/{pool}/popscle/freemuxlet/Individual_genotypes_subset.vcf.gz", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/scSplit/Individual_genotypes_subset.vcf.gz", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/matrix_out/matrix.mtx", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/demuxlet_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/freemuxlet_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/scSplit_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/souporcell_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/vireo_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/DoubletFinder_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/scds_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/scDblFinder_results_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/scDblFinder_known_doublets_results_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/solo_results_temp.txt", pool=samples.Pool))
all_files.append(expand(outdir + "/{pool}/CombinedResults/CombinedDemuxletResults.tsv", pool=samples.Pool))
all_files.append(expand(outdir + "/SimulatedOverlap/{pool}/PoolKeys.rds", pool=samples.Pool)) 
all_files.append(expand(outdir + "/benchmarks/{pool}_dropulation_call.qstat", pool=samples.Pool))
all_files.append(expand(outdir + "/benchmarks/{pool}_demuxalot.qstat", pool=samples.Pool))
all_files.append(expand(outdir + "/SimulatedOverlap/{pool}/single_software_metrics.tsv", pool=samples.Pool)),
all_files.append(expand(outdir + "/SimulatedOverlap/{pool}/doublet_detecting_metrics.tsv", pool=samples.Pool)),
all_files.append(expand(outdir + "/SimulatedOverlap/{pool}/demultiplexing_doublet_detecting_metrics_w_single_softs.tsv", pool=samples.Pool)),
all_files.append(expand(outdir + "/SimulatedOverlap/{pool}/demultiplexing_metrics.tsv", pool=samples.Pool)),



if os.path.exists(outdir + "/DoubletDetection_PASS_FAIL.tsv"):
    print("DoubletDetection_PASS_FAIL.tsv exists")
    DoubletDetection_rules.append(expand(outdir + "/{pool}/CombinedResults/DoubletDetection_temp.txt", pool=samples.Pool))
elif os.path.exists(outdir + "/DoubletDetection_rerun.tsv"):
    print("DoubletDetection_rerun.tsv exists")
    DoubletDetection_rerun_df = pd.read_csv(outdir + "/DoubletDetection_rerun.tsv", sep = "\t")
    dd_list = expand(outdir + "/{pool_DD}/DoubletDetection/DoubletDetection_predicted_doublet.txt", pool_DD=DoubletDetection_rerun_df.Pool)
    DoubletDetection_rules.append(expand(outdir + "/{pool_DD}/DoubletDetection/DoubletDetection_predicted_doublet.txt", pool_DD=DoubletDetection_rerun_df.Pool))
else:
    print("running doublet detection baseline")
    DoubletDetection_rules.append(expand(outdir + "/{pool}/DoubletDetection/DoubletDetection_predicted_doublet.txt", pool=samples.Pool))


# If the scrublet output is present => all the contents are there that are needed to move past 
if os.path.exists(outdir + "/scrublet_fix.tsv"):
    print("scrublet fix exists")
    scrublet_fix_df = pd.read_csv(outdir + "/scrublet_fix.tsv", sep = "\t")
    scrublet_rules.append(expand(outdir + "/{scrub_fix_pool}/scrublet_{fix_pctl}/scrublet_results.txt", zip, scrub_fix_pool=scrublet_fix_df["Pool"], fix_pctl=scrublet_fix_df.Percentile))
    number_doublets = [get_expected_Ndoublet(pool) for pool in scrublet_fix_df.Pool]
    doublet_rate = [get_expected_doubletRate(pool) for pool in scrublet_fix_df.Pool]
    scrublet_fix_df['Doublets'] = number_doublets
    scrublet_fix_df['DoubletRate'] = doublet_rate

if os.path.exists(outdir + "/scrublet_pctl.tsv"):
    print("scrublet pctl runs")
    scrublet_decisions = pd.read_csv(outdir + "/scrublet_pctl.tsv", sep = "\t", dtype = {'Pool': str,"Percentile": str})
    scrublet_rules.append(expand(outdir + "/{pool}/CombinedResults/scrublet_temp.txt", pool=scrublet_decisions.Pool))
else:
    print("scurblet results")
    scrublet_rules.append(expand(outdir + "/{pool}/scrublet_{pctl}/scrublet_results.txt", pool=samples.Pool, pctl=config["percentile"]))




# If the DoubletDetection output is present => all the contents are there that are needed to move past 
if os.path.exists(outdir + "/DoubletDecon_rhops.tsv"):
    print("DoubletDecon_rhops.tsv exists")
    DoubletDecon_decisions = pd.read_csv(outdir + "/DoubletDecon_rhops.tsv", sep = "\t", dtype = {'Pool': str,"rhop": str})
    DoubletDecon_rhops.append(expand(outdir + "/{pool}/CombinedResults/DoubletDecon_temp.txt", pool=DoubletDecon_decisions.Pool))
else:
    print("DoubletDecon_rhops.tsv doesn't exists")
    DoubletDecon_rhops.append(expand(outdir + "/{pool}/DoubletDecon_rhop{rhop}/Final_doublets_groups_DoubletDecon_results.txt", pool=samples.Pool, rhop=config["rhop"]))


if os.path.exists(outdir + "/scrublet_pctl.tsv"):
    if os.path.exists(outdir + "/DoubletDecon_rhops.tsv"):
        if os.path.exists(outdir + "/DoubletDetection_PASS_FAIL.tsv"):
            print("all exist")
            all_files.append(expand(outdir + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv", pool=samples.Pool))
            all_files.append(expand(outdir + "/SimulatedOverlap/{pool}/proportion_most_common_cell_type_per_barcode_single_soft.rds", pool=samples.Pool))
        else:
            print("none exist")
    else:
        print("DoubletDecon doesn't exist")
else:
    print("DoubletDecon doesn't exist")

    


rule all:
    input:
        all_files,
        DoubletDecon_rhops,
        scrublet_rules,
        DoubletDetection_rules




#########################################################
########## CREATE MATRICES FOR SIMULATED POOLS ##########
#########################################################
##### Assign genes to bam reads #####
rule feature_counts:
    input:
        gtf=GTF,
        bam=outdir + "/{pool}/pooled.sorted.bam"
    output:
        counts = temp(outdir + "/{pool}/pooled.sorted.bam.featureCounts.bam"),
        qstat = outdir + "/benchmarks/{pool}_feature_counts.qstat"
    params:
        o=outdir + "/{pool}/gene_assigned",
        q = "-q short.q"
    threads: 10
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 8,
    conda:
        simulation_scripts_dir + "/../../../conda_environments/subread.yaml"
    benchmark:
        outdir + "/benchmarks/{pool}_feature_counts.benchmark"
    shell:
        """
        featureCounts -a {input.gtf} -o {params.o} -R BAM -T {threads} {input.bam}
        qstat -j $JOB_ID > {output.qstat}
        """

##### Sort and index bam #####
rule sort_index:
    input:
        outdir + "/{pool}/pooled.sorted.bam.featureCounts.bam"
    threads: 2
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10,
    output:
        bam = outdir + "/{pool}/assigned_sorted.bam",
        qstat = outdir + "/benchmarks/{pool}_sort_index.qstat",
    params:
        q = "-q short.q",
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
        ancient(outdir + "/{pool}/assigned_sorted.bam")
    output:
        counts = outdir + "/{pool}/counts.tsv.gz",
        qstat = outdir + "/benchmarks/{pool}_umi_tools.qstat",
    threads: 4
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 8
    params: 
        q = "-q short.q",
        sif=demuxafy_sif,
    conda:
        simulation_scripts_dir + "/../../../conda_environments/umi_tools.yaml"
    benchmark:
        outdir + "/benchmarks/{pool}_umi_tools.benchmark"
    shell:
        """
        singularity exec {params.sif} umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --umi-tag=UB --cell-tag=CB --extract-umi-method=tag --cell-tag-split=None -I {input} -S {output.counts}
        qstat -j $JOB_ID > {output.qstat}
        """

##### create matrix files from counts files simulated bams for scrublet #####
rule umi_tools_2_mtx:
    input:
        ancient(outdir + "/{pool}/counts.tsv.gz")
    output:
        matrix = outdir + "/{pool}/matrix_out/matrix.mtx.gz",
        features = outdir + "/{pool}/matrix_out/features.tsv.gz",
        barcodes = outdir + "/{pool}/matrix_out/barcodes.tsv.gz",
        variables = outdir + "/{pool}/matrix_out/variables.txt",
        qstat = outdir + "/benchmarks/{pool}_umi_tools_2_mtx.qstat"
    params:
        q = "-q short.q",
        pool_dir = outdir + "/{pool}/",
        script = simulation_scripts_dir + "/../../../tool_scripts/umi_counts_mtx.R"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 12,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 12,
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


###################################
############# SCSPLIT #############
###################################

###### scSplit Preprocessing ######
rule scSplit_sam_header:
    input:
        bam=ancient(outdir + "/{pool}/assigned_sorted.bam"),
    threads: 8
    output:
        bam = temp(outdir + "/{pool}/scSplit/SAM_header"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_sam_header.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_sam_header.benchmark"
    shell:
        """
        singularity exec {params.sif} samtools view -@ {threads} -H {input.bam} > {output.bam}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_sam_body:
    input:
        bam=ancient(outdir + "/{pool}/assigned_sorted.bam"),
        barcodes= ancient(outdir + "/{pool}/matrix_out/barcodes.tsv"),
    threads: 8
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    output:
        bam = temp(outdir + "/{pool}/scSplit/filtered_SAM_body"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_sam_body.qstat",
    params:
        q = "-q short.q",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_sam_body.benchmark"
    shell:
        """
        singularity exec {params.sif} samtools view -@ {threads} -S -q 10 -F 3844 {input.bam} | grep -F -f {input.barcodes} > {output.bam}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_sam_combine:
    input:
        header=outdir + "/{pool}/scSplit/SAM_header",
        body=outdir + "/{pool}/scSplit/filtered_SAM_body"
    threads: 8
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q",
        sif=demuxafy_sif
    output:
        bam = temp(outdir + "/{pool}/scSplit/filtered.bam"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_sam_combine.qstat"
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_sam_combine.benchmark"
    shell:
        """
        singularity exec {params.sif} cat {input.header} {input.body} | singularity exec {params.sif} samtools view -@ {threads} -b - > {output.bam}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_reassign_RG:
    input:
        ancient(outdir + "/{pool}/scSplit/filtered.bam")
    threads: 4
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 256,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 256,
        java_mem = lambda wildcards, attempt: attempt * 640,
    output:
        bam = temp(outdir + "/{pool}/scSplit/filtered_RG.bam"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_reassign_RG.qstat"
    params:
        q = "-q short.q",
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_reassign_RG.benchmark"
    conda:
        simulation_scripts_dir + "/../../../conda_environments/picard.yaml"
    shell:
        """
        java -Xmx{resources.java_mem}g -Xms{resources.java_mem}g -jar picard.jar AddOrReplaceReadGroups I={input} O={output.bam} RGSM={wildcards.pool} RGLB=MissingLibrary.1 RGPL=ILLUMINA RGPU={wildcards.pool}:MissingLibrary:1:reassigned1
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_rmdupe:
    input:
        bam=outdir + "/{pool}/scSplit/filtered_RG.bam"
    output:
        bam = temp(outdir + "/{pool}/scSplit/dedup_filtered.bam"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_rmdupe.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 20,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 20,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_rmdupe.benchmark"
    shell:
        """
        singularity exec {params.sif} samtools rmdup {input.bam} {output.bam}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_sort:
    input:
        outdir + "/{pool}/scSplit/dedup_filtered.bam"
    threads: 1
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 15,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 15,
    output:
        bam = temp(outdir + "/{pool}/scSplit/possort_dedup_filtered.bam"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_sort.qstat"
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
        out = outdir + "/{pool}/"
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_sort.benchmark"
    shell:
        """
        singularity exec {params.sif} samtools sort -@ {threads} -T {params.out} -o {output.bam} {input}
        singularity exec {params.sif} samtools index {output.bam}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_regions:
    input:
        fai=FAI,
    output:
        regions = temp(outdir + "/{pool}/scSplit/regions_file"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_regions.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 5,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_regions.benchmark"
    shell:
        """
        singularity exec {params.sif} fasta_generate_regions.py {input.fai} 100000 > {output.regions}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_freebayes:
    input:
        fasta=FASTA,
        bam=outdir + "/{pool}/scSplit/possort_dedup_filtered.bam",
        regions=outdir + "/{pool}/scSplit/regions_file"
    output:
        vcf = temp(outdir + "/{pool}/scSplit/freebayes_var.vcf"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_freebayes.qstat"
    threads: 15
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 2,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 2,
    params:
        q = "-q short.q",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_freebayes.benchmark"
    shell:
        """
        export TMPDIR=/tmp
        singularity exec {params.sif} freebayes-parallel {input.regions} {threads} -f {input.fasta} -iXu -C 2 -q 1 {input.bam} > {output.vcf}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_vcf_qual_filt:
    input:
        vcf=outdir + "/{pool}/scSplit/freebayes_var.vcf"
    output:
        vcf = temp(outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf"),
        qstat = outdir + "/benchmarks/{pool}_scSplit_vcf_qual_filt.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 2,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 2,
    threads: 1
    params:
        q = "-q short.q",
        out=outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_vcf_qual_filt.benchmark"
    shell:
        """
        singularity exec {params.sif} vcftools --gzvcf {input.vcf} --minQ 30 --recode --recode-INFO-all --out {params.out}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.vcf} ]]
        echo $?
        """      

rule scSplit_bgzip:
    input:
        ancient(outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf")
    output:
        gz=outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf.gz",
        index=outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf.gz.tbi",
        qstat = outdir + "/benchmarks/{pool}_scSplit_bgzip.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 2,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 2,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_bgzip.benchmark"
    shell:
        """
        singularity exec {params.sif} bgzip  -c {input} > {output.gz}
        singularity exec {params.sif} tabix -p vcf {output.gz}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.index} ]]
        echo $?
        """

##### This is how it should be done but forgot before pileup #####
rule scSplit_subset_vcf:
    input:
        pileup=outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf.gz",
        snps=SNP_GENOTYPES + ".gz"
    output:
        vcf = outdir + "/{pool}/scSplit/frebayes_var_qual30_subset.vcf",
        qstat = outdir + "/benchmarks/{pool}_scSplit_subset_vcf.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_subset_vcf.benchmark"
    shell:
        """
        singularity exec {params.sif} bcftools view {input.pileup} -R {input.snps} -Ov -o {output.vcf}
        qstat -j $JOB_ID > {output.qstat}
        """


##### scSplit Allele Counting #####
rule scSplit_allele_matrices:
    input:
        snvs=SNP_GENOTYPES,
        vcf=outdir + "/{pool}/scSplit/frebayes_var_qual30_subset.vcf",
        bam=ancient(outdir + "/{pool}/scSplit/possort_dedup_filtered.bam"),
        barcodes= ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
    output:
        alt=outdir + "/{pool}/scSplit/alt_filtered.csv",
        ref=outdir + "/{pool}/scSplit/ref_filtered.csv",
        qstat = outdir + "/benchmarks/{pool}_scSplit_allele_matrices.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 24,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 24,
    threads: 4
    params:
        q = "-q short.q",
        out=outdir + "/{pool}/scSplit/",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_allele_matrices.benchmark"
    shell:
        """
        singularity exec {params.sif} scSplit count -c {input.snvs} -v {input.vcf} -i {input.bam} -b {input.barcodes} -r {output.ref} -a {output.alt} -o {params.out}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.alt} ]]
        echo $?
        """

##### scSplit Demultiplexing #####
rule scSplit_demultiplex:
    input:
        alt=outdir + "/{pool}/scSplit/alt_filtered.csv",
        ref=outdir + "/{pool}/scSplit/ref_filtered.csv"
    output:
        Psc=outdir + "/{pool}/scSplit/scSplit_P_s_c.csv",
        result=outdir + "/{pool}/scSplit/scSplit_result.csv",
        qstat = outdir + "/benchmarks/{pool}_scSplit_demultiplex.qstat"
    threads: 10
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 36,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 36,
    params:
        q = "-q short.q",
        out=outdir + "/{pool}/scSplit/",
        sif=demuxafy_sif,
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_demultiplex.benchmark"
    shell:
        """  
        singularity exec {params.sif} scSplit run -r {input.ref} -a {input.alt} -n {params.N} -o {params.out}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.Psc} ]]
        [[ -s {output.result} ]]
        echo $?
        """

##### scSplit Get Genotypes #####
rule scSplit_genotypes:
    input:
        alt=outdir + "/{pool}/scSplit/alt_filtered.csv",
        ref=outdir + "/{pool}/scSplit/ref_filtered.csv",
        demultiplex=outdir + "/{pool}/scSplit/scSplit_P_s_c.csv"
    output:
        vcf = outdir + "/{pool}/scSplit/scSplit.vcf",
        qstat = outdir + "/benchmarks/{pool}_scSplit_genotypes.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 96,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 96,
    threads: 4
    params:
        q = "-q short.q",
        out=outdir + "/{pool}/scSplit/",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_scSplit_genotypes.benchmark"
    shell:
        """
        singularity exec {params.sif} scSplit genotype -r {input.ref} -a {input.alt} -p {input.demultiplex} -o {params.out}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.vcf} ]] 
        echo $?
        """

###################################
############# POPSCLE #############
###################################

###### popscle Preprocessing ######
rule popscle_pileup:
    input:
        vcf=SNP_GENOTYPES,
        bam=ancient(outdir + "/{pool}/assigned_sorted.bam"),
        barcodes= outdir + "/{pool}/matrix_out/barcodes.tsv.gz",
    output:
        qstat = outdir + "/benchmarks/{pool}_popscle_pileup.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 96,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 96,
    threads: 8
    params:
        pileup = directory(outdir + "/{pool}/popscle/pileup/"),
        q = "-q short.q",
        sif=demuxafy_sif
    benchmark:
        outdir + "/benchmarks/{pool}_popscle_pileup.benchmark"
    shell:
        """
        singularity exec {params.sif} popscle dsc-pileup --sam {input.bam} --vcf {input.vcf} --group-list {input.barcodes} --out {params.pileup}pileup
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {params.pileup}pileup.var.gz ]]
        echo $?
        """

##### Popscle Freemuxlet Demultiplexing #####
rule popscle_freemuxlet:
    input:
        qstat = outdir + "/benchmarks/{pool}_popscle_pileup.qstat",
        barcodes= outdir + "/{pool}/matrix_out/barcodes.tsv.gz",
    output:
        results = outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.vcf.gz",
        qstat = outdir + "/benchmarks/{pool}_popscle_freemuxlet.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 96,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 96,
    threads: 8
    params:
        pileup=outdir + "/{pool}/popscle/pileup/",
        # q = "-q short.q",
        q = "-q long.q",
        out=outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT",
        sif=demuxafy_sif,
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "/benchmarks/{pool}_popscle_freemuxlet.benchmark"
    shell:
        """
        singularity exec {params.sif} popscle freemuxlet --plp {params.pileup}pileup --out {params.out} --group-list {input.barcodes} --nsample {params.N}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.results} ]]
        echo $?
        """

##### Popscle Demuxlet Individual File Generation #####
rule popscle_demuxlet_ind_files:
    input:
        outdir + "/{pool}/popscle/pileup/"
    output:
        indivs = outdir + "/{pool}/popscle/Individuals.txt",
        qstat = outdir + "/benchmarks/{pool}_popscle_demuxlet_ind_files.qstat"
    resources:
        mem_per_thread_gb=5,
        disk_per_thread_gb=5,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "/benchmarks/{pool}_popscle_demuxlet_ind_files.benchmark"
    shell:
        """
        singularity exec {params.sif} echo {params.individuals} | tr "," "\n" > {output.indivs}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.indivs} ]]
        echo $?
        """

##### Popscle Demuxlet Demultiplexing #####
rule popscle_demuxlet:
    input:
        pileup=outdir + "/{pool}/popscle/pileup/",
        snps=SNP_GENOTYPES,
        barcodes= ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
        individuals=ancient(outdir + "/{pool}/popscle/Individuals.txt"),
        qstat = ancient(outdir + "/benchmarks/{pool}_popscle_pileup.qstat"),
    output:
        best = outdir + "/{pool}/popscle/demuxlet/demuxletOUT_impute_vars.best",
        qstat = outdir + "/benchmarks/{pool}_popscle_demuxlet.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 128,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 128,
    threads: 4
    params:
        q = "-q short.q",
        out=outdir + "/{pool}/popscle/demuxlet/",
        sif=demuxafy_sif,
        field="GP"
    benchmark:
        outdir + "/benchmarks/{pool}_popscle_demuxlet.benchmark"
    shell:
        """
        singularity exec {params.sif} popscle demuxlet --plp {input.pileup}pileup --vcf {input.snps} --field {params.field} --group-list {input.barcodes} --geno-error-coeff 1.0 --geno-error-offset 0.05 --out {params.out}demuxletOUT_impute_vars --sm-list {input.individuals}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.best} ]]
        echo $?
        """

###################################
############## VIREO ##############
###################################

####### vireo Preprocessing #######
rule cellSNP:
    input:
        vcf=SNP_GENOTYPES,
        bam=ancient(outdir + "/{pool}/assigned_sorted.bam"),
        barcodes= ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
    output:
        vcf = outdir + "/{pool}/vireo/cellSNP.base.vcf.gz",
        qstat = outdir + "/benchmarks/{pool}_cellSNP.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 24,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 24,
    threads: 2
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
        p=20,
        maf=0.1,
        counts=20,
        out = outdir + "/{pool}/vireo/"
    benchmark:
        outdir + "/benchmarks/{pool}_cellSNP.benchmark"
    shell:
        """
        singularity exec {params.sif} cellsnp-lite -s {input.bam} -b {input.barcodes} -O {params.out} -R {input.vcf} -p {params.p} --minMAF {params.maf} --minCOUNT {params.counts} --gzip
        qstat -j $JOB_ID > {output.qstat}
        """

##### Subset the imputed genotype files by the individuals in the pools #####
rule subset_vcf:
    input:
        pileup=outdir + "/{pool}/vireo/cellSNP.base.vcf.gz",
        snps=SNP_GENOTYPES + ".gz"
    output:
        vcf = outdir + "/{pool}/vireo/Merged_MAF0.01.dose_GeneFiltered_hg38_individualSubset.vcf.gz",
        qstat = outdir + "/benchmarks/{pool}_subset_vcf.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 8,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "/benchmarks/{pool}_subset_vcf.benchmark"
    shell:
        """
        singularity exec {params.sif} bcftools view -R {input.pileup} -s {params.individuals} -Oz -o {output.vcf} {input.snps}
        qstat -j $JOB_ID > {output.qstat}
        """

##### Vireo demultiplexing #####
rule vireo:
    input:
        pileup= outdir + "/{pool}/vireo/cellSNP.base.vcf.gz",
        snps=outdir + "/{pool}/vireo/Merged_MAF0.01.dose_GeneFiltered_hg38_individualSubset.vcf.gz"
    output:
        tsv = outdir + "/{pool}/vireo/results/donor_ids.tsv",
        qstat = outdir + "/benchmarks/{pool}_vireo.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 36,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 36,
    threads: 5
    params:
        q = "-q short.q",
        # q = "-q long.q",
        out=outdir + "/{pool}/vireo/results/",
        sif=demuxafy_sif,
        field="GP",
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0],
        cellsnp_dir = outdir + "/{pool}/vireo/"
    benchmark:
        outdir + "/benchmarks/{pool}_vireo.benchmark"
    shell:
        """
        singularity exec {params.sif} vireo -c {params.cellsnp_dir} -d {input.snps} -o {params.out} -t {params.field} -N {params.N} --callAmbientRNAs
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.tsv} ]]
        echo $?
        """


####################################
############ SOUPORCELL ############
####################################

###### souporcell pipeline ########
rule souporcell:
    input:
        bam=ancient(outdir + "/{pool}/assigned_sorted.bam"),
        barcodes= ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
        fasta=FASTA,
        snps=SNP_GENOTYPES
    threads: 8
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 192,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 192,
    output:
        clusters=outdir + "/{pool}/souporcell/clusters.tsv",
        qstat = outdir + "/benchmarks/{pool}_souporcell.qstat",
        cluster_geno=outdir + "/{pool}/souporcell/cluster_genotypes.vcf",
    params:
        q = "-q short.q",
        out=outdir + "/{pool}/souporcell/",
        sif=demuxafy_sif,
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0],
    benchmark:
        outdir + "/benchmarks/{pool}_souporcell.benchmark"
    shell:
        """
        singularity exec {params.sif} souporcell_pipeline.py -i {input.bam} -b {input.barcodes} -f {input.fasta} -t {threads} -o {params.out} -k {params.N} --common_variants {input.snps}
        rm {params.out}/souporcell_minimap_tagged_sorted.bam
        rm {params.out}/souporcell_minimap_tagged_sorted.bam.bai
        qstat -j $JOB_ID >> {output.qstat}
        [[ -s {output.clusters} ]]
        echo $?
        """


#####################################
############ DROPULATION ############
#####################################
rule dropulation_tag:
    input:
        gtf = gtf,
        bam = ancient(outdir + "/{pool}/assigned_sorted.bam")
    threads: 8
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt *248,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 96,
        java_mem=lambda wildcards, attempt: attempt * 96
    benchmark:
        outdir + "/benchmarks/{pool}.dropulation_tag.txt"
    output:
        bam = temp(outdir + "/{pool}/dropulation/tagged_bam.bam"),
        qstat = outdir + "/benchmarks/{pool}_dropulation_tag.qstat"
    params:
        sif=demuxafy_sif,
        q = "-q short.q",
    benchmark:
        outdir + "/benchmarks/{pool}_dropulation_tag.benchmark"
    shell:
        """
        export _JAVA_OPTIONS="-Xmx{resources.java_mem}g -XX:+UseCompressedOops -XX:MaxPermSize=1G"

        singularity exec {params.sif} TagReadWithGeneFunction \
            -m {resources.java_mem}g \
            --ANNOTATIONS_FILE {input.gtf} \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} 

        qstat -j $JOB_ID > {output.qstat}
        """


rule dropulation_assign:
    input:
        snps=SNP_GENOTYPES,
        bam=outdir + "/{pool}/dropulation/tagged_bam.bam",
        barcodes= outdir + "/{pool}/matrix_out/barcodes.tsv",
        individuals = outdir + "/{pool}/popscle/Individuals.txt"
    threads: 16
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        java_mem=lambda wildcards, attempt: attempt * 8
    benchmark:
        outdir + "/benchmarks/{pool}_dropulation_assign.benchmark"
    output:
        qstat = outdir + "/benchmarks/{pool}.dropulation_assign.qstat",
        vcf = outdir + "/{pool}/dropulation/out_vcf.vcf",
        assignments = outdir + "/{pool}/dropulation/assignments.tsv.gz"
    params:
        sif=demuxafy_sif,
        q = "-q short.q",
    shell:
        """
        export _JAVA_OPTIONS="-Xmx{resources.java_mem}g"

        singularity exec {params.sif} AssignCellsToSamples --CELL_BC_FILE {input.barcodes} \
            --INPUT_BAM {input.bam} \
            --OUTPUT {output.assignments} \
            --VCF {input.snps} \
            --SAMPLE_FILE {input.individuals} \
            --CELL_BARCODE_TAG 'CB' \
            --MOLECULAR_BARCODE_TAG 'UB' \
            --VCF_OUTPUT {output.vcf} \
            --MAX_ERROR_RATE 0.05
        qstat -j $JOB_ID > {output.qstat}
        """


rule dropulation_doublet:
    input:
        assignments = outdir + "/{pool}/dropulation/assignments.tsv.gz",
        bam = outdir + "/{pool}/dropulation/tagged_bam.bam",
        barcodes = outdir + "/{pool}/matrix_out/barcodes.tsv",
        vcf = outdir + "/{pool}/dropulation/out_vcf.vcf",
        individuals = outdir + "/{pool}/popscle/Individuals.txt",
    threads: 16
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        java_mem=lambda wildcards, attempt: attempt * 8
    benchmark:
        outdir + "/benchmarks/{pool}_dropulation_doublet.benchmark"
    output:
        qstat = outdir + "/benchmarks/{pool}.dropulation_doublet.qstat",
        likelihoods = outdir + "/{pool}/dropulation/likelihoods.tsv.gz"
    params:
        sif=demuxafy_sif,
        q = "-q short.q",
    shell:
        """
        export _JAVA_OPTIONS="-Xmx{resources.java_mem}g"

        singularity exec {params.sif} DetectDoublets --CELL_BC_FILE {input.barcodes} \
            --INPUT_BAM {input.bam} \
            --OUTPUT {output.likelihoods} \
            --VCF {input.vcf} \
            --CELL_BARCODE_TAG 'CB' \
            --MOLECULAR_BARCODE_TAG 'UB' \
            --SINGLE_DONOR_LIKELIHOOD_FILE {input.assignments} \
            --SAMPLE_FILE {input.individuals} \
            --MAX_ERROR_RATE 0.05

        qstat -j $JOB_ID > {output.qstat}
        """


rule dropulation_call:
    input:
        assignments = outdir + "/{pool}/dropulation/assignments.tsv.gz",
        likelihoods = outdir + "/{pool}/dropulation/likelihoods.tsv.gz"
    threads: 16
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 4,
    benchmark:
        outdir + "/benchmarks/{pool}_dropulation_doublet.benchmark"
    output:
        results = outdir + "/{pool}/dropulation/updated_assignments.tsv.gz",
        qstat = outdir + "/benchmarks/{pool}_dropulation_call.qstat"
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
        script = simulation_scripts_dir + "/../../../tool_scripts/dropulation_call.R"
    shell:
        """
        Rscript {params.script} {input.assignments} {input.likelihoods} {output}
        qstat -j $JOB_ID > {output.qstat}
        """


###################################
############ DEMUXALOT ############
###################################
rule demuxalot:
    input:
        bam = ancient(outdir + "/{pool}/assigned_sorted.bam"),
        barcodes = ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
        fasta=FASTA,
        snps=ancient(SNP_GENOTYPES),
        individuals = ancient(outdir + "/{pool}/popscle/Individuals.txt"),
    threads: 8
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 32,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 32,
    benchmark:
        outdir + "/benchmarks/{pool}_demuxalot.benchmark"
    output:
        qstat = outdir + "/benchmarks/{pool}_demuxalot.qstat",
        vcf = outdir + "/{pool}/demuxalot/assignments.tsv.gz",
        vcf_refined = outdir + "/{pool}/demuxalot/assignments_refined.tsv.gz"
    params:
        q = "-q short.q",
        out = outdir + "/{pool}/demuxalot/"
    shell:
        """
        singularity exec {params.sif} Demuxalot.py \
            -b {input.barcodes} \
            -a {input.bam} \
            -n {input.individuals} \
            -v {input.snps} \
            -o {params.out}

        qstat -j $JOB_ID > {output.qstat}
        """


rule join_demultiplexing_results:
    input:
        demuxlet=outdir + "/{pool}/CombinedResults/demuxlet_temp.txt",
        freemuxlet=outdir + "/{pool}/CombinedResults/freemuxlet_temp.txt",
        scSplit=ancient(outdir + "/{pool}/CombinedResults/scSplit_temp.txt"),
        souporcell=ancient(outdir + "/{pool}/CombinedResults/souporcell_temp.txt"),
        vireo=outdir + "/{pool}/CombinedResults/vireo_temp.txt",
        dropulation = outdir + "/{pool}/CombinedResults/dropulation_temp.txt",
        demuxalot = outdir + "/{pool}/CombinedResults/demuxalot_temp.txt",
        demuxalot_refined = outdir + "/{pool}/CombinedResults/demuxalot_refined_temp.txt",
    output:
        outdir + "/{pool}/CombinedResults/CombinedDemuxletResults.tsv",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 4,
    threads: 1
    params:
        q = "-q short.q",
        sif = demuxafy_sif,
    shell:
        """
        singularity exec {params.sif} join -a1 -a2 -1 1 -2 1  -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7" {input.demuxlet} {input.freemuxlet} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
        singularity exec {params.sif} join -a1 -a2 -1 1 -2 1  -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3" - {input.scSplit} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
        singularity exec {params.sif} join -a1 -a2 -1 1 -2 1  -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.2,2.3,2.4,2.5" - {input.souporcell} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
        singularity exec {params.sif} join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.2,2.3,2.4,2.5,2.6" - {input.vireo} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
        singularity exec {params.sif} join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,2.2,2.3,2.4,2.5,2.6" - {input.dropulation} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
        singularity exec {params.sif} join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,2.2,2.3" - {input.demuxalot} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
        singularity exec {params.sif} join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,2.2,2.3" - {input.demuxalot_refined} | sed "s/ /\t/g" > {output}
        """
    

rule demultiplexing_doublets:
    input:
        demultiplex = ancient(outdir + "/{pool}/CombinedResults/CombinedDemuxletResults.tsv"),
    output:
        scDblFinder = outdir + "/{pool}/scDblFinder_known_doublets/scDblFinder_known_doublets.tsv",
        variables = outdir + "/{pool}/variables.tsv",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 4,
    threads: 1
    params:
        script = simulation_scripts_dir + "/demultiplex_doublets.R",
        q = "-q short.q",
        sif = demuxafy_sif,
        basedir = outdir,
        pool = "{pool}",
        outdir = outdir + "/{pool}/"
    conda:
        simulation_scripts_dir + "/../../../conda_environments/generalR.yaml"
    shell:
        """
        singularity exec {params.sif} echo {params.basedir} > {output.variables}
        singularity exec {params.sif} echo {params.pool} >> {output.variables}
        singularity exec {params.sif} echo {params.outdir} >> {output.variables}
        singularity exec {params.sif} Rscript {input.script} {output.variables}
        """


##################################
############ SCRUBLET ############
##################################
##### unzip files for scrublet #####
rule unzip:
    input:
        matrix = ancient(outdir + "/{pool}/matrix_out/matrix.mtx.gz"),
        genes = ancient(outdir + "/{pool}/matrix_out/features.tsv.gz"),
        barcodes = ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz")
    output:
        matrix = outdir + "/{pool}/matrix_out/matrix.mtx",
        genes = outdir + "/{pool}/matrix_out/genes.tsv",
        barcodes = outdir + "/{pool}/matrix_out/barcodes.tsv",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 4,
    threads: 2
    params:
        q = "-q short.q",
        sif=demuxafy_sif
    shell:
        """
        singularity exec {params.sif} gunzip -c {input.matrix} > {output.matrix}
        singularity exec {params.sif} gunzip -c {input.genes} > {output.genes}
        singularity exec {params.sif} gunzip -c {input.barcodes} > {output.barcodes}
        """

if os.path.exists(outdir + "/scrublet_fix.tsv"):
    rule scrublet_fix:
        input:
            fix = outdir + "/scrublet_fix.tsv",
            matrix=outdir + "/{scrub_fix_pool}/matrix_out/matrix.mtx",
            barcodes=outdir + "/{scrub_fix_pool}/matrix_out/barcodes.tsv"
        output:
            scrub = outdir + "/{scrub_fix_pool}/scrublet_{fix_pctl}/scrublet_results.txt",
            qstat = outdir + "/benchmarks/{scrub_fix_pool}_{fix_pctl}pct_scrublet_fix.qstat"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 16,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 16,
        threads: 1
        params:
            q = "-q short.q",
            thresh=lambda wildcards: scrublet_fix_df.Threshold[scrublet_fix_df.Pool == wildcards.scrub_fix_pool].iloc[0],
            sif = demuxafy_sif,
            out=outdir + "/{scrub_fix_pool}/scrublet_{fix_pctl}/",
            pipeline_dir =  simulation_scripts_dir + "/../../../tool_scripts/","
            script = simulation_scripts_dir + "/../../../tool_scripts/scrublet_pipeline.py",
            doublet_rate=lambda wildcards: scrublet_fix_df.DoubletRate[scrublet_fix_df.Pool == wildcards.scrub_fix_pool].iloc[0]
        benchmark:
            outdir + "/benchmarks/{scrub_fix_pool}_{fix_pctl}pct_scrublet_fix.benchmark"
        shell:
            """
            singularity exec {params.sif} python {params.script} --counts_matrix {input.matrix} --barcodes {input.barcodes} --min_gene_variability_pctl {wildcards.fix_pctl} -o {params.out} -d {params.pipeline_dir} -t {params.thresh} --dbl_rate {params.doublet_rate}
            [[ -s {output.scrub} ]]
            qstat -j $JOB_ID > {output.qstat}
            echo $?
            """

else:
    ###### scrublet pipeline ########
    rule scrublet:
        input:
            matrix= ancient(outdir + "/{pool}/matrix_out/matrix.mtx"),
            barcodes= ancient(outdir + "/{pool}/matrix_out/barcodes.tsv")
        output:
            scrub = outdir + "/{pool}/scrublet_{pctl}/scrublet_results.txt",
            qstat = outdir + "/benchmarks/{pool}_{pctl}pct_scrublet.qstat"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 16,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 16,
        threads: 2
        params:
            q = "-q short.q",
            qstat = outdir + "/benchmarks/{pool}_scrublet_{pctl}pctl_qstat.txt",
            pctl="{pctl}",
            sif=demuxafy_sif,
            out=outdir + "/{pool}/scrublet_{pctl}/",
            pipeline_dir =  simulation_scripts_dir + "/../../../tool_scripts/","
            script = simulation_scripts_dir + "/../../../tool_scripts/scrublet_pipeline.py",
            doublet_rate=lambda wildcards: samples.DoubletRate[samples.Pool == wildcards.pool].iloc[0]
        benchmark:
            outdir + "/benchmarks/{pool}_{pctl}pct_scrublet.benchmark"
        shell:
            """
            singularity exec {params.sif} python {params.script} --counts_matrix {input.matrix} --barcodes {input.barcodes} --min_gene_variability_pctl {params.pctl} -o {params.out} -d {params.pipeline_dir} --dbl_rate {params.doublet_rate}
            qstat -j $JOB_ID > {output.qstat}
            [[ -s {output.scrub} ]]
            echo $?
            """


rule scds:
    input:
        matrix=ancient(outdir + "/{pool}/matrix_out/matrix.mtx.gz")
    output: 
        doublets=outdir + "/{pool}/scds/scds_doublets.txt",
        variables=outdir + "/{pool}/scds/scds_variables.txt",
        qstat = outdir + "/benchmarks/{pool}_scds.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 48,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 48,
    threads: 2
    params:
    
        script = simulation_scripts_dir + "/../../../tool_scripts/scds.R",
        q = "-q short.q",
        out=outdir + "/{pool}/scds/",
        sif=demuxafy_sif,
        matrix_dir = outdir + "/{pool}/matrix_out/"
    benchmark:
        outdir + "/benchmarks/{pool}_scds.benchmark"
    shell:
        """
        singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} Rscript {input.script} {output.variables}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.doublets} ]]
        echo $?
        """


#######################################
############ DOUBLET DECON ############
#######################################
rule DoubletDecon:
    input:
        barcodes=ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz")
    output:
        doublets = outdir + "/{pool}/DoubletDecon_rhop{rhop}/Final_doublets_groups_DoubletDecon_results.txt",
        singlets = outdir + "/{pool}/DoubletDecon_rhop{rhop}/Final_nondoublets_groups_DoubletDecon_results.txt",
        variables=outdir + "/{pool}/DoubletDecon_rhop{rhop}/DoubletDecon_variables.txt",
        qstat = outdir + "/benchmarks/{pool}_rhop{rhop}_DoubletDecon.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 512,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 512,
    threads: 2
    params:
        q = "-q short.q",
        script=simulation_scripts_dir + "/../../../tool_scripts/DoubletDecon.R",
        matrix_dir=outdir + "/{pool}/matrix_out/",
        out=outdir + "/{pool}/DoubletDecon_rhop{rhop}/",
        sif=demuxafy_sif,
        res = "0.2"
    benchmark:
        outdir + "benchmarks/{pool}_rhop{rhop}_DoubletDecon.benchmark"
    shell:
        """
        singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} echo {wildcards.rhop} >> {output.variables}
        singularity exec {params.sif} echo {params.res} >> {output.variables}
        singularity exec {params.sif} Rscript {params.script} {output.variables}
        qstat -j $JOB_ID > {output.qstat}
        """

########################################
############ DOUBLET FINDER ############
########################################
rule DoubletFinder:
    input:
        matrix = ancient(outdir + "/{pool}/matrix_out/matrix.mtx")
    output:
        doublets = outdir + "/{pool}/DoubletFinder/DoubletFinder_doublets.txt",
        variables = outdir + "/{pool}/DoubletFinder/DoubletFinder_variables.txt",
        qstat = outdir + "/benchmarks/{pool}_DoubletFinder.qstat"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 160,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 160,
    threads: 2
    params:
        q = "-q short.q",
        qstat = outdir + "/benchmarks/{pool}_DoubletFinder_qstat.txt",
        matrix_dir = outdir + "/{pool}/matrix_out/",
        out = outdir + "/{pool}/DoubletFinder/",
        script = simulation_scripts_dir + "/../../../tool_scripts/DoubletFinder.R",
        sif=demuxafy_sif,
        QCdir = outdir + "/{pool}/QCfiltered_gene_bc_matrices/",
        ndbl = lambda wildcards: samples.Ndoublets[samples.Pool == wildcards.pool].iloc[0],
        filedir = tool_scripts + "../files/"
    benchmark:
        outdir + "benchmarks/{pool}_DoubletFinder.benchmark"
    shell:
        """
        singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.QCdir} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} echo {params.ndbl} >> {output.variables}
        singularity exec {params.sif} echo {params.filedir} >> {output.variables}
        singularity exec {params.sif} Rscript {params.script} {output.variables}
        qstat -j $JOB_ID > {output.qstat}
        [[ -s {output.doublets} ]]
        echo $?
        """

###########################################
############ DOUBLET DETECTION ############
###########################################
if os.path.exists(outdir + "/DoubletDetection_rerun.tsv"):
    iterations = 150
else:
    iterations = 50

rule DoubletDetection:
    input:
        matrix=ancient(outdir + "/{pool}/matrix_out/matrix.mtx.gz")
    output:
        doublets = outdir + "/{pool}/DoubletDetection/DoubletDetection_predicted_doublet.txt",
        variables= outdir + "/{pool}/DoubletDetection/DoubletDetection_variables.txt",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 96,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 96,
    threads: 2
    params:
        qstat = outdir + "/benchmarks/{pool}_DoubletDetection.qstat",
        q = "-q short.q",
        matrix_dir=outdir + "/{pool}/matrix_out/",
        out=outdir + "/{pool}/DoubletDetection/",
        script=simulation_scripts_dir + "/../../../tool_scripts/DoubletDetection.py",
        sif=demuxafy_sif,
        iterations = iterations
    benchmark:
        outdir + "benchmarks/{pool}_DoubletDetection.benchmark"
    shell:
        """
        singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} python {params.script} --counts_matrix {input.matrix} -o {params.out} --n_iterations {params.iterations}
        ls -l {output.doublets}
        qstat -j $JOB_ID > {params.qstat}
        [[ -s {output.doublets} ]]
        echo $?
        """


############################################
################### SOLO ###################
############################################
rule make_anndata:
    input:
        matrix = outdir + "/{pool}/matrix_out/matrix.mtx",
        barcodes = outdir + "/{pool}/matrix_out/barcodes.tsv",
        genes = outdir + "/{pool}/matrix_out/genes.tsv"
    output:
        outdir + "/{pool}/matrix_out/anndata.h5ad"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 8,
    threads: 1
    params:
        q = "-q short.q",
        sif = demuxafy_sif,
        matrix_dir = outdir + "/{pool}/matrix_out/",
        script = simulation_scripts_dir + "/../../../tool_scripts/create_anndata4solo.py"
    benchmark:
        outdir + "/benchmarks/{pool}.make_anndata.txt"
    shell:
        """
        singularity exec {params.sif} python {params.script} -m {input.matrix} -b {input.barcodes} -g {input.genes} -d {params.matrix_dir}
        """

rule solo:
    input:
        ancient(outdir + "/{pool}/matrix_out/anndata.h5ad")
    output:
        results = outdir + "/{pool}/solo/is_doublet.npy",
        score = outdir + "/{pool}/solo/logit_scores.npy",
        qstat = outdir + "/benchmarks/{pool}_solo.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 16,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 16,
    threads: 8
    params:
        # q = "-l nvgpu=1",
        q = "-l short.q",
        indir = outdir + "/{pool}/matrix_out/",
        out=outdir + "/{pool}/solo/",
        doublets=lambda wildcards: samples.Ndoublets[samples.Pool == wildcards.pool].iloc[0],
        sif=demuxafy_sif,
        json = simulation_scripts_dir + "/../../../tool_scripts/solo_model.json", 
        qstat = outdir + "/benchmarks/{pool}_solo_qstat.txt"
    benchmark:
        outdir + "/benchmarks/{pool}_solo.benchmark"
    shell:
        """
        singularity exec --nv {params.sif} solo -o {params.out} -e {params.doublets} -j {params.json} -d {input} -g
        qstat -j $JOB_ID > {output.qstat}
        """

rule solo_convert_results:
    input:
        doublets = outdir + "/{pool}/solo/is_doublet.npy",
        score = outdir + "/{pool}/solo/logit_scores.npy",
        barcodes = outdir + "/{pool}/matrix_out/barcodes.tsv"
    output:
        result = outdir + "/{pool}/solo/solo_results.txt",
        qstat = outdir + "/benchmarks/{pool}_solo_convert_results.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 4,
    threads: 1
    params:
        q = "-q short.q",
        out = outdir + "/{pool}/solo",
        solo = outdir + "/{pool}/solo",
        sif=demuxafy_sif,
        script = simulation_scripts_dir + "/../../../tool_scripts/Convert_solo_output.py"
    benchmark:
        outdir + "/benchmarks/{pool}_solo_convert_results.benchmark"
    shell:
        """
        singularity exec {params.sif} python {params.script} -b {input.barcodes} -s {params.solo} -o {params.out}
        qstat -j $JOB_ID > {output.qstat}
        """

###################################################
################### SCDBLFINDER ###################
###################################################
rule scdblfinder:
    input:
        matrix = ancient(outdir + "/{pool}/matrix_out/matrix.mtx")
    output:
        results = outdir + "/{pool}/scDblFinder/scDblFinder_results.txt",
        variables = outdir + "/{pool}/scDblFinder/scDblFinder_variables.txt",
        qstat = outdir + "/benchmarks/{pool}_scdblfinder.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 12,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 12,
    threads: 2
    params:
        q = "-q short.q",
        out = outdir + "/{pool}/scDblFinder/",
        matrix_dir = outdir + "/{pool}/matrix_out/",
        script = simulation_scripts_dir + "/../../../tool_scripts/scDblFinder.R",
        sif = demuxafy_sif,
        doublet_ratio = lambda wildcards: samples.DoubletRate[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "benchmarks/{pool}_scdblfinder.benchmark"
    shell:
        """
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} echo {params.doublet_ratio} >> {output.variables}
        singularity exec {params.sif} Rscript {params.script} {output.variables}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scdblfinder_known_doublets:
    input:
        matrix = ancient(outdir + "/{pool}/matrix_out/matrix.mtx"),
        known_doublets = ancient(outdir + "/{pool}/scDblFinder_known_doublets/scDblFinder_known_doublets.tsv")
    output:
        results = outdir + "/{pool}/scDblFinder_known_doublets/scDblFinder_results.txt",
        variables = outdir + "/{pool}/scDblFinder_known_doublets/scDblFinder_variables.txt",
        qstat = outdir + "/benchmarks/{pool}_scdblfinder_known_doublets.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 8,
    threads: 2
    params:
        q = "-q short.q",
        qstat = outdir + "/benchmarks/{pool}_scdblfinder_known_doublets_qstat.txt",
        out=outdir + "/{pool}/scDblFinder_known_doublets/",
        matrix_dir= outdir + "/{pool}/matrix_out/",
        script = simulation_scripts_dir + "/../../../tool_scripts/scDblFinder.R",
        sif=demuxafy_sif,
        doublet_ratio = lambda wildcards: samples.DoubletRate[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "benchmarks/{pool}_scdblfinder_known_doublets.benchmark"
    shell:
        """
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} echo {params.doublet_ratio} >> {output.variables}
        singularity exec {params.sif} echo {input.known_doublets} >> {output.variables}
        singularity exec {params.sif} Rscript {params.script} {output.variables}
        qstat -j $JOB_ID > {output.qstat}
        """


#################################
############ COMBINE ############
#################################
rule demuxlet_results_temp:
    input:
        demuxlet=outdir + "/{pool}/popscle/demuxlet/demuxletOUT_impute_vars.best"
    output:
        demuxlet_temp=outdir + "/{pool}/CombinedResults/demuxlet_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    threads: 1
    params:
        q = "-q short.q"
    shell:
        """
        awk 'BEGIN{{OFS=FS="\\t"}}{{print $2,$3,$5,$6,$14,$19,$20}}' {input.demuxlet} | sed "s/SNG/singlet/g" | sed "s/DBL/doublet/g" | awk 'BEGIN{{FS=OFS="\t"}} $3=="doublet" {{$4="doublet"}}1' | sed -E "s/,[0-9]+_[0-9]+,[0-9].[0-9]+\t/\t/g" | sed "s/NUM.SNPS/nSNP/g" | sed "s/DROPLET.TYPE/DropletType/g" | sed "s/BEST.GUESS/Assignment/g" | sed "s/singlet.BEST.LLK/SingletLLK/g" | sed "s/doublet.BEST.LLK/DoulbetLLK/g" | sed "s/DIFF.LLK.singlet.doublet/DiffLLK/g" | sed "1s/\t/\tdemuxlet_/g" | sed "s/BARCODE/Barcode/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}'  > {output.demuxlet_temp}
        """

rule freemuxlet_results_temp:
    input:
        freemuxlet = outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.vcf.gz"
    output:
        freemuxlet_temp = outdir + "/{pool}/CombinedResults/freemuxlet_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1
    threads: 1
    params:
        q = "-q short.q"
    shell:
        """
        gunzip -c {input.freemuxlet} | awk 'BEGIN{{OFS=FS="\\t"}}{{print $2,$3,$5,$6,$14,$19,$20 }}' | sed "s/SNG/singlet/g" | sed "s/DBL/doublet/g" | awk 'BEGIN{{FS=OFS="\\t"}} $3=="doublet" {{$4="doublet"}}1' | sed -E "s/,[0-9]+\t/\t/g" | sed "s/NUM.SNPS/nSNP/g" | sed "s/DROPLET.TYPE/DropletType/g" | sed "s/BEST.GUESS/Assignment/g" | sed "s/singlet.BEST.LLK/SingletLLK/g" | sed "s/doublet.BEST.LLK/DoulbetLLK/g" | sed "s/DIFF.LLK.singlet.doublet/DiffLLK/g" | sed "s/BARCODE/Barcode/g" | sed "1s/\t/\tfreemuxlet_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.freemuxlet_temp}
        """

rule scSplit_results_temp:
    input:
        scSplit=ancient(outdir + "/{pool}/scSplit/scSplit_result.csv")
    output:
        scSplit_temp=outdir + "/{pool}/CombinedResults/scSplit_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1
    threads: 1
    params:
        q = "-q short.q"
    shell:
        """
        sed -E 's/\tDBL-[0-9]+/\tdoublet\tdoublet/g' {input.scSplit} | sed 's/SNG-/singlet\t/g' | sed 's/Cluster/DropletType\tAssignment/g' | sed "1s/\t/\tscSplit_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.scSplit_temp}
        """

rule souporcell_results_temp:
    input:
        souporcell=outdir + "/{pool}/souporcell/clusters.tsv"
    output:
        souporcell_temp=outdir + "/{pool}/CombinedResults/souporcell_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q"
    threads: 1
    shell:
        """
        awk 'BEGIN{{OFS=FS="\\t"}}{{print $1,$2,$3,$4,$5}}' {input.souporcell} | awk 'BEGIN{{FS=OFS="\t"}} $2=="doublet" {{$3="doublet"}}1' | awk 'BEGIN{{FS=OFS="\t"}} $2=="unassigned" {{$4="unassigned"}}1' | sed "s/status/DropletType/g" | sed "s/assignment/Assignment/g" | sed "s/log_prob_singleton/LogProbSinglet/g" | sed "s/log_prob_doublet/LogProbDoublet/g" | sed "s/barcode/Barcode/g" | sed "1s/\t/\tsouporcell_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.souporcell_temp}
        """

rule vireo_results_temp:
    input:
        vireo=outdir + "/{pool}/vireo/results/donor_ids.tsv"
    output:
        vireo_temp=outdir + "/{pool}/CombinedResults/vireo_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q"
    threads: 1
    shell:
        """
        awk 'BEGIN{{OFS=FS="\\t"}}{{print $1,$2,$2,$3,$4,$5}}' {input.vireo} | awk 'BEGIN{{FS=OFS="\\t"}}{{gsub("donor[0-9]+","singlet",$3)}}1' | awk 'BEGIN{{FS=OFS="\\t"}}{{gsub("[0-9]+_[0-9]+","singlet",$3)}}1' | sed "s/donor_id\tdonor_id/Assignment\tDropletType/g" | sed "s/prob_max/ProbSinglet/g" | sed "s/prob_doublet/ProbDoublet/g" | sed "s/n_vars/nSNP/g" | sed "s/cell/Barcode/g" | sed "1s/\t/\tvireo_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.vireo_temp}
        """


rule dropulation_results_temp:
    input:
        dropulation = outdir + "/{pool}/dropulation/updated_assignments.tsv.gz",
    output:
        dropulation_temp=outdir + "/{pool}/CombinedResults/dropulation_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q"
    threads: 1
    shell:
        """
        gunzip -c {input.dropulation} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}'  > {output.dropulation_temp}
        """



rule demuxalot_results_temp:
    input:
        demuxalot=outdir + "/{pool}/demuxalot/assignments.tsv.gz"
    output:
        demuxalot_temp=outdir + "/{pool}/CombinedResults/demuxalot_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q"
    threads: 1
    benchmark:
        outdir + "benchmarks/{pool}.demuxalot_results_temp.txt"
    shell:
        """
        gunzip -c {input.demuxalot} | awk 'BEGIN{{OFS=FS="\\t"}}{{print $1, $2, $2}}' |  awk 'BEGIN{{FS=OFS="\\t"}}{{gsub("[0-9]+_[0-9]+\\\\+[0-9]+_[0-9]+","doublet",$3)}}1' | awk 'BEGIN{{FS=OFS="\\t"}}{{gsub("[0-9]+_[0-9]+","singlet",$3)}}1' | sed "s/BARCODE/Barcode/g" | awk 'BEGIN{{FS=OFS="\\t"}} $3=="doublet" {{$2="doublet"}}1' | sed "s/0\t0/Assignment\tDropletType/g" | sed "1s/\t/\tdemuxalot_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.demuxalot_temp}
        """


rule demuxalot_refined_results_temp:
    input:
        demuxalot=outdir + "/{pool}/demuxalot/assignments_refined.tsv.gz"
    output:
        demuxalot_temp=outdir + "/{pool}/CombinedResults/demuxalot_refined_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q"
    threads: 1
    benchmark:
        outdir + "benchmarks/{pool}.demuxalot_refined_results_temp.txt"
    shell:
        """
        gunzip -c {input.demuxalot} | awk 'BEGIN{{OFS=FS="\\t"}}{{print $1, $2, $2}}' |  awk 'BEGIN{{FS=OFS="\\t"}}{{gsub("[0-9]+_[0-9]+\\\\+[0-9]+_[0-9]+","doublet",$3)}}1' | awk 'BEGIN{{FS=OFS="\\t"}}{{gsub("[0-9]+_[0-9]+","singlet",$3)}}1' | sed "s/BARCODE/Barcode/g" | awk 'BEGIN{{FS=OFS="\\t"}} $3=="doublet" {{$2="doublet"}}1' | sed "s/0\t0/Assignment\tDropletType/g" | sed "1s/\t/\tdemuxalot_refined_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.demuxalot_temp}
        """


if os.path.exists(outdir + "/scrublet_pctl.tsv"):
    rule scrublet_results_temp:
        input:
            barcode=outdir + "/{pool}/matrix_out/barcodes.tsv"
        output:
            scrublet_temp=outdir + "/{pool}/CombinedResults/scrublet_temp.txt"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        params:
            q = "-q short.q",
            path=outdir + "/{pool}/scrublet_",
            pctl=lambda wildcards: scrublet_decisions.Percentile[scrublet_decisions.Pool == wildcards.pool].iloc[0]
        threads: 1
        shell:
            """
            awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {params.path}{params.pctl}/scrublet_results.txt > {output.scrublet_temp}
            """

rule scds_results_temp:
    input:
        outdir + "/{pool}/scds/scds_doublets.txt"
    output:
        outdir + "/{pool}/CombinedResults/scds_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q"
    threads: 1
    shell:
        """
        awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input}  > {output}
        """

rule DoubletFinder_results_temp:
    input:
        outdir + "/{pool}/DoubletFinder/DoubletFinder_doublets.txt"
    output:
        outdir + "/{pool}/CombinedResults/DoubletFinder_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q"
    threads: 1
    shell:
        """
        awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input}  > {output}
        """

rule DoubletDetection_results_temp:
    input:
        barcode=outdir + "/{pool}/matrix_out/barcodes.tsv",
        results=outdir + "/{pool}/DoubletDetection/DoubletDetection_predicted_doublet.txt"
    output:
        scrublet_temp=outdir + "/{pool}/CombinedResults/DoubletDetection_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    params:
        q = "-q short.q"
    threads: 1
    shell:
        """
        sed "s/0.0/singlet/g" {input.results} | sed "s/1.0/doublet/g" | paste {input.barcode} - | sed "1 i Barcode\tDoubletDetection_DropletType" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.scrublet_temp}
        """

if os.path.exists(outdir + "/DoubletDecon_rhops.tsv"):
    rule DoubletDecon_results_temp:
        input:
            barcode=outdir + "/{pool}/matrix_out/barcodes.tsv"
        output:
            temp=temp(outdir + "/{pool}/CombinedResults/DoubletDecon_temp_temp.txt"),
            final=outdir + "/{pool}/CombinedResults/DoubletDecon_temp.txt"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        threads: 1
        params:
            q = "-q short.q",
            path=outdir + "/{pool}/",
            rhop=lambda wildcards: DoubletDecon_decisions.rhop[DoubletDecon_decisions.Pool == wildcards.pool].iloc[0]
        shell:
            """
            tail -n +2 {params.path}DoubletDecon_rhop{params.rhop}/Final_doublets_groups_DoubletDecon_results.txt | awk '{{print $1}}' | sed 's/\"//g' | sed "s/$/\tdoublet/" | tr "." "-" | sed "1 i Barcode\tDoubletDecon_DropletType" > {output.temp}
            tail -n +2 {params.path}DoubletDecon_rhop{params.rhop}/Final_nondoublets_groups_DoubletDecon_results.txt | awk '{{print $1}}' | sed 's/\"//g' | sed "s/$/\tsinglet/" | tr "." "-" >> {output.temp}
            awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {output.temp} > {output.final}
            """

rule solo_results_temp:
    input:
        outdir + "/{pool}/solo/solo_results.txt"
    output:
        outdir + "/{pool}/CombinedResults/solo_results_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
    shell:
        """
        awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input} > {output}
        """ 



rule scDblFinder_results_temp:
    input:
        outdir + "/{pool}/scDblFinder/scDblFinder_results.txt"
    output:
        outdir + "/{pool}/CombinedResults/scDblFinder_results_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
    shell:
        """
        awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input} > {output}
        """ 

rule scDblFinder_known_doublets_results_temp:
    input:
        ancient(outdir + "/{pool}/scDblFinder_known_doublets/scDblFinder_results.txt")
    output:
        outdir + "/{pool}/CombinedResults/scDblFinder_known_doublets_results_temp.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
    shell:
        """
        sed 's/scDblFinder/scDblFinder_known_doublets/g' {input} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output}
        """ 


if os.path.exists(outdir + "/scrublet_pctl.tsv"):
    if os.path.exists(outdir + "/DoubletDecon_rhops.tsv"):
        if os.path.exists(outdir + "/DoubletDetection_PASS_FAIL.tsv"):
            print("running because exist")
            rule join_results:
                input:
                    demultiplexing_combined = outdir + "/{pool}/CombinedResults/CombinedDemuxletResults.tsv",
                    demuxlet=outdir + "/{pool}/CombinedResults/demuxlet_temp.txt",
                    freemuxlet=outdir + "/{pool}/CombinedResults/freemuxlet_temp.txt",
                    scSplit=ancient(outdir + "/{pool}/CombinedResults/scSplit_temp.txt"),
                    souporcell=ancient(outdir + "/{pool}/CombinedResults/souporcell_temp.txt"),
                    vireo=outdir + "/{pool}/CombinedResults/vireo_temp.txt",
                    dropulation = outdir + "/{pool}/CombinedResults/dropulation_temp.txt",
                    demuxalot = outdir + "/{pool}/CombinedResults/demuxalot_temp.txt",
                    demuxalot_refined = outdir + "/{pool}/CombinedResults/demuxalot_refined_temp.txt",
                    scrublet=outdir + "/{pool}/CombinedResults/scrublet_temp.txt",
                    scds=outdir + "/{pool}/CombinedResults/scds_temp.txt",
                    DoubletDetection= ancient(outdir + "/{pool}/CombinedResults/DoubletDetection_temp.txt"),
                    DoubletFinder=outdir + "/{pool}/CombinedResults/DoubletFinder_temp.txt",
                    DoubletDecon=outdir + "/{pool}/CombinedResults/DoubletDecon_temp.txt",
                    solo = outdir + "/{pool}/CombinedResults/solo_results_temp.txt",
                    scDblFinder = outdir + "/{pool}/CombinedResults/scDblFinder_results_temp.txt",
                    scDblFinder_known_doublets = outdir + "/{pool}/CombinedResults/scDblFinder_known_doublets_results_temp.txt"
                output:
                    outdir + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
                resources:
                    mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
                    disk_per_thread_gb=lambda wildcards, attempt: attempt * 5,
                params:
                    q = "-q short.q"
                threads: 1
                shell:
                    """
                    join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,2.2" {input.demultiplexing_combined} {input.scrublet} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
                    join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,2.2,2.3" - {input.scds} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
                    join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,2.2" - {input.DoubletDetection} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
                    join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,2.2,2.3" - {input.DoubletFinder} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
                    join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,2.2" - {input.DoubletDecon} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
                    join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,2.2,2.3" - {input.solo} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
                    join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,2.2,2.3" - {input.scDblFinder} | sed "s/ /\t/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' | \
                    join -a1 -a2 -1 1 -2 1 -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,2.2,2.3" - {input.scDblFinder_known_doublets} | sed "s/ /\t/g" > {output}
                    """


            rule singlet_doublets:
                input:
                    key = ancient( outdir + "/SimulatedOverlap/{pool}/PoolKeys.rds"),
                    assign = outdir + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv",
                    meta = sample_file
                output:
                    variables = outdir + "/SimulatedOverlap/{pool}/sing_dbl_variables.tsv",
                    result = outdir + "/SimulatedOverlap/{pool}/proportion_most_common_cell_type_per_barcode.rds"
                resources:
                    mem_per_thread_gb=lambda wildcards, attempt: attempt * 360,
                    disk_per_thread_gb=lambda wildcards, attempt: attempt * 360,
                threads: 1
                params:
                    q = "-q short.q",
                    sif = demuxafy_sif,
                    directory = outdir,
                    outdir = outdir + "/SimulatedOverlap/{pool}/",
                    script = simulation_scripts_dir + "/PullSingletsDoublets_separate_pools.R"
                conda:
                    simulation_scripts_dir + "/../../../conda_environments/generalR.yaml"
                shell:
                    """
                    echo {wildcards.pool} > {output.variables}
                    echo {resources.mem_per_thread_gb} >> {output.variables}
                    echo {threads} >> {output.variables}
                    echo {params.directory} >> {output.variables}
                    echo {params.outdir} >> {output.variables}
                    echo {input.meta} >> {output.variables}
                    Rscript {params.script} {output.variables}
                    """


            rule singlet_doublets_individual:
                input:
                    key = ancient( outdir + "/SimulatedOverlap/{pool}/PoolKeys.rds"),
                    assign = outdir + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv",
                    meta = sample_file,
                output:
                    variables = outdir + "/SimulatedOverlap/{pool}/sing_dbl_variables_single_soft.tsv",
                    common = outdir + "/SimulatedOverlap/{pool}/singlet_counts_per_barcode_single_soft.rds",
                    result = outdir + "/SimulatedOverlap/{pool}/proportion_most_common_cell_type_per_barcode_single_soft.rds"
                resources:
                    mem_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                    disk_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                threads: 1
                params:
                    q = "-q short.q",
                    sif = demuxafy_sif,
                    directory = outdir,
                    outdir = outdir + "/SimulatedOverlap/{pool}/",
                    script = simulation_scripts_dir + "/PullSingletsDoublets_separate_pools_single_software.R"
                conda:
                    simulation_scripts_dir + "/../../../conda_environments/generalR.yaml"
                shell:
                    """
                    singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
                    singularity exec {params.sif} echo {resources.mem_per_thread_gb} >> {output.variables}
                    singularity exec {params.sif} echo {threads} >> {output.variables}
                    singularity exec {params.sif} echo {params.directory} >> {output.variables}
                    singularity exec {params.sif} echo {params.outdir} >> {output.variables}
                    singularity exec {params.sif} echo {input.meta} >> {output.variables}
                    Rscript {params.script} {output.variables}
                    """


            rule singlet_doublets_individual_multisoft:
                input:
                    key = ancient( outdir + "/SimulatedOverlap/{pool}/PoolKeys.rds"),
                    assign = outdir + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv",
                    meta = sample_file,
                output:
                    variables = outdir + "/SimulatedOverlap/{pool}/sing_dbl_variables_multi_soft.tsv",
                    common = outdir + "/SimulatedOverlap/{pool}/singlet_counts_per_barcode_w_ref_data.rds",
                    result = outdir + "/SimulatedOverlap/{pool}/proportion_most_common_individual_per_barcode_w_ref_data.rds"
                resources:
                    mem_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                    disk_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                threads: 1
                params:
                    script = simulation_scripts_dir + "/PullSingletsDoublets_separate_pools.R",
                    q = "-q short.q",
                    directory = outdir,
                    outdir = outdir + "/SimulatedOverlap/{pool}/",
                conda:
                    simulation_scripts_dir + "/../../../conda_environments/generalR.yaml"
                shell:
                    """
                    echo {wildcards.pool} > {output.variables}
                    echo {resources.mem_per_thread_gb} >> {output.variables}
                    echo {threads} >> {output.variables}
                    echo {params.directory} >> {output.variables}
                    echo {params.outdir} >> {output.variables}
                    echo {input.meta} >> {output.variables}
                    Rscript {params.script} {output.variables}
                    """


            rule best_single_soft:
                input:
                    common = ancient(outdir + "/SimulatedOverlap/{pool}/singlet_counts_per_barcode_single_soft.rds"),
                    barcodes = ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
                    meta = sample_file,
                output:
                    variables = outdir + "/SimulatedOverlap/{pool}/single_software_metrics.tsv",
                resources:
                    mem_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                    disk_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                threads: 1
                params:
                    script = simulation_scripts_dir + "/Best_Single_Softwares_per_pool.R",
                    q = "-q short.q",
                    outdir = outdir,
                conda:
                    simulation_scripts_dir + "/../../../conda_environments/generalR.yaml"
                shell:
                    """
                    Rscript {params.script} -o {params.outdir} -p {wildcards.pool}
                    """


            rule best_doublet_detecting:
                input:
                    common = ancient(outdir + "/SimulatedOverlap/{pool}/singlet_counts_per_barcode_single_soft.rds"),
                    annotation = ancient(outdir + "/SimulatedOverlap/DropletAnnotation/simulated_droplet_types.rds"),
                    barcodes = ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
                    meta = sample_file
                output:
                    variables = outdir + "/SimulatedOverlap/{pool}/doublet_detecting_metrics.tsv",
                    classifications = outdir + "/SimulatedOverlap/{pool}/doublet_detecting_combos_classifications.rds"
                resources:
                    mem_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                    disk_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                threads: 1
                params:
                    q = "-q short.q",
                    outdir = outdir,
                    script = simulation_scripts_dir + "/Best_DoubletDetecting_Combination_per_pool.R"
                conda:
                    simulation_scripts_dir + "/../../../conda_environments/generalR.yaml"
                shell:
                    """
                    Rscript {params.script} -o {params.outdir} -p {wildcards.pool}
                    """



            rule best_demultiplexing:
                input:
                    common = ancient(outdir + "/SimulatedOverlap/{pool}/singlet_counts_per_barcode_single_soft.rds"),
                    annotation = ancient(outdir + "/SimulatedOverlap/DropletAnnotation/simulated_droplet_types.rds"),
                    barcodes = ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
                    meta = sample_file
                output:
                    variables = outdir + "/SimulatedOverlap/{pool}/demultiplexing_metrics.tsv",
                resources:
                    mem_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                    disk_per_thread_gb=lambda wildcards, attempt: attempt * 12,
                threads: 1
                params:
                    q = "-q short.q",
                    outdir = outdir,
                    script = simulation_scripts_dir + "/Best_Demultiplexing_Combination_per_pool.R"
                shell:
                    """
                    Rscript {params.script} -o {params.outdir} --pool {wildcards.pool}
                    """


            rule combination:
                input:
                    common = ancient(outdir + "/SimulatedOverlap/{pool}/singlet_counts_per_barcode_single_soft.rds"),
                    annotation = ancient(outdir + "/SimulatedOverlap/DropletAnnotation/simulated_droplet_types.rds"),
                    barcodes = ancient(outdir + "/{pool}/matrix_out/barcodes.tsv.gz"),
                    meta = sample_file
                output:
                    variables = outdir + "/SimulatedOverlap/{pool}/demultiplexing_doublet_detecting_metrics_w_single_softs.tsv",
                    metrics = outdir + "/SimulatedOverlap/{pool}/demultiplexing_doublet_detecting_metrics.tsv"
                resources:
                    mem_per_thread_gb=lambda wildcards, attempt: attempt * 128,
                    disk_per_thread_gb=lambda wildcards, attempt: attempt * 128,
                threads: 1
                params:
                    q = "-q short.q",
                    outdir = outdir,
                    script = simulation_scripts_dir + "/Best_Combination_per_pool.R"
                shell:
                    """
                    Rscript {params.script} -o {params.outdir} --pool {wildcards.pool}
                    """


#####################################
############ SUBSET VCFS ############
#####################################

### Subset vcfs based on individuals in each pool and locations in cluster-level vcfs
rule freemuxlet_pool_vcf:
    input:
        genotypes=SNP_GENOTYPES + ".gz",
        cluster_geno=outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.vcf.gz"
    output:
        vcf = outdir + "/{pool}/popscle/freemuxlet/Individual_genotypes_subset.vcf.gz",
        qstat = outdir + "/benchmarks/{pool}_freemuxlet_pool_vcf.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 5,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0],
    benchmark:
        outdir + "benchmarks/{pool}_freemuxlet_pool_vcf.benchmark"
    shell:
        """
        singularity exec {params.sif} bcftools view -s {params.individuals} -R {input.cluster_geno} -Oz -o {output.vcf} {input.genotypes}
        qstat -j $JOB_ID > {output.qstat}
        """

rule scSplit_pool_vcf:
    input:
        genotypes=SNP_GENOTYPES + ".gz",
        cluster_geno=ancient(outdir + "/{pool}/scSplit/scSplit.vcf"),
    output:
        vcf = outdir + "/{pool}/scSplit/Individual_genotypes_subset.vcf.gz",
        qstat = outdir + "/benchmarks/{pool}_scSplit_pool_vcf.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 5,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "benchmarks/{pool}_scSplit_pool_vcf.benchmark"
    shell:
        """
        singularity exec {params.sif} bcftools view -s {params.individuals} -R {input.cluster_geno} -Oz -o {output.vcf} {input.genotypes}
        qstat -j $JOB_ID > {output.qstat}
        """

rule souporcell_pool_vcf:
    input:
        genotypes=SNP_GENOTYPES + ".gz",
        cluster_geno=outdir + "/{pool}/souporcell/cluster_genotypes.vcf",
    output:
        vcf = outdir + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        qstat = outdir + "/benchmarks/{pool}_souporcell_pool_vcf.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 5,
    threads: 1
    params:
        q = "-q short.q",
        sif=demuxafy_sif,
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    benchmark:
        outdir + "benchmarks/{pool}_souporcell_pool_vcf.benchmark"
    shell:
        """
        singularity exec {params.sif} bcftools view -s {params.individuals} -R {input.cluster_geno} -Oz -o {output.vcf} {input.genotypes}
        qstat -j $JOB_ID > {output.qstat}
        """

rule individual_keys:
    input:
        freemuxlet = outdir + "/{pool}/popscle/freemuxlet/Individual_genotypes_subset.vcf.gz",
        scSplit = ancient(outdir + "/{pool}/scSplit/Individual_genotypes_subset.vcf.gz"),
        souporcell = ancient(outdir + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz"),
        meta = sample_file
    output:
        variables = outdir + "/SimulatedOverlap/{pool}/key_variables.tsv",
        key = outdir + "/SimulatedOverlap/{pool}/PoolKeys.rds",
        qstat = outdir + "/benchmarks/{pool}_individual_keys.qstat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 4,
    threads: 2
    params:
        q = "-q short.q",
        sif = demuxafy_sif,
        directory = outdir,
        outdir = outdir + "/SimulatedOverlap/{pool}/",
        script = simulation_scripts_dir + "/Assign_Indiv_by_Geno_separate_pools.R"
    benchmark:
        outdir + "benchmarks/{pool}_individual_keys.benchmark"
    shell:
        """
        singularity exec {params.sif} echo {input.meta} > {output.variables}
        singularity exec {params.sif} echo {params.directory} >> {output.variables}
        singularity exec {params.sif} echo {params.outdir} >> {output.variables}
        singularity exec {params.sif} echo {wildcards.pool} >> {output.variables}
        singularity exec {params.sif} Rscript {params.script} {output.variables}

        qstat -j $JOB_ID > {output.qstat}
        """