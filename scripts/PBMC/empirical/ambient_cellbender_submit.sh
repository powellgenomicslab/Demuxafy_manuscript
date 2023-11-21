#!/bin/bash

# Define path to the CellBender singularity image
cellbender_path="/path/to/Demuxafy_manuscript/files/cellbender.sif" ### singularity image for cellbender, can be found on zenodo
meta="/path/to/Demuxafy_manuscript/files/PBMC/PBMC_sample_meta.tsv" ### metadatafile that contains information for each of the fibroblast pools; provided on zenodo


for pool in `awk '{print $1}' $meta | tail -n +1`
do

	echo $pool

	# Define directories
	script="/path/to/Demuxafy_manuscript/scripts/fibroblast/ambient_cellbender.sh"
	data_path="/path/to/10x/data/dir/$pool/outs"
	output_path="/path/to/output/PBMC/$pool/ambient"
	logs=$output_path/logs
	mkdir -p $logs


	qsub -S /bin/bash \
		-cwd \
		-N cellbender \
		-l mem_requested=60G \
		-l tmp_requested=60G \
		-e $output_path/logs \
		-o $output_path/logs \
		-r yes \
		-j y \
		-V \
		-v pool=$pool,output_path=$output_path,data_path=$data_path,cellbender_path=$cellbender_path \
		-l nvgpu=3 \
		-C '' $script

done