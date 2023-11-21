#!/bin/bash
## This script is to rusn th ambient_cellbender.sh script using qsub to submit each pool as a separate job in parallel on an sge system

# Define path to the CellBender singularity image
cellbender_path="/path/to/Demuxafy_manuscript/files/cellbender.sif" ### singularity image for cellbender, can be found on zenodo
meta="/path/to/Demuxafy_manuscript/files/fibroblast/fibroblast_sample_meta.tsv" ### metadatafile that contains information for each of the fibroblast pools; provided on zenodo


for pool in `awk '{print $1}' $meta | tail -n +1`
do

	# Define directories
	script="/path/to/Demuxafy_manuscript/scripts/fibroblast/ambient_cellbender.sh"
	data_path="/path/to/10x/data/dir/$pool/outs"
	output_path="/path/to/output/fibroblast/$pool/ambient"
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