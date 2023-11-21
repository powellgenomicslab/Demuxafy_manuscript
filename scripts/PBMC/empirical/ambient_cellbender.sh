#!/bin/bash



# Define the number of expected cells and the number of droplets to analyse
cell_num=20000
droplets_num=50000

# Run CellBender's remove-background tool
singularity exec --nv $cellbender_path cellbender remove-background \
				--input ${data_path}/raw_gene_bc_matrices_h5.h5 \
				--output ${output_path}/output.h5 \
				--cuda \
				--expected-cells $cell_num \
				--total-droplets-included $droplets_num \
				--fpr 0.01 \
				--epochs 150
				# --posterior-batch-size 5 \
				# --cells-posterior-reg-calc

