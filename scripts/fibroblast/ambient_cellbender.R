## this script takes the results from ambient_cellbender.sh and provides ambient RNA 

# Import libraries --------------------------------------------------------

library(tidyverse)
library(Seurat)
library(rhdf5)
library(data.table)




# Define directories ------------------------------------------------------
meta <- "/path/to/Demuxafy_manuscript/scripts/fibroblast/fibroblast_sample_meta.tsv"
dir <- "/path/to/output/" ### This is the path to the output used for cellbender
scdir <- "/path/to/cellranger/results/" # path to folder containing a folder for each pool that has the raw_feature_bc_matrix.h5 from 10x cellranger
outdir <- paste0(dir,"ambient/")
dir.create(outdir, recursive = TRUE)


## Set up variables ##
metadata <- fread(meta)
results <- metadata[, 'Pool']


# Process CellBender output -----------------------------------------------
# Get the path to the h5 file output by cellbender
results$cell_bender_results <- lapply(results$Pool, function(x){
	Read10X_h5(paste0(dir, x,"/ambient/output.h5"), use.names = TRUE)
})

# Import unfiltered counts (cellranger raw)
results$original_results <- lapply(results$Pool, function(x){
	print(x)
	Read10X_h5(paste0(scdir, x,"/outs/raw_gene_bc_matrices_h5.h5"), use.names = TRUE)
})


# Calculate the percentage of ambient RNA for each cell barcode 
results$cell_bender_results_filtered <- lapply(results$Pool, function(x){
	results[Pool == x]$cell_bender_results[[1]][,match(colnames(results[Pool == x]$original_results[[1]]), colnames(results[Pool == x]$cell_bender_results[[1]]))]
})

results$original_results_filtered <- lapply(results$Pool, function(x){
	results[Pool == x]$original_results[[1]][,match(colnames(results[Pool == x]$cell_bender_results[[1]]), colnames(results[Pool == x]$original_results[[1]]))]
})


# Add to the as metadata to the Seurat object
ambient_rna_list <- list()
ambient_rna_list <- lapply(results$Pool, function(x){
	ambient_rna_list[[x]] <- data.table(Barcode = colnames(results[Pool == x]$original_results_filtered[[1]]), percent_ambient_RNA = 100*(colSums(results[Pool == x]$original_results_filtered[[1]]) - colSums(results[Pool == x]$cell_bender_results_filtered[[1]]))/colSums(results[Pool == x]$original_results_filtered[[1]]), pool = x)
})
names(ambient_rna_list) <- results$Pool

### Read in cell barcodes ###
cell_barcodes <- lapply(results$Pool, function(pool){
	tmp <- fread(paste0(scdir,pool,"/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/barcodes.tsv"), sep = "\t",header = FALSE)
	setnames(tmp, "Barcode")
})
names(cell_barcodes) <- results$Pool


ambient_rna_list_cells <- lapply(names(ambient_rna_list) , function(x){
	ambient_rna_list[[x]][cell_barcodes[[x]], on= "Barcode"]
})
names(ambient_rna_list_cells) <- names(ambient_rna_list)

lapply(ambient_rna_list_cells, function(x) max(x$percent_ambient_RNA))
lapply(names(cell_barcodes), function(x) {nrow(cell_barcodes[[x]]) == nrow(ambient_rna_list_cells[[x]])})


saveRDS(ambient_rna_list_cells, paste0(outdir,"ambient_rna_per_barcode.rds"))



### Get percent of barcodes that are not in cellbender
cellbender_barcodes <- lapply(results$Pool, function(pool){
	tmp <- fread(paste0(dir,pool,"/ambient/output_cell_barcodes.csv"), sep = "\t",header = FALSE)
	setnames(tmp, "Barcode")
})
names(cellbender_barcodes) <- results$Pool


cellbender_comparison <- lapply(names(ambient_rna_list_cells), function(pool){
    bcs <- ambient_rna_list_cells[[pool]]$Barcode[!(cell_barcodes[[pool]]$Barcode %in% cellbender_barcodes[[pool]]$Barcode)]
    # print(ambient_rna_list_cells[[pool]][!(Barcode %in% cellbender_barcodes[[pool]]$Barcode),])
    # print(ambient_rna_list_cells[[pool]][which(ambient_rna_list_cells[[pool]]$Barcode %in% cellbender_barcodes[[pool]]$Barcode),])
    # subset(data.table(ambient_rna_list_cells[[pool]]), Barcode %in% cellbender_barcodes[[pool]]$Barcode)
    data.table(Percent_cellbender_cells = (nrow(ambient_rna_list_cells[[pool]]) - length(bcs))/nrow(ambient_rna_list_cells[[pool]]) * 100)
    # length(subset(ambient_rna_list_cells[[pool]], Barcode %in% bcs))
})

cellbender_comparison_dt <- do.call(rbind, cellbender_comparison)

mean(cellbender_comparison_dt$Percent_cellbender_cells)

