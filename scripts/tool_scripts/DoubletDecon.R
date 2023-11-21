.libPaths("/usr/local/lib/R/site-library")
library(DoubletDecon)
library(tidyverse)
library(Seurat)
library(ggplot2)

args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
pool <- arguments[1,]
out <- arguments[2,]
tenX <- as.character(arguments[3,])
rhop <- arguments[4,]
rhop <- as.numeric(as.character(rhop[1]))
res <- arguments[5,]
res <- as.numeric(as.character(res[1]))
message(paste0("The 10x file directory is:", tenX))

if (file.exists(paste0(out, "seurat_",res,".rds"))){
  seurat <- readRDS(paste0(out, "seurat_",res,".rds"))
} else {
  ## Read in data
  counts <- Read10X(as.character(tenX[1]), gene.column = 1)

  ## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
  seurat <- CreateSeuratObject(counts, min.features = 150)
  seurat <- SCTransform(seurat)
  seurat <- RunPCA(seurat)
  seurat <- RunUMAP(seurat, dims = 1:30)
  seurat <- FindNeighbors(seurat, dims = 1:30, verbose = FALSE)
  seurat <- FindClusters(seurat, verbose = FALSE, resolution = res)
  UMAP <- DimPlot(seurat, label = TRUE) + NoLegend()
  ggsave(paste0(out, "UMAP.png"), plot = UMAP)
  saveRDS(seurat, paste0(out, "seurat_",res,".rds"))
}

## Preprocess ##
processed <- Improved_Seurat_Pre_Process(seurat, num_genes=50, write_files=FALSE)

## Run Doublet Decon ##
results <- Main_Doublet_Decon(rawDataFile = processed$newExpressionFile, 
  groupsFile = processed$newGroupsFile, 
  filename = "DoubletDecon_results",
  location = out,
  fullDataFile = NULL, 
  removeCC = FALSE, 
  species = "hsa", 
  rhop = rhop,
  write = TRUE, 
  PMF = TRUE, 
  useFull = FALSE, 
  heatmap = FALSE, 
  centroids=TRUE, 
  num_doubs=100, 
  only50=FALSE, 
  min_uniq=4, 
  nCores=16)

