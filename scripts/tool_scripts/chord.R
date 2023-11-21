library(data.table)
library(adabag)
library(gbm)
library(Chord)
library(Seurat)



args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
pool <- arguments[1,]
out <- arguments[2,]
tenX <- arguments[3,]
dbl_rate <- as.numeric(as.character(arguments[4,]))

setwd(out)

## Read in data
counts <- Read10X(as.character(tenX[1]), gene.column = 1)
rownames(counts) <- gsub("\\.\\d+","", rownames(counts))

if (!file.exists(paste0(out,"seurat.rds"))){

    ## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
    seurat <- CreateSeuratObject(counts)

    saveRDS(seurat,paste0(out,"seurat.rds"))
} else {
    seurat <- readRDS(paste0(out,"seurat.rds"))
}

chord(seu=seurat,doubletrate=dbl_rate,overkill=T,outname="chord_results")