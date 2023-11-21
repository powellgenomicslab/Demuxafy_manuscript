.libPaths("/usr/local/lib/R/site-library")
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(tidyr)
library(tidyverse)
library(DropletUtils)

args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
pool <- arguments[1,]
out <- arguments[2,]
QCout <- arguments[3,]
QCout <- as.character(QCout[1])
tenX <- arguments[4,]
dblN <- as.numeric(as.character(arguments[5,]))
filedir <- as.character(arguments[6,])
message(paste0("The 10x file directory is:", tenX))
RB_genes <- read.delim(file = paste0(filedir,"/RibosomalGeneList_GeneID_ENSG.txt", header = T))
MT_genes <- read.delim(file = paste0(filedir,"/MtGeneList_GeneID_ENSG.txt", header = T))


## Add max future globals size for large pools
options(future.globals.maxSize=(850*1024^2))

## Define Functions 
mad_function <- function(seurat, column, number_mad){
    mad <- mad(seurat@meta.data[,column])
    low <- median(seurat@meta.data[,column]) - number_mad*mad
    high <- median(seurat@meta.data[,column]) + number_mad*mad
    print("The mad is:")
    print(mad)
    print("The lower bound is:")
    print(low)
    print("The upper bound is:")
    print(high)
    seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low & seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
    return(seurat)
}

## Read in data
counts <- Read10X(as.character(tenX[1]), gene.column = 1)
rownames(counts) <- gsub("\\.\\d+","", rownames(counts))

if (!file.exists(paste0(out,"seurat_processed_QCed.rds"))){

    ## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
    seurat <- CreateSeuratObject(counts)



    ##### Get the mitochondiral and ribosomal percentage QC metrics for each cell #####
    if ((sum(which(MT_genes$GeneID %in% rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat)))) & (sum(which(MT_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
        mt_features <- MT_genes$GeneID[MT_genes$GeneID %in% rownames(seurat)]
        seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mt_features)
    } else if ((sum(which(MT_genes$ENSG %in% rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (sum(which(MT_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
        mt_features <- MT_genes$ENSG[MT_genes$ENSG %in% rownames(seurat)]
        seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mt_features)
    } else if ((length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat))))){
        seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = rownames(seurat)[grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))])
    } else {
        message("Either you do not have mitochondrial genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
    }

    if ((sum(which(RB_genes$GeneID %in% rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat)))) & (sum(which(RB_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
        rb_features <- RB_genes$GeneID[RB_genes$GeneID %in% rownames(seurat)]
        seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rb_features)
    } else if ((sum(which(RB_genes$ENSG %in% rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (sum(which(RB_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
        rb_features <- RB_genes$ENSG[RB_genes$ENSG %in% rownames(seurat)]
        seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rb_features)
    } else if ((length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat))))){
        seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rownames(seurat)[grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))])
    } else {
        message("Either you do not have ribosomal genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
    }


    ## Calculate MADs
    seurat <- mad_function(seurat, "percent.mt",3)
    seurat <- mad_function(seurat, "percent.rb",3)
    seurat <- mad_function(seurat, "nCount_RNA", 3)
    seurat <- mad_function(seurat, "nFeature_RNA", 3)


    ##### Remove the outliers #####
    print(seurat)
    seurat <- subset(seurat, subset = percent.mt_mad == "NotOutlier") 
    seurat <- subset(seurat, subset = percent.rb_mad == "NotOutlier")
    seurat <- subset(seurat, subset = nCount_RNA_mad == "NotOutlier")
    seurat <- subset(seurat, subset = nFeature_RNA_mad == "NotOutlier")
    print(seurat)

    seurat <- SCTransform(seurat)
    seurat <- RunPCA(seurat)
    seurat <- RunUMAP(seurat, dims = 1:10)
    seurat <- FindNeighbors(seurat, dims = 1:10, verbose = FALSE)
    seurat <- FindClusters(seurat, verbose = FALSE, resolution = 0.3)
    umap <- DimPlot(seurat, label = TRUE) + NoLegend()
    ggsave(umap, filename = paste0(out,"umap.png"))


    write10xCounts(path = QCout, x = seurat[["RNA"]]@counts, version="2", type = "sparse", overwrite = TRUE)
    saveRDS(seurat,paste0(out,"seurat_processed_QCed.rds"))
} else {
    seurat <- readRDS(paste0(out,"seurat_processed_QCed.rds"))
}

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
plot <- ggplot(bcmvn, aes(pK, BCmetric)) +
    geom_point()
ggsave(plot, filename = paste0(out,"pKvBCmetric.png"))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Idents(seurat)
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- dblN
print(paste0("Expected number of doublets: ", dblN))
nExp_poi.adj <- round(dblN*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurat <- doubletFinder_v3(seurat, PCs = 1:10, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))])), nExp = dblN, reuse.pANN = FALSE, sct = TRUE)
doublets <- as.data.frame(cbind(colnames(seurat), seurat@meta.data[,grepl(paste0("pANN_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(seurat@meta.data))], seurat@meta.data[,grepl(paste0("DF.classifications_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(seurat@meta.data))]))
colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
doublets$DoubletFinder_DropletType <- gsub("Singlet","singlet",doublets$DoubletFinder_DropletType) %>% gsub("Doublet","doublet",.)

write_delim(doublets, path = paste0(out,"DoubletFinder_doublets.txt"), delim = "\t")
