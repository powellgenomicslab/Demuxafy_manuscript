.libPaths("/usr/local/lib/R/site-library")
library(scDblFinder)
# library(DropletUtils)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)


args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
out <- arguments[1,]
tenX <- arguments[2,]
doublet_ratio <- as.numeric(as.character(arguments[3,]))
dbl_file <- arguments[4,]

print(tenX)
### read in the known doublets if running with known doublets
# if (nrow(arguments > 2)){
if (!is.na(dbl_file)){
    known_doublets <- read_delim(dbl_file, col_names = c("Doublets"), delim = "\t")
} else {
    known_doublets <- NA
}
print(known_doublets)


### Read in data as an sce object ###
counts <- Read10X(as.character(tenX[1]), gene.column = 1)
sce <- SingleCellExperiment(list(counts=counts))


### Calculate doublet ratio ###
# doublet_ratio <- ncol(sce)/1000*0.008
print(doublet_ratio)

### Calculate Singlets and Doublets ###
if (all(is.na(known_doublets))){

    print("Running without known doublets")
    sce <- scDblFinder(sce, dbr=doublet_ratio)

    ### Make a dataframe of the results ###
    results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score)

} else if (all(known_doublets$Doublets == FALSE)) {

    print("No common doublets across all softwares")
    results <- data.frame(matrix(nrow = 0, ncol = 3))
    colnames(results) <- c("Barcode", "scDblFinder_DropletType", "scDblFinder_Score")

} else {

    print("Running with known doublets")
    sce <- scDblFinder(sce, dbr=doublet_ratio, knownDoublets = known_doublets$Doublets)

    ### Make a dataframe of the results ###
    results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score)
}


write_delim(results, path = paste0(out,"scDblFinder_results.txt"), delim = "\t")