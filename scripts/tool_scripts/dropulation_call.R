library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
singletFile <- args[1]
doubletFile <- args[2]
outMapFile <- args[3]



createDonorMap<-function (singletFile, doubletFile, outMapFile, fdrThreshold=0.05, doubletPvalueThreshold=0.9) {
    a=read.table(singletFile, header=T, stringsAsFactors = F, sep="\t")
    b=read.table(doubletFile, header=T, stringsAsFactors = F, sep="\t")
    
    confidentAssignmentCells=a[a$FDR_pvalue<=fdrThreshold,"cell"]
    singletCells=b[b$doublet_pval<doubletPvalueThreshold, "cell"]
    confidentAssignmentSingletCells=intersect(confidentAssignmentCells, singletCells)
    doublets = a$cell[!(a$cell %in% confidentAssignmentSingletCells)]

    mapa=a[match(confidentAssignmentSingletCells, a$cell), c("cell", "bestLikelihood", "bestSample")]
    colnames(mapa) = c("Barcode","dropulation_Likelihood", "dropulation_Assignment")
    mapa$dropulation_DropletType = "singlet"

    mapb= b[match(doublets, b$cell), c("cell", "mixedSampleLikelihood", "mixedSample")]
    colnames(mapb) = c("Barcode", "dropulation_Likelihood","dropulation_Assignment")
    mapb$dropulation_Assignment = "doublet"
    mapb$dropulation_DropletType = "doublet"

    map = rbind(mapa,mapb)

    map = left_join(map, a[,c("cell", "num_snps", "num_umis")], by = c("Barcode" = "cell"))
    colnames(map) = gsub("num_snps", "dropulation_Nsnps", colnames(map)) %>%
                        gsub("num_umis", "dropulation_Numis", .)


    fwrite(map, outMapFile, row.names=F, col.names = T, quote=F, sep="\t")
}


createDonorMap(singletFile, doubletFile, outMapFile)
