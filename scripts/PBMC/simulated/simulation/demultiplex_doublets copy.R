#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")
library("tidyverse")
library(data.table)


##### Set directories #####
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
dir <- arguments[1,]
pool <- arguments[2,]
outdir <- arguments[3,]

##### Read in Files #####
results <- read_delim(paste0(dir,"/",pool,"/CombinedResults/CombinedDemuxletResults.tsv"), delim = "\t")


##### Demuxlet may have not run for some cells, so fill the entries #####
results[is.na(results$demuxlet_nSNP), "demuxlet_nSNP"] <- 0
results[is.na(results$demuxlet_DropletType), "demuxlet_DropletType"] <- "unassigned"
results[is.na(results$demuxlet_Assignment), "demuxlet_Assignment"] <- "unassigned"
results[is.na(results$scSplit_DropletType), "scSplit_DropletType"] <- "unassigned"
results[is.na(results$scSplit_Assignment), "scSplit_Assignment"] <- "unassigned"

##### make a list of all softwares #####
softwares <- c("demuxalot", "demuxalot_refined", "demuxlet", "dropulation", "freemuxlet", "scSplit", "souporcell", "vireo")
demultiplexing_combn <- t(combn(softwares, 8, simplify = TRUE)) %>% apply(. , 1 , paste , collapse = "_" )


##### Make a dataframe to provide assignments for intersection of doublets #####
intersection_doublet_demultiplex <- results[,c("Barcode")]


### Get a df of the droplet type ###
temp_DropletType <- results[,paste0(softwares,"_DropletType")]

# Account for the fact that some softwares didn't succeed and needed to use 'fake' files with just header
temp_DropletType <- data.table(Filter(function(x)(length(unique(x))>1), temp_DropletType))


# Create column called "DropletType_temp"
# In this column, if number of singlet assignments per cell are all singlets, then label as singlet
# Otherwise, label as doublets
intersection_doublet_demultiplex <- intersection_doublet_demultiplex %>% mutate(DropletType_temp = if_else(rowSums(temp_DropletType == "doublet") == ncol(temp_DropletType), "doublet", "singlet"))


# ##### Create dataframe to be used as input for solo for common doublets - requires True False in single column saved as tsv, no header
# intersection_doublet_demultiplex_solo <- data.frame("Doublet" = ifelse(intersection_doublet_demultiplex$DropletType_temp == "doublet", "True", "False"))
# write_delim(intersection_doublet_demultiplex_solo, paste0(outdir,"/solo_known_doublets/solo_known_doublets.tsv"), col_names = FALSE)


##### Create dataframe to be used as input for scDblFinder #####
intersection_doublet_demultiplex_scDblFinder <- data.frame("Doublet" = ifelse(intersection_doublet_demultiplex$DropletType_temp == "doublet", TRUE, FALSE))
write_delim(intersection_doublet_demultiplex_scDblFinder, paste0(outdir,"/scDblFinder_known_doublets/scDblFinder_known_doublets.tsv"), col_names = FALSE)
