print("loading libraries")

.libPaths("/usr/local/lib/R/site-library")
library("tidyr")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("ComplexHeatmap")

library(circlize)
library(ggpubr)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(knitr)
library(future.apply)
library(data.table)


set.seed(79)


### Bring in arguments
args <- commandArgs(trailingOnly = TRUE)
arguments <- read.table(args, header = F)
pool <- as.character(arguments[1,])
mem <- as.numeric(as.character(arguments[2,]))
threads <- as.numeric(as.character(arguments[3,]))
dir <- as.character(arguments[4,])
out <- as.character(arguments[5,])
meta_file <- as.character(arguments[6,])
print(meta_file)

dir <- paste0(dir, "/")
out <- paste0(out, "/")

options(future.globals.maxSize= 1024^2*as.numeric(mem))
plan(multicore(workers = as.numeric(threads))) ## causes errors with 

##### Set directories #####
sim_pool_meta <- read_delim(meta_file, delim = "\t")
pool <- pool
print(pool)


##### Set up universal 
demultiplexing_list <- c("demuxalot", "demuxalot_refined","demuxlet", "dropulation", "freemuxlet","scSplit","souporcell","vireo")
names(demultiplexing_list) <- c("Demuxalot", "Demuxalot (refined)","Demuxlet", "Dropulation", "Freemuxlet","ScSplit","Souporcell","Vireo")
doublet_detection_list <- c("DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder", "scDblFinder_known_doublets", "scds",  "scrublet", "solo")
names(doublet_detection_list) <- c("DoubletDecon", "DoubletDetection","DoubletFinder",  "ScDblFinder", "ScDblFinder (known doublets)", "Scds", "Scrublet", "Solo")
softwares <- c(demultiplexing_list, doublet_detection_list)


##### Read in data
result_files <- paste0(dir,pool,"/CombinedResults/CombinedDropletAssignments.tsv")


##### Read in the Pool Metadata #####
sim_pool_meta <- sim_pool_meta[,c("Individuals","Pool" )]
sim_pool_indivs <- sim_pool_meta[which(sim_pool_meta$Pool == pool), "Individuals"]

sim_pool_indivs <- separate(sim_pool_indivs, col = Individuals, into = paste0("Individual_",1:(str_count(sim_pool_indivs$Individuals[1], pattern = ",") + 1)), sep = ",") 

sim_pool_indivs <- pivot_longer(sim_pool_indivs, cols = colnames(sim_pool_indivs), names_to = "Individual_Number", values_to = "Individual")

sim_pool_indivs$Individual_Number <- as.numeric(as.character(gsub("Individual_","",sim_pool_indivs$Individual_Number)))


##### Read in the cell_info barcode files #####
simulated_barcodes <- fread(paste0(dir,pool,"/matrix_out/barcodes.tsv.gz"), sep = "\t", header = FALSE, col.names = c("Barcode"))

simulated_barcodes$DropletType <- ifelse(grepl(":", simulated_barcodes$Barcode), "doublet","singlet")
simulated_barcodes$DoubletType_DoubletType <- ifelse(gsub("[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", gsub(":.+", "", simulated_barcodes$Barcode)) == gsub(".+:[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", simulated_barcodes$Barcode), "homogenic_doublet", ifelse(grepl(":", simulated_barcodes$Barcode), "heterogenic_doublet","singlet"))


##### Join barcodes with individual identifiers #####
# simulated_barcodes_ref <- left_join(simulated_barcodes, sim_pool_indivs, by = c("Individual_Number"))
simulated_barcodes$Individual <- gsub("[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", gsub(":.+", "", simulated_barcodes$Barcode))
simulated_barcodes$Individual2 <- ifelse(!grepl(":", simulated_barcodes$Barcode), NA, gsub(".+:[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", simulated_barcodes$Barcode))

simulated_barcodes$Individual <- ifelse(simulated_barcodes$DropletType == "doublet", "doublet", simulated_barcodes$Individual)
simulated_barcodes <- simulated_barcodes[match(unique(simulated_barcodes$Barcode), simulated_barcodes$Barcode),]


##### Read in Files #####
results_list <- read_delim(result_files, delim = "\t")

##### Read in the key tables for freemuxlet, scSplit and souporcell #####
key <- readRDS(paste0(out,"PoolKeys.rds"))
key$Cluster_ID <- gsub("CLUST","",key$Cluster_ID)
key$Software <- paste0(key$Software, "_Assignment")


##### Pivot the dataframes longer for joining#####
results_list$scSplit_Assignment <- as.character(results_list$scSplit_Assignment)

results_list_long <- pivot_longer(results_list,cols = paste0(demultiplexing_list, "_Assignment"), names_to = "Software")

##### Left_join the common assignments to the dataframe #####
results_list_long <- left_join(results_list_long,key, by = c("value" = "Cluster_ID", "Software" = "Software"))
results_list_long$Genotype_ID <- ifelse(results_list_long$Software == "demuxalot_Assignment", results_list_long$value, results_list_long$Genotype_ID)
results_list_long$Genotype_ID <- ifelse(results_list_long$Software == "demuxalot_refined_Assignment", results_list_long$value, results_list_long$Genotype_ID)
results_list_long$Genotype_ID <- ifelse(results_list_long$Software == "demuxlet_Assignment", results_list_long$value, results_list_long$Genotype_ID)
results_list_long$Genotype_ID <- ifelse(results_list_long$Software == "dropulation_Assignment", results_list_long$value, results_list_long$Genotype_ID)
results_list_long$Genotype_ID <- ifelse(results_list_long$Software == "vireo_Assignment", results_list_long$value, results_list_long$Genotype_ID)
results_list_long$Genotype_ID <- ifelse(results_list_long$value == "doublet", "doublet", results_list_long$Genotype_ID)
results_list_long$Genotype_ID <- ifelse(is.na(results_list_long$Genotype_ID), "unassigned", results_list_long$Genotype_ID)

rm(key)


##### Pivot wider to get the counts of singlet, doublets... and intersection across softwares #####
results_list_wide <- results_list_long
results_list_wide$value <- NULL
results_list_wide$Correlation <- NULL
results_list_wide <- pivot_wider(results_list_wide,names_from = "Software", values_from = "Genotype_ID")

rm(results_list_long)

##### Add in DoubletDecon_DropletType for those where missing #####
if(!any(grep("DoubletDecon_DropletType", colnames(results_list_wide)))){
  results_list_wide$DoubletDecon_DropletType <- "unassigned"
}

for (method in softwares){
  cols <- paste0(method,"_DropletType")
  results_list_wide[,cols] <- replace_na(unlist(results_list_wide[,cols]), "unassigned")
}


##### Order the ref barcodes to be in same order as results
simulated_barcodes_ordered <- simulated_barcodes[match(results_list_wide$Barcode,simulated_barcodes$Barcode),]



message("Creating List of Softwares")
##### make a list of all softwares #####
ones <- softwares

all_combn <- c(ones)

all_combn_filt <- c()
for (soft in demultiplexing_list){
  all_combn_filt <- unique(c(all_combn_filt, all_combn[grep(soft, all_combn)]))
}


##### Make a loop to make a dataframe to provide assignments for intersection of doublets #####
singlets <- results_list_wide[,c("Barcode")]

percent_correct <- results_list_wide[,c("Barcode")]

individual <- results_list_wide[,c("Barcode")]

individual_proportions <- results_list_wide[,c("Barcode")]

cell_type <- results_list_wide[,c("Barcode")]

cell_type_proportion <- results_list_wide[,c("Barcode")]

names <- names(individual_proportions)

message("Creating Intersection and Union Assignments for each combination")


for (comparison in 1:length(ones)){
  print(comparison)
  combination <- ones[comparison]
  softwares <- strsplit(ones[comparison],"-")[[1]]

  ### Add the number of singlets for the combination in question ###
  singlets[,combination] <- rowSums(results_list_wide[,paste0(softwares,"_DropletType")] == "singlet", na.rm = TRUE)

  if (any(softwares %in% demultiplexing_list)){
	### Add the most common ID assignment
	individual[,combination] <- future_apply(results_list_wide[,paste0(softwares[(softwares %in% demultiplexing_list)], "_Assignment")], 1, function(y) names(which.max(table(y))))


	### Add the proportion of the most common ID assignment
	individual_proportions[,combination] <- rowSums(results_list_wide[,paste0(softwares[(softwares %in% demultiplexing_list)], "_Assignment")] == separate(individual[,combination], sep = "-", into = "assign", col = combination)$assign, na.rm = TRUE)/length(softwares[(softwares %in% demultiplexing_list)])
  }


  ### Add the percent correct for each combination
  if (any(softwares %in% doublet_detection_list)){
    percent_correct[,combination] <- (rowSums(results_list_wide[,paste0(softwares[(softwares %in% doublet_detection_list)],"_DropletType")] == simulated_barcodes_ordered$DropletType, na.rm = TRUE))/length(softwares)
  } else {
    percent_correct[,combination] <- rowSums(results_list_wide[,paste0(softwares[(softwares %in% demultiplexing_list)],"_DropletType")] == simulated_barcodes_ordered$DropletType & results_list_wide[,paste0(softwares[(softwares %in% demultiplexing_list)],"_Assignment")] == simulated_barcodes_ordered$Individual, na.rm = TRUE)/length(softwares)
  }


    ### Add the most common ID assignment
  cell_type[,combination] <- ifelse(singlets[,combination] >= length(softwares)/2, "singlet", "doublet")



    ### Add the proportion of the most common ID assignment or most common cell type
  if (any(softwares %in% doublet_detection_list)){
    cell_type_proportion[,combination] <- rowSums(results_list_wide[,paste0(softwares[(softwares %in% doublet_detection_list)], "_DropletType")] == c(cell_type[,combination]), na.rm = TRUE)
  } else {
    cell_type_proportion[,combination] <- individual_proportions[,combination][,1]
  }
}


saveRDS(singlets, paste0(out, "/singlet_counts_per_barcode_single_soft.rds"))
saveRDS(individual, paste0(out, "/most_common_individual_per_barcode_single_soft.rds"))
saveRDS(individual_proportions, paste0(out, "/proportion_most_common_individual_per_barcode_single_soft.rds"))
saveRDS(percent_correct, paste0(out, "/percent_correct_per_barcode_single_soft.rds"))
saveRDS(cell_type, paste0(out, "/proportion_most_common_cell_type_barcode_single_soft.rds"))
saveRDS(cell_type_proportion, paste0(out, "/proportion_most_common_cell_type_per_barcode_single_soft.rds"))



##### Left join to add in the reference individual assignments #####
message("Add reference information to dataframes")
singlets_ref <- left_join(singlets, simulated_barcodes, by = c("Barcode"))
singlets_ref$Individual <- ifelse(singlets_ref$DropletType == "singlet", singlets_ref$Individual, singlets_ref$DropletType)


individual_ref <- left_join(individual, simulated_barcodes, by = c("Barcode"))
individual_ref$Individual <- ifelse(individual_ref$DropletType == "singlet", individual_ref$Individual, individual_ref$DropletType)


individual_proportions_ref <- left_join(individual_proportions, simulated_barcodes, by = c("Barcode"))
individual_proportions_ref$Individual <- ifelse(individual_proportions_ref$DropletType == "singlet", individual_proportions_ref$Individual, individual_proportions_ref$DropletType)


saveRDS(singlets_ref, paste0(out, "/singlet_counts_per_barcode_w_ref_data_single_soft.rds"))
saveRDS(individual_ref, paste0(out, "/most_common_individual_per_barcode_w_ref_data_single_soft.rds"))
saveRDS(individual_proportions_ref, paste0(out, "/proportion_most_common_individual_per_barcode_w_ref_data_single_soft.rds"))



#### Next steps:
##### 1. Make heatmap of singlets results
##### 1.b. Make heatmap of singlets results - pull out single softwares
##### 2. Make PCA and/or UMAP of results
