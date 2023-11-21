#!/usr/bin/env Rscript
library(tidyverse)
library(circlize)
library(plyr)
library(scales)
library(ggforce)
library(data.table)
suppressMessages(suppressWarnings(library(argparse)))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-p", "--pool", required = FALSE, type = "character", help = "the pool ID")

args <- parser$parse_args()

##### Set up Functions #####
save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}


# ###### Set up Directories #####
dir <- paste0(args$out,"/SimulatedOverlap/")
dir.create(dir, recursive = TRUE)


##### Set up universal 
demultiplexing_list <- c("demuxalot", "demuxalot_refined","demuxlet", "dropulation", "freemuxlet","scSplit","souporcell","vireo")
names(demultiplexing_list) <- c("Demuxalot", "Demuxalot (refined)","Demuxlet", "Dropulation", "Freemuxlet","ScSplit","Souporcell","Vireo")
doublet_detection_list <- c("DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder", "scDblFinder_known_doublets", "scds",  "scrublet", "solo")
names(doublet_detection_list) <- c("DoubletDecon", "DoubletDetection","DoubletFinder",  "ScDblFinder", "ScDblFinder (known doublets)", "Scds", "Scrublet", "Solo")
softwares <- c(demultiplexing_list, doublet_detection_list)

softwares_dt <- data.table(software = softwares, Software = names(softwares))


software_colors <- c("#575E57", "#008751", "#5DBB4D", "#EFE305", "#F8BE33", "#F97919", "#D21F09", "#8D132B", "#F980BF", "#EAA1D1", "#BAAAFF", "#9B64E0", "#6B82E3", "#63D4FE", "#39BCBC", "#C4E4DF")
names(software_colors) <- names(softwares)


##### Read in Results ##### 
individual_assignment <- as.data.table(readRDS(paste0(dir, args$pool, "/most_common_individual_per_barcode_single_soft.rds")))

singlets <- readRDS(paste0(dir, args$pool, "/singlet_counts_per_barcode_single_soft.rds"))


##### read in barcodes 2 remove #####
barcodes2remove <- readRDS(paste0(dir,"DropletAnnotation/barcodes2remove.rds"))

print("updating names")
if (args$pool %in% names(barcodes2remove)){
	individual_assignment <- individual_assignment[!(Barcode %in% barcodes2remove[[args$pool]])]
	singlets <- singlets[!(singlets$Barcode %in% barcodes2remove[[args$pool]]),]
} else {
	base <- gsub("_ambient.+", "", args$pool) %>%
			gsub("_mt.+", "", .) %>%
			gsub("_subsampled", "", .) %>%
			gsub("_unevenN.+", "", .)
	individual_assignment <- individual_assignment[!(Barcode %in% barcodes2remove[[base]])]
	singlets <- singlets[!(singlets$Barcode %in% barcodes2remove[[base]]),]
}


singlets <- singlets %>% mutate_all(as.character)
singlets[singlets == "0"] <- "doublet"
singlets[singlets == "1"] <- "singlet"




########## Assess combinations of softwares ##########
##### Read in the Pool Metadata #####
sim_pool_meta <- read_delim(dir,"/updated_metadata_unique.tsv", delim = "\t")

sim_pool_meta <- sim_pool_meta[,c("Individuals","Pool" )]

sim_pool_indivs <- sim_pool_meta[which(sim_pool_meta$Pool == args$pool), "Individuals"]


sim_pool_indivs <- separate(sim_pool_indivs, col = Individuals, into = paste0("Individual_",1:(str_count(sim_pool_indivs$Individuals[1], pattern = ",") + 1)), sep = ",") 
sim_pool_indivs <- pivot_longer(sim_pool_indivs, cols = colnames(sim_pool_indivs), names_to = "Individual_Number", values_to = "Individual")
sim_pool_indivs$Individual_Number <- as.numeric(as.character(gsub("Individual_","",sim_pool_indivs$Individual_Number)))



##### Read in the cell_info barcode files #####
simulated_barcodes <- fread(paste0(dir,args$pool,"/matrix_out/barcodes.tsv.gz"), sep = "\t", header = FALSE, col.names = c("Barcode"))

simulated_barcodes$DropletType <- ifelse(grepl(":", simulated_barcodes$Barcode), "doublet","singlet")
simulated_barcodes$DoubletType_DoubletType <- ifelse(gsub("[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", gsub(":.+", "", simulated_barcodes$Barcode)) == gsub(".+:[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", simulated_barcodes$Barcode), "homogenic_doublet", ifelse(grepl(":", simulated_barcodes$Barcode), "heterogenic_doublet","singlet"))



##### Join barcodes with individual identifiers #####
simulated_barcodes$Individual <- gsub("[A,C,G,T,L,M,N,O,P,Q,R]+-", "", gsub(":.+", "", simulated_barcodes$Barcode))
simulated_barcodes$Individual2 <- ifelse(!grepl(":", simulated_barcodes$Barcode), NA, gsub(".+:[A,C,G,T,L,M,N,O,P,Q,R]+-", "", simulated_barcodes$Barcode))

simulated_barcodes$Individual <- ifelse(simulated_barcodes$DropletType == "doublet", "doublet", simulated_barcodes$Individual)
simulated_barcodes <- simulated_barcodes[match(unique(simulated_barcodes$Barcode), simulated_barcodes$Barcode),]

simulated_barcodes_ordered <- simulated_barcodes[match(individual_assignment$Barcode,simulated_barcodes$Barcode),]



demultiplexing_doublet_detection_list <- c(demultiplexing_list,doublet_detection_list)



##### Calculate necessary metrics for individual softwares ######
## True singlet (TP)
## False singlet (FP)
## True doublet (TN)
## False doublet (FN)
## Accuracy
## Balanced Accuracy
## Sensitivity (TPR)
## Specificity (TNR)
## Precision
## MCC

##### Set up dataframes #####
Metrics <- list()


Metrics <- data.frame(matrix(nrow = length(demultiplexing_doublet_detection_list), ncol = 13))
colnames(Metrics) <- c("software", "true_singlet", "false_singlet", "true_doublet", "false_doublet", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "npv", "ppv", "precision", "mcc")
Metrics$software <- demultiplexing_doublet_detection_list
rownames(Metrics) <- demultiplexing_doublet_detection_list
Metrics$pool <- args$pool


for (soft in demultiplexing_doublet_detection_list){
	## True singlet (TP)
	if (soft %in% demultiplexing_list){
		Metrics[soft,]$true_singlet <- as.numeric(length(which(singlets[, soft] == "singlet" & singlets[, soft] == simulated_barcodes_ordered$DropletType & individual_assignment[, ..soft] == simulated_barcodes_ordered$Individual)))
	} else {
		Metrics[soft,]$true_singlet <- as.numeric(length(which(singlets[, soft] == "singlet" & singlets[, soft] == simulated_barcodes_ordered$DropletType)))
	}
	# Metrics[soft,]$true_singlet[Metrics[soft,]$true_singlet == 0] <- NA

	## False singlet (FP)
	Metrics[soft,]$false_singlet <- as.numeric(length(which(singlets[, soft] == "singlet")) - Metrics[soft,]$true_singlet)
	# Metrics[soft,]$false_singlet[Metrics[soft,]$false_singlet == 0] <- NA

	### Update NA for pools where could not run (ie DoubletDecon)
	if (Metrics[soft,]$true_singlet == 0 & Metrics[soft,]$false_singlet == 0){
		Metrics[soft,]$true_doublet <- NA
		Metrics[soft,]$false_doublet <- NA
		Metrics[soft,]$accuracy
		Metrics[soft,]$sensitivity <- NA
		Metrics[soft,]$specificity <- NA
		Metrics[soft,]$npv <- NA
		Metrics[soft,]$ppv <- NA
		Metrics[soft,]$balanced_accuracy <- NA
		Metrics[soft,]$precision <- NA
		Metrics[soft,]$mcc <- NA
	} else{

		## True doublet (TN)
		Metrics[soft,]$true_doublet <- as.numeric(length(which(singlets[, soft] == "doublet" & singlets[, soft] == simulated_barcodes_ordered$DropletType)))


		## False doublet (FN)
		Metrics[soft,]$false_doublet <- as.numeric(length(which(singlets[, soft] == "doublet")) - Metrics[soft,]$true_doublet)


		## Accuracy
		Metrics[soft,]$accuracy <- as.numeric((Metrics[soft,]$true_singlet + Metrics[soft,]$true_doublet)/nrow(singlets))


		## Sensitivity (TPR)
		Metrics[soft,]$sensitivity <- as.numeric(Metrics[soft,]$true_singlet/length(which(simulated_barcodes_ordered$DropletType == "singlet")))


		## Specificity (TNR)
		Metrics[soft,]$specificity <- as.numeric(Metrics[soft,]$true_doublet/length(which(simulated_barcodes_ordered$DropletType == "doublet")))


		## Negative Predictive Value
		Metrics[soft,]$npv <- as.numeric(Metrics[soft,]$true_doublet/length(which(singlets[, soft] == "doublet")))


		## Positive Predictive Value
		Metrics[soft,]$ppv <- as.numeric(Metrics[soft,]$true_singlet/length(which(singlets[, soft] == "singlet")))


		## Balanced Accuracy
		Metrics[soft,]$balanced_accuracy <- as.numeric((Metrics[soft,]$sensitivity + Metrics[soft,]$specificity)/2)


		## Precision
		Metrics[soft,]$precision <- as.numeric(Metrics[soft,]$true_singlet/(Metrics[soft,]$true_singlet + Metrics[soft,]$false_singlet))

		
		## MCC
		Metrics[soft,]$mcc <- as.numeric(((Metrics[soft,]$true_singlet * Metrics[soft,]$true_doublet) - (Metrics[soft,]$false_singlet * Metrics[soft,]$false_doublet))/sqrt((Metrics[soft,]$true_singlet + Metrics[soft,]$false_singlet) * (Metrics[soft,]$true_singlet + Metrics[soft,]$false_doublet) * (Metrics[soft,]$true_doublet + Metrics[soft,]$false_singlet) * (Metrics[soft,]$true_doublet + Metrics[soft,]$false_doublet)))
	}
}


rownames(Metrics) <- NULL
Metrics <- as.data.table(Metrics)

Metrics$size <- factor(as.numeric(gsub("size","", Metrics$pool) %>% gsub("_.+", "", .)), levels = c(2,4,8,16,32,64,128))
Metrics$mt_percent <- ifelse(grepl("mt", Metrics$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_mt_","", Metrics$pool) %>% gsub("pctl", "", .)), 0)
Metrics$ambient_percent <- ifelse(grepl("ambient", Metrics$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_ambient_","", Metrics$pool) %>% gsub("pctl", "", .)), 0)
Metrics$subsampled <- ifelse(grepl("subsampled", Metrics$pool), "subsampled", "normal")
Metrics$uneven <- ifelse(grepl("unevenN", Metrics$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_unevenN_","", Metrics$pool) %>% gsub("pctl", "", .)), 0)


Metrics <- softwares_dt[Metrics, on = "software"]


fwrite(Metrics, paste0(dir, args$pool,"/single_software_metrics.tsv"), sep = "\t")
