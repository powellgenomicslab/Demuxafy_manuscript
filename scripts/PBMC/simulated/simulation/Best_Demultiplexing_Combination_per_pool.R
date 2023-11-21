library(dsLib)
library(tidyverse)
library(circlize)
library(plyr)
library(scales)
library(data.table)
library(future.apply)
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


###### Set up Directories #####
dir <- paste0(args$out,"/SimulatedOverlap/")
dir.create(dir, recursive = TRUE)
outdir <- paste0(dir, "/", args$pool)
droplet_type_dir <- paste0(dir,"/DropletAnnotation/")


##### Set up colors for plotting #####
downsample_colors <- c("#C7EAE5", "#35978F")
mt_percent_colors <- c("0" = "white", "5" = "#D0E0D3", "10" = "#A3C4A6", "25" = "#1E6D2D")
ambient_percent_colors <- c("0" = "white", "10" = "#FDEDD2", "20" = "#FBD491", "50" = "#E6AB00")
uneven_n_colors <- c("0" = "", "50" = "", "75" = "", "95" = "")


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

individual_assignment <- individual_assignment[!(Barcode %in% barcodes2remove[[args$pool]])]
singlets <- singlets[!(singlets$Barcode %in% barcodes2remove[[args$pool]]),]





droplet_type <- singlets %>% mutate_all(as.character)
droplet_type[droplet_type == "0"] <- "doublet"
droplet_type[droplet_type == "1"] <- "singlet"




##### Read in the cell_info barcode files #####
simulated_barcodes <- fread(paste0(args$out,"/",args$pool,"/matrix_out/barcodes.tsv.gz"), sep = "\t", header = FALSE, col.names = c("Barcode"))

simulated_barcodes$DropletType <- ifelse(grepl(":", simulated_barcodes$Barcode), "doublet","singlet")
simulated_barcodes$DoubletType_DoubletType <- ifelse(gsub("[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", gsub(":.+", "", simulated_barcodes$Barcode)) == gsub(".+:[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", simulated_barcodes$Barcode), "homogenic_doublet", ifelse(grepl(":", simulated_barcodes$Barcode), "heterogenic_doublet","singlet"))



##### Join barcodes with individual identifiers #####
simulated_barcodes$Individual <- gsub("[A,C,G,T,L,M,N,O,P,Q,R]+-", "", gsub(":.+", "", simulated_barcodes$Barcode))
simulated_barcodes$Individual2 <- ifelse(!grepl(":", simulated_barcodes$Barcode), NA, gsub(".+:[A,C,G,T,L,M,N,O,P,Q,R]+-", "", simulated_barcodes$Barcode))

simulated_barcodes$Individual <- ifelse(simulated_barcodes$DropletType == "doublet", "doublet", simulated_barcodes$Individual)
simulated_barcodes <- simulated_barcodes[match(unique(simulated_barcodes$Barcode), simulated_barcodes$Barcode),]

simulated_barcodes_ordered <- simulated_barcodes[match(individual_assignment$Barcode,simulated_barcodes$Barcode),]



# ### This table indicated that using a combination of softwares is better than any one
# ##### Test combinations of demultiplexing softwares and compare to the original ones to demonstrate advantage #####
# ### 1. Make calls based on at least half of them ###
# ## if >half say singlet => singlet, otherwise, majority
demultiplexing_twos <- apply(t(combn(demultiplexing_list, 2, simplify = TRUE)) , 1 , paste , collapse = "-" )
demultiplexing_threes <- apply(t(combn(demultiplexing_list, 3, simplify = TRUE)) , 1 , paste , collapse = "-" )
demultiplexing_fours <- apply(t(combn(demultiplexing_list, 4, simplify = TRUE)) , 1 , paste , collapse = "-" )
demultiplexing_fives <- apply(t(combn(demultiplexing_list, 5, simplify = TRUE)) , 1 , paste , collapse = "-" )
demultiplexing_sixes <- apply(t(combn(demultiplexing_list, 6, simplify = TRUE)) , 1 , paste , collapse = "-" )
demultiplexing_sevens <- apply(t(combn(demultiplexing_list, 7, simplify = TRUE)) , 1 , paste , collapse = "-" )
demultiplexing_eights <- apply(t(combn(demultiplexing_list, 8, simplify = TRUE)) , 1 , paste , collapse = "-" )

demultiplexing_combinations <- c(demultiplexing_twos, demultiplexing_threes, demultiplexing_fours, demultiplexing_fives, demultiplexing_sixes, demultiplexing_sevens, demultiplexing_eights)


## Initialize lists ##
demultiplexing_combos_results <- list()

common_assignment <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))
colnames(common_assignment) <- c("Barcode", demultiplexing_combinations)
common_assignment$Barcode <- singlets$Barcode


demultiplexing_combos_results[["majority_singlet"]][["droplet_type"]] <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))
demultiplexing_combos_results[["majority_singlet"]][["individual"]] <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))

colnames(demultiplexing_combos_results[["majority_singlet"]][["droplet_type"]]) <- c("Barcode", demultiplexing_combinations)
colnames(demultiplexing_combos_results[["majority_singlet"]][["individual"]]) <- c("Barcode", demultiplexing_combinations)

demultiplexing_combos_results[["majority_singlet"]][["droplet_type"]]$Barcode <- singlets$Barcode
demultiplexing_combos_results[["majority_singlet"]][["individual"]]$Barcode <- singlets$Barcode

demultiplexing_combos_results[["majority_doublet"]][["droplet_type"]] <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))
demultiplexing_combos_results[["majority_doublet"]][["individual"]] <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))

colnames(demultiplexing_combos_results[["majority_doublet"]][["droplet_type"]]) <- c("Barcode", demultiplexing_combinations)
colnames(demultiplexing_combos_results[["majority_doublet"]][["individual"]]) <- c("Barcode", demultiplexing_combinations)

demultiplexing_combos_results[["majority_doublet"]][["droplet_type"]]$Barcode <- singlets$Barcode
demultiplexing_combos_results[["majority_doublet"]][["individual"]]$Barcode <- singlets$Barcode


for (combo in demultiplexing_combinations){
	print(combo)
	softwares <- strsplit(combo,"-")[[1]]

	common_assignment[,combo] <- future_apply(individual_assignment[,..softwares], 1, function(y) names(which.max(table(y)))) ### This only works because all the individuals are numbers so will be first in table (before doublet)

	demultiplexing_combos_results[["majority_singlet"]][["droplet_type"]][,combo] <- ifelse(rowSums(singlets[,softwares]) >= length(softwares)/2 & rowSums(individual_assignment[,..softwares] == common_assignment[,combo]) >= length(softwares)/2, "singlet", "doublet")

	demultiplexing_combos_results[["majority_doublet"]][["droplet_type"]][,combo] <- ifelse(rowSums(singlets[,softwares]) > length(softwares)/2  & rowSums(individual_assignment[,..softwares] == common_assignment[,combo]) > length(softwares)/2, "singlet", "doublet")


	demultiplexing_combos_results[["majority_singlet"]][["individual"]][,combo] <- ifelse(demultiplexing_combos_results[["majority_singlet"]][["droplet_type"]][,combo] == "singlet", as.vector(t(common_assignment[,combo])), "doublet")

	demultiplexing_combos_results[["majority_doublet"]][["individual"]][,combo] <- ifelse(demultiplexing_combos_results[["majority_doublet"]][["droplet_type"]][,combo] == "singlet", as.vector(t(common_assignment[,combo])), "doublet")
}


which.max.simple=function(x,na.rm=TRUE,tie_value="NA"){
	if(na.rm){
		x=x[!is.na(x)]
	}
	if(length(x)==0){
		return(NA)
	}
	maxval=max(x)
	if(is.na(maxval)){
		return(NA)
	}
	if(sum(x %in% maxval) > 1){
		# Ties exist, figure out what to do with them.
		if(tie_value=="NA"){
			return(NA)
		}
		if(tie_value=="random"){
			tie_postions=which(x==maxval)
			return(sample(tie_postions,size=1))
		}
		if(tie_value=="first"){
			tie_postions=which(x==maxval)
			return(tie_postions[1])
		}
	} else{
		return(which.max(x))
	}
}

common_assignment2 <- list()
demultiplexing_combos_results2 <- list()


common_assignment2 <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))

colnames(common_assignment2) <- c("Barcode", demultiplexing_combinations)

common_assignment2$Barcode <- singlets$Barcode


demultiplexing_combos_results2[["majority_singlet"]][["individual"]] <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))

colnames(demultiplexing_combos_results2[["majority_singlet"]][["individual"]]) <- c("Barcode", demultiplexing_combinations)

demultiplexing_combos_results2[["majority_singlet"]][["individual"]]$Barcode <- singlets$Barcode

demultiplexing_combos_results2[["majority_doublet"]][["individual"]] <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))

colnames(demultiplexing_combos_results2[["majority_doublet"]][["individual"]]) <- c("Barcode", demultiplexing_combinations)

demultiplexing_combos_results2[["majority_doublet"]][["individual"]]$Barcode <- singlets$Barcode


for (combo in demultiplexing_combinations){
	print(combo)
	softwares <- strsplit(combo,"-")[[1]]

		common_assignment2[,combo] <- unlist(future_apply(individual_assignment[,..softwares], 1, function(y) ifelse(is.null(names(which.max.simple(table(y)[rownames(table(y)) != "doublet"][names(table(y)[rownames(table(y)) != "doublet"]) != "doublet"]))), NA, names(which.max.simple(table(y)[rownames(table(y)) != "doublet"][names(table(y)[rownames(table(y)) != "doublet"]) != "doublet"])))))

		demultiplexing_combos_results2[["majority_singlet"]][["individual"]][,combo] <- ifelse(demultiplexing_combos_results[["majority_singlet"]][["droplet_type"]][,combo] == "singlet", as.vector(t(common_assignment2[,combo])), "doublet")

		demultiplexing_combos_results2[["majority_doublet"]][["individual"]][,combo] <- ifelse(demultiplexing_combos_results[["majority_doublet"]][["droplet_type"]][,combo] == "singlet", as.vector(t(common_assignment2[,combo])), "doublet")
}




### 2. Binary classification for correct or not, then calcultae % correct, MCC and balanced accuracy
demultiplexing_combos_results[["majority_singlet"]][["correct"]] <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))
demultiplexing_combos_results[["majority_doublet"]][["correct"]] <- data.frame(matrix(ncol = length(demultiplexing_combinations) + 1, nrow = nrow(singlets)))

colnames(demultiplexing_combos_results[["majority_singlet"]][["correct"]]) <- c("Barcode", demultiplexing_combinations)
colnames(demultiplexing_combos_results[["majority_doublet"]][["correct"]]) <- c("Barcode", demultiplexing_combinations)

demultiplexing_combos_results[["majority_singlet"]][["correct"]]$Barcode <- singlets$Barcode
demultiplexing_combos_results[["majority_doublet"]][["correct"]]$Barcode <- singlets$Barcode

for (combo in demultiplexing_combinations){
	print(combo)
	demultiplexing_combos_results[["majority_singlet"]][["correct"]][,combo] <- ifelse(demultiplexing_combos_results[["majority_singlet"]][["droplet_type"]][,combo] == simulated_barcodes_ordered$DropletType & demultiplexing_combos_results[["majority_singlet"]][["individual"]][,combo] == simulated_barcodes_ordered$Individual, "Correct", "Incorrect")
	demultiplexing_combos_results[["majority_doublet"]][["correct"]][,combo] <- ifelse(demultiplexing_combos_results[["majority_doublet"]][["droplet_type"]][,combo] == simulated_barcodes_ordered$DropletType & demultiplexing_combos_results[["majority_doublet"]][["individual"]][,combo] == simulated_barcodes_ordered$Individual, "Correct", "Incorrect")
}



# pool <- "size8_SimulatedPool5"
# method <- "majority_singlet"
# soft <- "demuxlet-souporcell-vireo"
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


for (method in names(demultiplexing_combos_results)){
	Metrics[[method]] <- data.frame(matrix(nrow = length(demultiplexing_combinations), ncol = 13))
	colnames(Metrics[[method]]) <- c("software", "true_singlet", "false_singlet", "true_doublet", "false_doublet", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "npv", "ppv", "precision", "mcc")
	Metrics[[method]]$software <- demultiplexing_combinations
	rownames(Metrics[[method]]) <- demultiplexing_combinations
	Metrics[[method]]$pool <- args$pool
	Metrics[[method]]$method <- method

	for (soft in demultiplexing_combinations){
		## True singlet (TP)
		# Metrics[soft,]$true_singlet <- as.numeric(length(which(droplet_type[, soft] == "singlet" & droplet_type[, soft] == simulated_barcodes_ordered$DropletType & individual_assignment[, ..soft] == simulated_barcodes_ordered$Individual)))
		Metrics[[method]][soft,]$true_singlet <- as.numeric(length(which(demultiplexing_combos_results[[method]][["droplet_type"]][, soft] == "singlet" & demultiplexing_combos_results[[method]][["droplet_type"]][, soft] == simulated_barcodes_ordered$DropletType & demultiplexing_combos_results[[method]][["individual"]][, soft] == simulated_barcodes_ordered$Individual)))
		# Metrics[soft,]$true_singlet[Metrics[soft,]$true_singlet == 0] <- NA

		## False singlet (FP)
		Metrics[[method]][soft,]$false_singlet <- as.numeric(length(which(demultiplexing_combos_results[[method]][["droplet_type"]][, soft] == "singlet")) - Metrics[[method]][soft,]$true_singlet)
		# Metrics[soft,]$false_singlet[Metrics[soft,]$false_singlet == 0] <- NA

		### Update NA for pools where could not run (ie DoubletDecon)
		if (Metrics[[method]][soft,]$true_singlet == 0 & Metrics[[method]][soft,]$false_singlet == 0){
			Metrics[[method]][soft,]$true_doublet <- NA
			Metrics[[method]][soft,]$false_doublet <- NA
			Metrics[[method]][soft,]$accuracy
			Metrics[[method]][soft,]$sensitivity <- NA
			Metrics[[method]][soft,]$specificity <- NA
			Metrics[[method]][soft,]$npv <- NA
			Metrics[[method]][soft,]$ppv <- NA
			Metrics[[method]][soft,]$balanced_accuracy <- NA
			Metrics[[method]][soft,]$precision <- NA
			Metrics[[method]][soft,]$mcc <- NA
		} else{

			## True doublet (TN)
			Metrics[[method]][soft,]$true_doublet <- as.numeric(length(which(demultiplexing_combos_results[[method]][["droplet_type"]][, soft] == "doublet" & demultiplexing_combos_results[[method]][["droplet_type"]][, soft] == simulated_barcodes_ordered$DropletType)))


			## False doublet (FN)
			Metrics[[method]][soft,]$false_doublet <- as.numeric(length(which(demultiplexing_combos_results[[method]][["droplet_type"]][, soft] == "doublet")) - Metrics[[method]][soft,]$true_doublet)


			## Accuracy
			Metrics[[method]][soft,]$accuracy <- as.numeric((Metrics[[method]][soft,]$true_singlet + Metrics[[method]][soft,]$true_doublet)/nrow(demultiplexing_combos_results[[method]][["droplet_type"]]))


			## Sensitivity (TPR)
			Metrics[[method]][soft,]$sensitivity <- as.numeric(Metrics[[method]][soft,]$true_singlet/length(which(simulated_barcodes_ordered$DropletType == "singlet")))


			## Specificity (TNR)
			Metrics[[method]][soft,]$specificity <- as.numeric(Metrics[[method]][soft,]$true_doublet/length(which(simulated_barcodes_ordered$DropletType == "doublet")))


			## Negative Predictive Value
			Metrics[[method]][soft,]$npv <- as.numeric(Metrics[[method]][soft,]$true_doublet/length(which(demultiplexing_combos_results[[method]][["droplet_type"]][, soft] == "doublet")))


			## Positive Predictive Value
			Metrics[[method]][soft,]$ppv <- as.numeric(Metrics[[method]][soft,]$true_singlet/length(which(demultiplexing_combos_results[[method]][["droplet_type"]][, soft] == "singlet")))


			## Balanced Accuracy
			Metrics[[method]][soft,]$balanced_accuracy <- as.numeric((Metrics[[method]][soft,]$sensitivity + Metrics[[method]][soft,]$specificity)/2)


			## Precision
			Metrics[[method]][soft,]$precision <- as.numeric(Metrics[[method]][soft,]$true_singlet/(Metrics[[method]][soft,]$true_singlet + Metrics[[method]][soft,]$false_singlet))

			
			## MCC
			Metrics[[method]][soft,]$mcc <- as.numeric(((Metrics[[method]][soft,]$true_singlet * Metrics[[method]][soft,]$true_doublet) - (Metrics[[method]][soft,]$false_singlet * Metrics[[method]][soft,]$false_doublet))/sqrt((Metrics[[method]][soft,]$true_singlet + Metrics[[method]][soft,]$false_singlet) * (Metrics[[method]][soft,]$true_singlet + Metrics[[method]][soft,]$false_doublet) * (Metrics[[method]][soft,]$true_doublet + Metrics[[method]][soft,]$false_singlet) * (Metrics[[method]][soft,]$true_doublet + Metrics[[method]][soft,]$false_doublet)))
		}
	}
}



Metrics_df <- do.call(rbind, Metrics)
rownames(Metrics_df) <- NULL
Metrics_df <- as.data.table(Metrics_df)

Metrics_df$size <- factor(as.numeric(gsub("size","", Metrics_df$pool) %>% gsub("_.+", "", .)), levels = c(2,4,8,16,32,64,128))
Metrics_df$mt_percent <- ifelse(grepl("mt", Metrics_df$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_mt_","", Metrics_df$pool) %>% gsub("pctl", "", .)), 0)
Metrics_df$ambient_percent <- ifelse(grepl("ambient", Metrics_df$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_ambient_","", Metrics_df$pool) %>% gsub("pctl", "", .)), 0)
Metrics_df$subsampled <- ifelse(grepl("subsampled", Metrics_df$pool), "subsampled", "normal")
Metrics_df$uneven <- ifelse(grepl("unevenN", Metrics_df$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_unevenN_","", Metrics_df$pool) %>% gsub("pctl", "", .)), 0)

fwrite(Metrics_df, paste0(outdir,"/demultiplexing_metrics.tsv"), sep = "\t")
