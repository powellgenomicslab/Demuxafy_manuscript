library(tidyverse)
library(circlize)
library(plyr)
library(scales)
library(data.table)
library(future.apply)
library(ggpubr)
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




print("generating combinations")
doublet_detecting_twos <- apply(t(combn(doublet_detection_list, 2, simplify = TRUE)) , 1 , paste , collapse = "-" )
doublet_detecting_threes <- apply(t(combn(doublet_detection_list, 3, simplify = TRUE)) , 1 , paste , collapse = "-" )
doublet_detecting_fours <- apply(t(combn(doublet_detection_list, 4, simplify = TRUE)) , 1 , paste , collapse = "-" )
doublet_detecting_fives <- apply(t(combn(doublet_detection_list, 5, simplify = TRUE)) , 1 , paste , collapse = "-" )
doublet_detecting_sixes <- apply(t(combn(doublet_detection_list, 6, simplify = TRUE)) , 1 , paste , collapse = "-" )
doublet_detecting_sevens <- apply(t(combn(doublet_detection_list, 7, simplify = TRUE)) , 1 , paste , collapse = "-" )
doublet_detecting_eights <- apply(t(combn(doublet_detection_list, 8, simplify = TRUE)) , 1 , paste , collapse = "-" )

doublet_detecting_combinations <- c(doublet_detecting_twos, doublet_detecting_threes, doublet_detecting_fours, doublet_detecting_fives, doublet_detecting_sixes, doublet_detecting_sevens, doublet_detecting_eights)


## Initialize lists ##
doublet_detecting_combos_results <- list()

print("majority")
### Make a dataframe to record the calls with majority singlet
doublet_detecting_combos_results[["majority_singlet"]][["droplet_type"]] <- data.frame(matrix(ncol = length(doublet_detecting_combinations) + 1, nrow = nrow(singlets)))
colnames(doublet_detecting_combos_results[["majority_singlet"]][["droplet_type"]]) <- c("Barcode", doublet_detecting_combinations)
doublet_detecting_combos_results[["majority_singlet"]][["droplet_type"]]$Barcode <- singlets$Barcode


### Make a dataframe to record the calls with majority doublet
doublet_detecting_combos_results[["majority_doublet"]][["droplet_type"]] <- data.frame(matrix(ncol = length(doublet_detecting_combinations) + 1, nrow = nrow(singlets)))
colnames(doublet_detecting_combos_results[["majority_doublet"]][["droplet_type"]]) <- c("Barcode", doublet_detecting_combinations)
doublet_detecting_combos_results[["majority_doublet"]][["droplet_type"]]$Barcode <- singlets$Barcode


for (combo in doublet_detecting_combinations){
	print(combo)
	softwares <- strsplit(combo,"-")[[1]]

		doublet_detecting_combos_results[["majority_singlet"]][["droplet_type"]][,combo] <- ifelse(rowSums(singlets[,softwares]) >= length(softwares)/2, "singlet", "doublet")
		doublet_detecting_combos_results[["majority_doublet"]][["droplet_type"]][,combo] <- ifelse(rowSums(singlets[,softwares]) > length(softwares)/2, "singlet", "doublet")
}


print("binary")
### 2. Binary classification for correct or not, then calcultae % correct, MCC and balanced accuracy
doublet_detecting_combos_results[["majority_singlet"]][["correct"]] <- data.frame(matrix(ncol = length(doublet_detecting_combinations) + 1, nrow = nrow(singlets)))
doublet_detecting_combos_results[["majority_doublet"]][["correct"]] <- data.frame(matrix(ncol = length(doublet_detecting_combinations) + 1, nrow = nrow(singlets)))

colnames(doublet_detecting_combos_results[["majority_singlet"]][["correct"]]) <- c("Barcode", doublet_detecting_combinations)
colnames(doublet_detecting_combos_results[["majority_doublet"]][["correct"]]) <- c("Barcode", doublet_detecting_combinations)

doublet_detecting_combos_results[["majority_singlet"]][["correct"]]$Barcode <- singlets$Barcode
doublet_detecting_combos_results[["majority_doublet"]][["correct"]]$Barcode <- singlets$Barcode

for (combo in doublet_detecting_combinations){
	print(combo)
	doublet_detecting_combos_results[["majority_singlet"]][["correct"]][,combo] <- ifelse(doublet_detecting_combos_results[["majority_singlet"]][["droplet_type"]][,combo] == simulated_barcodes_ordered$droplettype, "Correct", "Incorrect")
	doublet_detecting_combos_results[["majority_doublet"]][["correct"]][,combo] <- ifelse(doublet_detecting_combos_results[["majority_doublet"]][["droplet_type"]][,combo] == simulated_barcodes_ordered$droplettype, "Correct", "Incorrect")
}



saveRDS(doublet_detecting_combos_results, paste0(outdir, "/doublet_detecting_combos_classifications.rds"))

print("metrics")

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

for (method in names(doublet_detecting_combos_results)){
		Metrics[[method]] <- data.frame(matrix(nrow = length(doublet_detecting_combinations), ncol = 13))
		colnames(Metrics[[method]]) <- c("software", "true_singlet", "false_singlet", "true_doublet", "false_doublet", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "npv", "ppv", "precision", "mcc")
		Metrics[[method]]$software <- doublet_detecting_combinations
		rownames(Metrics[[method]]) <- doublet_detecting_combinations
		Metrics[[method]]$pool <- args$pool
		Metrics[[method]]$method <- method

		for (soft in doublet_detecting_combinations){
			## True singlet (TP)
			Metrics[[method]][soft,]$true_singlet <- as.numeric(length(which(doublet_detecting_combos_results[[method]][["droplet_type"]][, soft] == "singlet" & doublet_detecting_combos_results[[method]][["droplet_type"]][, soft] == simulated_barcodes_ordered$DropletType)))

			## False singlet (FP)
			Metrics[[method]][soft,]$false_singlet <- as.numeric(length(which(doublet_detecting_combos_results[[method]][["droplet_type"]][, soft] == "singlet")) - Metrics[[method]][soft,]$true_singlet)

			### Update NA for pools where could not run (ie DoubletDecon)
			if (Metrics[[method]][soft,]$true_singlet == 0 & Metrics[[method]][soft,]$false_singlet == 0){
				Metrics[[method]][soft,]$true_doublet <- NA
				Metrics[[method]][soft,]$false_doublet <- NA
				Metrics[[method]][soft,]$accuracy <- NA
				Metrics[[method]][soft,]$sensitivity <- NA
				Metrics[[method]][soft,]$specificity <- NA
				Metrics[[method]][soft,]$npv <- NA
				Metrics[[method]][soft,]$ppv <- NA
				Metrics[[method]][soft,]$balanced_accuracy <- NA
				Metrics[[method]][soft,]$precision <- NA
				Metrics[[method]][soft,]$mcc <- NA
			} else{

				## True doublet (TN)
				Metrics[[method]][soft,]$true_doublet <- as.numeric(length(which(doublet_detecting_combos_results[[method]][["droplet_type"]][, soft] == "doublet" & doublet_detecting_combos_results[[method]][["droplet_type"]][, soft] == simulated_barcodes_ordered$DropletType)))


				## False doublet (FN)
				Metrics[[method]][soft,]$false_doublet <- as.numeric(length(which(doublet_detecting_combos_results[[method]][["droplet_type"]][, soft] == "doublet")) - Metrics[[method]][soft,]$true_doublet)


				## Accuracy
				Metrics[[method]][soft,]$accuracy <- as.numeric((Metrics[[method]][soft,]$true_singlet + Metrics[[method]][soft,]$true_doublet)/nrow(doublet_detecting_combos_results[[method]][["droplet_type"]]))


				## Sensitivity (TPR)
				Metrics[[method]][soft,]$sensitivity <- as.numeric(Metrics[[method]][soft,]$true_singlet/length(which(simulated_barcodes_ordered$DropletType == "singlet")))


				## Specificity (TNR)
				Metrics[[method]][soft,]$specificity <- as.numeric(Metrics[[method]][soft,]$true_doublet/length(which(simulated_barcodes_ordered$DropletType == "doublet")))


				## Negative Predictive Value
				Metrics[[method]][soft,]$npv <- as.numeric(Metrics[[method]][soft,]$true_doublet/length(which(doublet_detecting_combos_results[[method]][["droplet_type"]][, soft] == "doublet")))


				## Positive Predictive Value
				Metrics[[method]][soft,]$ppv <- as.numeric(Metrics[[method]][soft,]$true_singlet/length(which(doublet_detecting_combos_results[[method]][["droplet_type"]][, soft] == "singlet")))


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

fwrite(Metrics_df, paste0(outdir, "/doublet_detecting_metrics.tsv"), sep = "\t")
