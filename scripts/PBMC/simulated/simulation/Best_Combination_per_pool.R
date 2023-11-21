library(dsLib)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(plyr)
library(scales)
library(viridis)
library(RColorBrewer)
library(ggforce)
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


singlets_updated <- singlets %>% mutate_all(as.character)
singlets_updated[singlets_updated == "0"] <- "doublet"
singlets_updated[singlets_updated == "1"] <- "singlet"


##### Read in the cell_info barcode files #####
simulated_barcodes <- fread(paste0(dir, "/",args$pool,"/matrix_out/barcodes.tsv.gz"), sep = "\t", header = FALSE, col.names = c("Barcode"))

simulated_barcodes$DropletType <- ifelse(grepl(":", simulated_barcodes$Barcode), "doublet","singlet")
simulated_barcodes$DoubletType_DoubletType <- ifelse(gsub("[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", gsub(":.+", "", simulated_barcodes$Barcode)) == gsub(".+:[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", simulated_barcodes$Barcode), "homogenic_doublet", ifelse(grepl(":", simulated_barcodes$Barcode), "heterogenic_doublet","singlet"))



##### Join barcodes with individual identifiers #####
simulated_barcodes$Individual <- gsub("[A,C,G,T,L,M,N,O,P,Q,R]+-", "", gsub(":.+", "", simulated_barcodes$Barcode))
simulated_barcodes$Individual2 <- ifelse(!grepl(":", simulated_barcodes$Barcode), NA, gsub(".+:[A,C,G,T,L,M,N,O,P,Q,R]+-", "", simulated_barcodes$Barcode))

simulated_barcodes$Individual <- ifelse(simulated_barcodes$DropletType == "doublet", "doublet", simulated_barcodes$Individual)
simulated_barcodes <- simulated_barcodes[match(unique(simulated_barcodes$Barcode), simulated_barcodes$Barcode),]

simulated_barcodes_ordered <- simulated_barcodes[match(individual_assignment$Barcode,simulated_barcodes$Barcode),]



# ### This table indicated that using a combination of softwares is better than any one
##### Test combinations of demultiplexing softwares and compare to the original ones to demonstrate advantage #####
### 1. Make calls based on at least half of them ###
## if >half say singlet => singlet, otherwise, majority
twos <- apply(t(combn(softwares, 2, simplify = TRUE)) , 1 , paste , collapse = "-" )
threes <- apply(t(combn(softwares, 3, simplify = TRUE)) , 1 , paste , collapse = "-" )
fours <- apply(t(combn(softwares, 4, simplify = TRUE)) , 1 , paste , collapse = "-" )
fives <- apply(t(combn(softwares, 5, simplify = TRUE)) , 1 , paste , collapse = "-" )
sixes <- apply(t(combn(softwares, 6, simplify = TRUE)) , 1 , paste , collapse = "-" )
sevens <- apply(t(combn(softwares, 7, simplify = TRUE)) , 1 , paste , collapse = "-" )
eights <- apply(t(combn(softwares, 8, simplify = TRUE)) , 1 , paste , collapse = "-" )
nines <- apply(t(combn(softwares, 9, simplify = TRUE)) , 1 , paste , collapse = "-" )
tens <- apply(t(combn(softwares, 10, simplify = TRUE)) , 1 , paste , collapse = "-" )
elevens <- apply(t(combn(softwares, 11, simplify = TRUE)) , 1 , paste , collapse = "-" )
twelves <- apply(t(combn(softwares, 12, simplify = TRUE)) , 1 , paste , collapse = "-" )
thirteens <- apply(t(combn(softwares, 13, simplify = TRUE)) , 1 , paste , collapse = "-" )
fourteens <- apply(t(combn(softwares, 14, simplify = TRUE)) , 1 , paste , collapse = "-" )
fifteens <- apply(t(combn(softwares, 15, simplify = TRUE)) , 1 , paste , collapse = "-" )
sixteens <- apply(t(combn(softwares, 16, simplify = TRUE)) , 1 , paste , collapse = "-" )

combinations_all <- c(twos, 
					threes, 
					fours,
					fives, 
					sixes, 
					sevens,
					eights,
					nines,
					tens,
					elevens,
					twelves,
					thirteens,
					fourteens,
					fifteens,
					sixteens)


combinations <- c()

for (combo in combinations_all){
	softs <- unlist(as.vector(strsplit(combo, "-")), use.names=FALSE)
	if (!all(softs %in% demultiplexing_list) & !all(softs %in% doublet_detection_list)){
		combinations <- c(combinations,combo)
	}
}


## Initialize lists ##
common_assignment <- list()
combos_results <- list()


common_assignment <- data.frame(matrix(ncol = length(combinations) + 1, nrow = nrow(singlets)))

colnames(common_assignment) <- c("Barcode", combinations)

common_assignment$Barcode <- singlets$Barcode


combos_results[["majority_singlet"]][["droplet_type"]] <- data.frame(matrix(ncol = length(combinations) + 1, nrow = nrow(singlets)))
combos_results[["majority_singlet"]][["individual"]] <- data.frame(matrix(ncol = length(combinations) + 1, nrow = nrow(singlets)))

colnames(combos_results[["majority_singlet"]][["droplet_type"]]) <- c("Barcode", combinations)
colnames(combos_results[["majority_singlet"]][["individual"]]) <- c("Barcode", combinations)

combos_results[["majority_singlet"]][["droplet_type"]]$Barcode <- singlets$Barcode
combos_results[["majority_singlet"]][["individual"]]$Barcode <- singlets$Barcode

combos_results[["majority_doublet"]][["droplet_type"]] <- data.frame(matrix(ncol = length(combinations) + 1, nrow = nrow(singlets)))
combos_results[["majority_doublet"]][["individual"]] <- data.frame(matrix(ncol = length(combinations) + 1, nrow = nrow(singlets)))

colnames(combos_results[["majority_doublet"]][["droplet_type"]]) <- c("Barcode", combinations)
colnames(combos_results[["majority_doublet"]][["individual"]]) <- c("Barcode", combinations)

combos_results[["majority_doublet"]][["droplet_type"]]$Barcode <- singlets$Barcode
combos_results[["majority_doublet"]][["individual"]]$Barcode <- singlets$Barcode


for (combo in combinations){
	print(combo)
	software_list <- strsplit(combo,"-")[[1]]
	software_demultiplexing_list <- software_list[software_list %in% demultiplexing_list]
	software_doublet_detection_list <- software_list[software_list %in% doublet_detection_list]


		common_assignment[,combo] <- future_apply(individual_assignment[,..software_demultiplexing_list], 1, function(y) names(which.max(table(y)))) ### This only works because all the individuals are numbers so will be first in table (before doublet)

		combos_results[["majority_singlet"]][["droplet_type"]][,combo] <- ifelse(rowSums(singlets[,software_list]) >= length(software_list)/2 & rowSums(individual_assignment[,..software_demultiplexing_list] == common_assignment[,combo]) >= length(software_demultiplexing_list)/2, "singlet", "doublet")

		combos_results[["majority_doublet"]][["droplet_type"]][,combo] <- ifelse(rowSums(singlets[,software_list]) > length(software_list)/2  & rowSums(individual_assignment[,..software_demultiplexing_list] == common_assignment[,combo]) > length(software_demultiplexing_list)/2, "singlet", "doublet")


		combos_results[["majority_singlet"]][["individual"]][,combo] <- ifelse(combos_results[["majority_singlet"]][["droplet_type"]][,combo] == "singlet", as.vector(t(common_assignment[,combo])), "doublet")

		combos_results[["majority_doublet"]][["individual"]][,combo] <- ifelse(combos_results[["majority_doublet"]][["droplet_type"]][,combo] == "singlet", as.vector(t(common_assignment[,combo])), "doublet")
}



### 2. Binary classification for correct or not, then calculate % correct, MCC and balanced accuracy
combos_results[["majority_singlet"]][["correct"]] <- data.frame(matrix(ncol = length(combinations) + 1, nrow = nrow(singlets)))
combos_results[["majority_doublet"]][["correct"]] <- data.frame(matrix(ncol = length(combinations) + 1, nrow = nrow(singlets)))

colnames(combos_results[["majority_singlet"]][["correct"]]) <- c("Barcode", combinations)
colnames(combos_results[["majority_doublet"]][["correct"]]) <- c("Barcode", combinations)

combos_results[["majority_singlet"]][["correct"]]$Barcode <- singlets$Barcode
combos_results[["majority_doublet"]][["correct"]]$Barcode <- singlets$Barcode

for (combo in combinations){
	print(combo)
	combos_results[["majority_singlet"]][["correct"]][,combo] <- ifelse(combos_results[["majority_singlet"]][["droplet_type"]][,combo] == simulated_barcodes_ordered$DropletType & combos_results[["majority_singlet"]][["individual"]][,combo] == simulated_barcodes_ordered$Individual, "Correct", "Incorrect")
	combos_results[["majority_doublet"]][["correct"]][,combo] <- ifelse(combos_results[["majority_doublet"]][["droplet_type"]][,combo] == simulated_barcodes_ordered$DropletType & combos_results[["majority_doublet"]][["individual"]][,combo] == simulated_barcodes_ordered$Individual, "Correct", "Incorrect")
}



saveRDS(combos_results, paste0(outdir, "/Combination_Results.rds"))
combos_results <- readRDS(paste0(outdir, "/Combination_Results.rds"))




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

for (method in names(combos_results)){
	Metrics[[method]] <- data.frame(matrix(nrow = length(combinations), ncol = 13))
	colnames(Metrics[[method]]) <- c("software", "true_singlet", "false_singlet", "true_doublet", "false_doublet", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "npv", "ppv", "precision", "mcc")
	Metrics[[method]]$software <- combinations
	rownames(Metrics[[method]]) <- combinations
	Metrics[[method]]$pool <- args$pool
	Metrics[[method]]$method <- method

	for (soft in combinations){
		## True singlet (TP)
		# Metrics[soft,]$true_singlet <- as.numeric(length(which(droplet_type[, soft] == "singlet" & droplet_type[, soft] == simulated_barcodes_ordered$DropletType & individual_assignment[, ..soft] == simulated_barcodes_ordered$Individual)))
		Metrics[[method]][soft,]$true_singlet <- as.numeric(length(which(combos_results[[method]][["droplet_type"]][, soft] == "singlet" & combos_results[[method]][["droplet_type"]][, soft] == simulated_barcodes_ordered$DropletType & combos_results[[method]][["individual"]][, soft] == simulated_barcodes_ordered$Individual)))
		# Metrics[soft,]$true_singlet[Metrics[soft,]$true_singlet == 0] <- NA

		## False singlet (FP)
		Metrics[[method]][soft,]$false_singlet <- as.numeric(length(which(combos_results[[method]][["droplet_type"]][, soft] == "singlet")) - Metrics[[method]][soft,]$true_singlet)
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
			Metrics[[method]][soft,]$true_doublet <- as.numeric(length(which(combos_results[[method]][["droplet_type"]][, soft] == "doublet" & combos_results[[method]][["droplet_type"]][, soft] == simulated_barcodes_ordered$DropletType)))


			## False doublet (FN)
			Metrics[[method]][soft,]$false_doublet <- as.numeric(length(which(combos_results[[method]][["droplet_type"]][, soft] == "doublet")) - Metrics[[method]][soft,]$true_doublet)


			## Accuracy
			Metrics[[method]][soft,]$accuracy <- as.numeric((Metrics[[method]][soft,]$true_singlet + Metrics[[method]][soft,]$true_doublet)/nrow(combos_results[[method]][["droplet_type"]]))


			## Sensitivity (TPR)
			Metrics[[method]][soft,]$sensitivity <- as.numeric(Metrics[[method]][soft,]$true_singlet/length(which(simulated_barcodes_ordered$DropletType == "singlet")))


			## Specificity (TNR)
			Metrics[[method]][soft,]$specificity <- as.numeric(Metrics[[method]][soft,]$true_doublet/length(which(simulated_barcodes_ordered$DropletType == "doublet")))


			## Negative Predictive Value
			Metrics[[method]][soft,]$npv <- as.numeric(Metrics[[method]][soft,]$true_doublet/length(which(combos_results[[method]][["droplet_type"]][, soft] == "doublet")))


			## Positive Predictive Value
			Metrics[[method]][soft,]$ppv <- as.numeric(Metrics[[method]][soft,]$true_singlet/length(which(combos_results[[method]][["droplet_type"]][, soft] == "singlet")))


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

fwrite(Metrics_df, paste0(outdir, "/demultiplexing_doublet_detecting_metrics.tsv"), sep = "\t")

setkey(Metrics_df, mcc)


##### Read in the single software results to  plot on the same figure
Metrics_df_single_soft <- fread(paste0(args$out,"/SimulatedOverlap/", args$pool, "/single_software_metrics.tsv"))

Metrics_df_single_soft$method <- "single_software"
Metrics_df_single_soft <- Metrics_df_single_soft[,colnames(Metrics_df), with = FALSE]
Metrics_df_combined <- rbind(Metrics_df, Metrics_df_single_soft)

fwrite(Metrics_df_combined, paste0(outdir, "/demultiplexing_doublet_detecting_metrics_w_single_softs.tsv"), sep = "\t")

