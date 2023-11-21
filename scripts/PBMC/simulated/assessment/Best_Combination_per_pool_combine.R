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


##### Set up Functions #####
save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}


###### Set up Directories #####
dir <- "/path/to/output/PBMC/simulation/SimulatedPools/SimulatedOverlap/"
droplet_type_dir <- paste0(dir,"/DropletAnnotation/")
outdir <- paste0(dir, "Best_Combination_Demultiplexing_DoubletDetecting/")
dir.create(outdir)


##### Set up colors for plotting #####
downsample_colors <- c("#C7EAE5", "#35978F")
mt_percent_colors <- c("0" = "white", "5" = "#D0E0D3", "10" = "#A3C4A6", "25" = "#1E6D2D")
ambient_percent_colors <- c("0" = "white", "10" = "#FDEDD2", "20" = "#FBD491", "50" = "#E6AB00")
uneven_n_colors <- c("0" = "", "50" = "", "75" = "", "95" = "")


##### Set up universal 
demultiplexing_list <- c("demuxalot", "demuxalot_refined","demuxlet", "dropulation", "freemuxlet","scSplit","souporcell","vireo")
names(demultiplexing_list) <- c("Demuxalot", "Demuxalot (refined)","Demuxlet", "Dropulation", "Freemuxlet","ScSplit","Souporcell","Vireo")
doublet_detection_list <- c("DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder","scds",  "scrublet", "solo")
names(doublet_detection_list) <- c("DoubletDecon", "DoubletDetection","DoubletFinder",  "ScDblFinder", "Scds", "Scrublet", "Solo")
softwares <- c(demultiplexing_list, doublet_detection_list)

softwares_dt <- data.table(software = softwares, Software = names(softwares))


software_colors <- c("#575E57", "#008751", "#5DBB4D", "#EFE305", "#F8BE33", "#F97919", "#D21F09", "#8D132B", "#F980BF", "#EAA1D1", "#BAAAFF", "#9B64E0", "#6B82E3", "#63D4FE", "#39BCBC", "#C4E4DF", "#E2E2E2")
names(software_colors) <- names(softwares)


##### Get list of files #####
### Make a list of the pools
pools <- dir(dir, pattern = "size\\d+_SimulatedPool\\d+")
pools <- pools[!grepl("size128_", pools)]
pools <- pools[!grepl("size64_", pools)]
pools <- pools[!grepl("size32_", pools)]
pools <- pools[!grepl("mt", pools)]
pools <- pools[!grepl("ambient_", pools)]
pools <- pools[!grepl("uneven", pools)]
pools <- pools[!grepl("subsampled", pools)]


### paste together pool list with batch numbers
files <- c()

for (pool in pools){
	for (batch in 1:651){
		files <- c(files, paste0(dir,pool, "/", batch, "/demultiplexing_doublet_detecting_metrics.tsv"))
	}
}

results_list <- lapply(files, function(file){
	fread(file)
})

Metrics_df <- do.call(rbind, results_list)





setkey(Metrics_df, mcc)

summary <- Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0,.(MCC_Mean=mean(mcc, na.rm = TRUE), Balanced_Accuracy_Mean = mean(balanced_accuracy, na.rm = TRUE)),.(size,software,method)];
head(summary[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 2], n = 5)
head(summary[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 4], n = 5)
head(summary[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 8], n = 5)
head(summary[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 16], n = 5)


Metrics_df[, max(mcc, na.rm = T), by = size]
Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal", .SD[mcc %in% max(mcc, na.rm = T)], by=size]

fwrite(Metrics_df_combined, paste0(outdir,"demultiplexing_doublet_detecting_metrics.tsv"), sep = "\t")
Metrics_df_combined <-fread(paste0(outdir,"demultiplexing_doublet_detecting_metrics.tsv"), sep = "\t")


##### Read in the single software results to  plot on the same figure
Metrics_df_single_soft <- fread(dir, "/single_software_metrics.tsv")

Metrics_df_single_soft$method <- "single_software"
Metrics_df_single_soft <- Metrics_df_single_soft[,colnames(Metrics_df), with = FALSE]
Metrics_df_combined <- rbind(Metrics_df, Metrics_df_single_soft)

fwrite(Metrics_df_combined, paste0(outdir,"demultiplexing_doublet_detecting_metrics_w_single_softs.tsv"), sep = "\t")


Metrics_df_combined_noknowndoublets <- Metrics_df_combined[!grepl("scDblFinder_known_doublets", software)]

fwrite(Metrics_df_combined[!grepl("scDblFinder_known_doublets", software) & method != "single_software"], paste0(outdir,"demultiplexing_doublet_detecting_metrics_no_knowndoublets.tsv"), sep = "\t")


summary_combined <- Metrics_df_combined_noknowndoublets[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0,.(MCC_Mean=mean(mcc, na.rm = TRUE), Balanced_Accuracy_Mean = mean(balanced_accuracy, na.rm = TRUE)),.(size,software,method)];
head(summary_combined[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 2], n = 5)
head(summary_combined[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 4], n = 5)
head(summary_combined[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 8], n = 5)
head(summary_combined[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 16], n = 5)


### without ref genotypes 
summary_combined_no_ref <- Metrics_df_combined_noknowndoublets[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0 & !grepl("demuxalot", software) & !grepl("demuxalot_refined", software) & !grepl("demuxlet", software) & !grepl("dropulation", software) & !grepl("dropulation", software),.(MCC_Mean=mean(mcc, na.rm = TRUE), Balanced_Accuracy_Mean = mean(balanced_accuracy, na.rm = TRUE)),.(size,software,method)];
head(summary_combined_no_ref[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 2], n = 5)
head(summary_combined_no_ref[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 4], n = 5)
head(summary_combined_no_ref[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 8], n = 5)
head(summary_combined_no_ref[!is.na(MCC_Mean)][rev(order(MCC_Mean))][size == 16], n = 5)




##### Make of MCC for all softwares across individuals #####
pNPV_single <- ggplot(Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal"], aes(size, npv, color = software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_wrap(~ method) +
	theme_classic() +
	ylab("NPV") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6)

save_figs(pNPV_single,  paste0(outdir, "npv_multi_software_w_single"), width =30, height = 8)


##### Make of MCC for all softwares across individuals #####
pPPV_single <- ggplot(Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal"], aes(size, ppv, color = software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_wrap(~ method, scales = ) +
	theme_classic() +
	ylab("PPV") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6)

save_figs(pPPV_single,  paste0(outdir, "ppv_multi_software_w_single"), width =30, height = 20)

##### Make of MCC for all softwares across individuals #####
pMCC_single <- ggplot(Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & software %in% c(demultiplexing_list, demultiplexing_combinations)], aes(size, mcc, color = software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_wrap(~ method, scales = ) +
	theme_classic() +
	ylab("MCC") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6)

save_figs(pMCC_single,  paste0(outdir, "mcc_multi_software_w_single"), width =30, height = 20)



Metrics_df_combined_noknowndoublets[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0 & software %in% c(demultiplexing_list, demultiplexing_combinations), .SD[ppv %in% max(ppv, na.rm = T)], by=size]
summary_w_single <- Metrics_df_combined_noknowndoublets[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0  ,.(ppv_mean = mean(ppv, na.rm = TRUE), mean_N = mean(true_singlet + false_singlet + true_doublet + false_doublet), min_N = min(true_singlet + false_singlet + true_doublet + false_doublet), max_N = max(true_singlet + false_singlet + true_doublet + false_doublet)),.(size,software,method)]
summary_w_single[, .SD[ppv_mean %in% max(ppv_mean, na.rm = T)], by=size]




Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" ,][size == 128 & software == "demuxlet-freemuxlet-vireo"]
summary_w_single_balance <- Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" ,.(balanced_accuracy_mean = mean(balanced_accuracy, na.rm = TRUE), mcc_mean = mean(mcc, na.rm = TRUE), mean_N = mean(true_singlet + false_singlet + true_doublet + false_doublet), min_N = min(true_singlet + false_singlet + true_doublet + false_doublet), max_N = max(true_singlet + false_singlet + true_doublet + false_doublet)),.(size,software,method)]
summary_w_single_balance[, .SD[balanced_accuracy_mean %in% max(balanced_accuracy_mean, na.rm = T)], by=size]
summary_w_single_balance[, .SD[mcc_mean %in% max(mcc_mean, na.rm = T)], by=size]
summary_w_single_balance[software == "freemuxlet-souporcell-vireo-scds" & method == "majority_doublet"]
summary_w_single_balance[software == "demuxlet-souporcell-vireo-scds" & method == "majority_doublet"]




##### Make images for final figures #####
## Will use the top software for each time and plot across a three sizes since both are better than any alone anyhow
summary_recommended <- list()
summary_recommended[["remove_doublets"]] <- unique(summary_w_single[, .SD[ppv_mean %in% max(ppv_mean, na.rm = T)], by=size][size %in% c(2,4,8,16),c("software", "method")])
summary_recommended[["balance"]] <- unique(summary_w_single_balance[, .SD[mcc_mean %in% max(mcc_mean, na.rm = T)], by=size][size %in% c(2,4,8,16),c("software", "method")])
## Note, odd number classifiers are the same either way (majority singlet and majority doublet), so both are in the dataframe so will remove one for those that are duplicated



### Get 
Metrics_df_recommended <- list()
Metrics_df_recommended[["remove_doublets"]] <- Metrics_df_combined[summary_recommended[["remove_doublets"]], on = c("software", "method")][size %in% c(2,4,8,16)]
Metrics_df_recommended[["balance"]] <- Metrics_df_combined[summary_recommended[["balance"]], on = c("software", "method")][size %in% c(2,4,8,16)]


for (software in unique(summary_recommended[["remove_doublets"]]$software)){
	softwares <- strsplit(software,"-")[[1]]
	Metrics_df_recommended[["remove_doublets"]] <- rbind(Metrics_df_recommended[["remove_doublets"]], Metrics_df_combined[software %in% softwares & size %in% c(2,4,8,16)])
}
Metrics_df_recommended[["remove_doublets"]]  <- unique(Metrics_df_recommended[["remove_doublets"]] )

for (software in unique(summary_recommended[["balance"]]$software)){
	softwares <- strsplit(software,"-")[[1]]
	Metrics_df_recommended[["balance"]] <- rbind(Metrics_df_recommended[["balance"]], Metrics_df_combined[software %in% softwares& size %in% c(2,4,8,16)])
}
Metrics_df_recommended[["balance"]] <- unique(Metrics_df_recommended[["balance"]])


summary_w_single_balance <- Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal",.(balanced_accuracy_mean = mean(balanced_accuracy, na.rm = TRUE), mcc_mean = mean(mcc, na.rm = TRUE), mean_N = mean(true_singlet + false_singlet + true_doublet + false_doublet), min_N = min(true_singlet + false_singlet + true_doublet + false_doublet), max_N = max(true_singlet + false_singlet + true_doublet + false_doublet)),.(size,software,method)][size %in% c(2,4,8,16)]
summary_w_single <- Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal",.(ppv_mean = mean(ppv, na.rm = TRUE), mean_N = mean(true_singlet + false_singlet + true_doublet + false_doublet), min_N = min(true_singlet + false_singlet + true_doublet + false_doublet), max_N = max(true_singlet + false_singlet + true_doublet + false_doublet)),.(size,software,method)][size %in% c(2,4,8,16)]

Metrics_df_recommended[["remove_doublets"]] <- summary_w_single[Metrics_df_recommended[["remove_doublets"]], on = c("size", "software", "method")]
Metrics_df_recommended[["balance"]] <- summary_w_single_balance[Metrics_df_recommended[["balance"]], on = c("size", "software", "method")]


Metrics_df_recommended[["balance"]]$software <- factor(Metrics_df_recommended[["balance"]]$software, levels = unique(Metrics_df_recommended[["balance"]][order(Metrics_df_recommended[["balance"]]$mcc_mean)]$software))
Metrics_df_recommended[["remove_doublets"]]$software <- factor(Metrics_df_recommended[["remove_doublets"]]$software, unique(Metrics_df_recommended[["remove_doublets"]][order(Metrics_df_recommended[["remove_doublets"]]$ppv_mean)]$software))


#### Statistically test between groups ####
ttest_dt_list <- list()

for (method in names(Metrics_df_recommended)){
	dt2 <- list()
	for (combo in unique(grep("-", Metrics_df_recommended[[method]]$software, value = TRUE))){
		print(combo)
		softwares <- strsplit(combo,"-")[[1]]
		dt <- list()
		for (n in unique(Metrics_df_recommended[[method]][software == combo]$size)){
			dt[[method]][[as.character(n)]] <- data.table(software1 = rep(combo, length(softwares)), software2 = softwares, size = n)
			dt[[method]][[as.character(n)]]$statistic <- 0
			dt[[method]][[as.character(n)]]$p <- 0

			for (row in 1:nrow(dt[[method]][[as.character(n)]])){
				soft <- dt[[method]][[as.character(n)]]$software2[row]
				if (method == "balance"){
					ttest <- t.test(Metrics_df_recommended[[method]][software == combo & size == n]$mcc, Metrics_df_recommended[[method]][software == soft & size == n]$mcc, alternative = "greater")
					dt[[method]][[as.character(n)]][row,"statistic"] <- ttest$statistic
					dt[[method]][[as.character(n)]][row,"p"] <- ttest$p.value

				} else if (method == "remove_doublets"){
					ttest <- t.test(Metrics_df_recommended[[method]][software == combo & size == n]$ppv, Metrics_df_recommended[[method]][software == soft & size == n]$ppv, alternative = "greater")
					dt[[method]][[as.character(n)]][row,"statistic"] <- ttest$statistic
					dt[[method]][[as.character(n)]][row,"p"] <- ttest$p.value
				}
			}
		}
	dt2[[method]][[combo]] <- rbindlist(dt[[method]])
	}
	ttest_dt_list[[method]] <- rbindlist(dt2[[method]])
}



#### Statistically test between groups ####
##### Ftest for all softwares vs combos #####
summary_recommended_n <- list()
summary_recommended_n[["remove_doublets"]] <- unique(summary_w_single[, .SD[ppv_mean %in% max(ppv_mean, na.rm = T)], by=size][,c("software", "method", "size")])
summary_recommended_n[["balance"]] <- unique(summary_w_single_balance[, .SD[mcc_mean %in% max(mcc_mean, na.rm = T)], by=size][,c("software", "method", "size")])

Metrics_df_ftest <- list()
Metrics_df_ftest[["remove_doublets"]] <- rbind(Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & method == "single_software" & software %in% c(demultiplexing_list, doublet_detection_list)], Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" ][summary_recommended_n[["remove_doublets"]], on = c("size", "software", "method")])
Metrics_df_ftest[["balance"]] <- rbind(Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & method == "single_software" & software %in% c(demultiplexing_list, doublet_detection_list)], Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" ][summary_recommended_n[["balance"]], on = c("size", "software", "method")])



ttest_dt_list_all <- list()
ftest_dt_list_all <- list()

dtf <- list()
dt <- list()

dt2 <- list()
dt2f <- list()



for (method in names(Metrics_df_ftest)){
	for (combo in unique(grep("-", Metrics_df_ftest[[method]]$software, value = TRUE))){
		print(combo)
		for (n in unique(Metrics_df_ftest[[method]][software == combo]$size)){
			softwares <- unique(Metrics_df_ftest[[method]][method == "single_software"]$software)
			dtf[[method]][[as.character(n)]] <- data.table(software1 = rep(combo, length(softwares)), software2 = softwares, size = n)
			dtf[[method]][[as.character(n)]]$statistic <- as.numeric(NA)
			dtf[[method]][[as.character(n)]]$p <- as.numeric(NA)

			dt[[method]][[as.character(n)]] <- dtf[[method]][[as.character(n)]]

			for (row in 1:nrow(dtf[[method]][[as.character(n)]])){
				soft <- dtf[[method]][[as.character(n)]]$software2[row]
				print(soft)
				if (method == "balance" & length(which(!is.na(Metrics_df_ftest[[method]][software == soft & size == n]$mcc))) > 2){
					ttest <- t.test(Metrics_df_ftest[[method]][software == combo & size == n]$mcc, Metrics_df_ftest[[method]][software == soft & size == n]$mcc, alternative = "greater")
					dt[[method]][[as.character(n)]][row,"statistic"] <- ttest$statistic
					dt[[method]][[as.character(n)]][row,"p"] <- ttest$p.value

					ftest <- Bonett.Seier.test(x = Metrics_df_ftest[[method]][software == combo & size == n]$mcc, y = Metrics_df_ftest[[method]][software == soft & size == n]$mcc, alternative = c("two.sided"))
					dtf[[method]][[as.character(n)]][row,"statistic"] <- ftest$Statistic
					dtf[[method]][[as.character(n)]][row,"p"] <- ftest$p.value

				} else if (method == "remove_doublets"& length(which(!is.na(Metrics_df_ftest[[method]][software == soft & size == n]$ppv))) > 2){
					ttest <- t.test(Metrics_df_ftest[[method]][software == combo & size == n]$ppv, Metrics_df_ftest[[method]][software == soft & size == n]$ppv, alternative = "greater")
					dt[[method]][[as.character(n)]][row,"statistic"] <- ttest$statistic
					dt[[method]][[as.character(n)]][row,"p"] <- ttest$p.value

					ftest <- Bonett.Seier.test(Metrics_df_ftest[[method]][software == combo & size == n]$ppv, Metrics_df_ftest[[method]][software == soft & size == n]$ppv, alternative = c("two.sided"))
					dtf[[method]][[as.character(n)]][row,"statistic"] <- ftest$Statistic
					dtf[[method]][[as.character(n)]][row,"p"] <- ftest$p.value
				}
			}
			dt[[method]][[as.character(n)]]$fdr_p <- p.adjust(dt[[method]][[as.character(n)]]$p, method = "fdr")
			dtf[[method]][[as.character(n)]]$fdr_p <- p.adjust(dtf[[method]][[as.character(n)]]$p, method = "fdr")
		}
	}
	ttest_dt_list_all[[method]] <- rbindlist(dt[[method]])
	ftest_dt_list_all[[method]] <- rbindlist(dtf[[method]])
}



lapply(ttest_dt_list_all, function(x){x[statistic < 0]})
lapply(ttest_dt_list_all, as.data.frame)
lapply(ttest_dt_list_all, function(x){length((which(x$fdr_p < 0.05)))})
lapply(ttest_dt_list_all, function(x){length((which(!is.na(x$statistic))))})
lapply(ttest_dt_list_all, function(x){length((which(x$fdr_p < 0.05)))/length((which(!is.na(x$statistic))))})
lapply(ttest_dt_list_all_free, function(x){x[fdr_p > 0.05]})




pBalanced <- ggplot() +
	geom_bar(data = unique(Metrics_df_recommended[["balance"]][mt_percent == 0 & ambient_percent == 0 & subsampled == "normal", c("size", "software", "mcc_mean")]), stat = "identity", aes(x = as.factor(size), y = mcc_mean, color = software, fill = software), position = position_dodge(0.75), width=0.4) +
	geom_point(data = Metrics_df_recommended[["balance"]][mt_percent == 0 & ambient_percent == 0 & subsampled == "normal"], aes(x = size, y = mcc, group = software), size=0.5, position=position_dodge(0.75)) +
	theme_classic() +
	ylab("MCC") +
	xlab("Number Multiplexed Individuals") +
	scale_color_manual(values = c(software_colors)) +
	scale_fill_manual(values = c(software_colors)) +
	theme(legend.position = "none") +
	scale_y_continuous(limits=c(0.25,0.95), breaks = seq(0.25, 0.95, 0.1), oob = rescale_none)


save_figs(pBalanced,  paste0(outdir, "recommendations_balanced"), width =12, height = 7)



pSing <- ggplot() +
	geom_bar(data = unique(Metrics_df_recommended[["remove_doublets"]][mt_percent == 0 & ambient_percent == 0 & subsampled == "normal", c("size", "software", "ppv_mean")]), stat = "identity", aes(x = as.factor(size), y = ppv_mean, color = software, fill = software), position = position_dodge(0.75), width=0.4) +
	geom_point(data = Metrics_df_recommended[["remove_doublets"]][mt_percent == 0 & ambient_percent == 0 & subsampled == "normal"], aes(x = size, y = ppv, group = software), size=0.5, position=position_dodge(0.75)) +
	theme_classic() +
	ylab("PPV") +
	xlab("Number Multiplexed Individuals") +
	scale_color_manual(values = c(software_colors)) +
	scale_fill_manual(values = c(software_colors)) +
	theme(legend.position = "none") +
	scale_y_continuous(limits=c(0.8,1), breaks = seq(0.8, 1, 0.05), oob = rescale_none) 


save_figs(pSing,  paste0(outdir, "recommendations_strict"), width =6, height = 7)



##### Recommendations for reference-free #####
Metrics_df_combined_ref_free <- Metrics_df_combined[!grep("demuxlet", software)]


summary_w_single_ref_free  <- Metrics_df_combined_ref_free[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" ,.(ppv_mean = mean(ppv, na.rm = TRUE), mean_N = mean(true_singlet + false_singlet + true_doublet + false_doublet), min_N = min(true_singlet + false_singlet + true_doublet + false_doublet), max_N = max(true_singlet + false_singlet + true_doublet + false_doublet)),.(size,software,method)][size %in% c(2,4,8,16)]
summary_w_single_balance_ref_free  <- Metrics_df_combined_ref_free[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" ,.(balanced_accuracy_mean = mean(balanced_accuracy, na.rm = TRUE), mcc_mean = mean(mcc, na.rm = TRUE), mean_N = mean(true_singlet + false_singlet + true_doublet + false_doublet), min_N = min(true_singlet + false_singlet + true_doublet + false_doublet), max_N = max(true_singlet + false_singlet + true_doublet + false_doublet)),.(size,software,method)][size %in% c(2,4,8,16)]


summary_recommended_ref_free <- list()
summary_recommended_ref_free[["remove_doublets"]] <- unique(summary_w_single_ref_free[, .SD[ppv_mean %in% max(ppv_mean, na.rm = T)], by=size][,c("software", "method")])
summary_recommended_ref_free[["balance"]] <- unique(summary_w_single_balance_ref_free[, .SD[mcc_mean %in% max(mcc_mean, na.rm = T)], by=size][,c("software", "method")])



### Get 
Metrics_df_recommended_ref_free <- list()
Metrics_df_recommended_ref_free[["remove_doublets"]] <- Metrics_df_combined_ref_free[summary_recommended_ref_free[["remove_doublets"]], on = c("software", "method")][size %in% c(2,4,8,16)]
Metrics_df_recommended_ref_free[["balance"]] <- Metrics_df_combined_ref_free[summary_recommended_ref_free[["balance"]], on = c("software", "method")][size %in% c(2,4,8,16)]


for (software in unique(summary_recommended_ref_free[["remove_doublets"]]$software)){
	softwares <- strsplit(software,"-")[[1]]
	Metrics_df_recommended_ref_free[["remove_doublets"]] <- rbind(Metrics_df_recommended_ref_free[["remove_doublets"]], Metrics_df_combined_ref_free[software %in% softwares & size %in% c(2,4,8,16)])
}
Metrics_df_recommended_ref_free[["remove_doublets"]]  <- unique(Metrics_df_recommended_ref_free[["remove_doublets"]] )

for (software in unique(summary_recommended_ref_free[["balance"]]$software)){
	softwares <- strsplit(software,"-")[[1]]
	Metrics_df_recommended_ref_free[["balance"]] <- rbind(Metrics_df_recommended_ref_free[["balance"]], Metrics_df_combined_ref_free[software %in% softwares& size %in% c(2,4,8,16)])
}
Metrics_df_recommended_ref_free[["balance"]] <- unique(Metrics_df_recommended_ref_free[["balance"]])


summary_w_single_balance <- Metrics_df_combined_ref_free[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal",.(balanced_accuracy_mean = mean(balanced_accuracy, na.rm = TRUE), mcc_mean = mean(mcc, na.rm = TRUE), mean_N = mean(true_singlet + false_singlet + true_doublet + false_doublet), min_N = min(true_singlet + false_singlet + true_doublet + false_doublet), max_N = max(true_singlet + false_singlet + true_doublet + false_doublet)),.(size,software,method)][size %in% c(2,4,8,16)]
summary_w_single <- Metrics_df_combined_ref_free[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal",.(ppv_mean = mean(ppv, na.rm = TRUE), mean_N = mean(true_singlet + false_singlet + true_doublet + false_doublet), min_N = min(true_singlet + false_singlet + true_doublet + false_doublet), max_N = max(true_singlet + false_singlet + true_doublet + false_doublet)),.(size,software,method)][size %in% c(2,4,8,16)]

Metrics_df_recommended_ref_free[["remove_doublets"]] <- summary_w_single[Metrics_df_recommended_ref_free[["remove_doublets"]], on = c("size", "software", "method")]
Metrics_df_recommended_ref_free[["balance"]] <- summary_w_single_balance[Metrics_df_recommended_ref_free[["balance"]], on = c("size", "software", "method")]


Metrics_df_recommended_ref_free[["balance"]]$software <- factor(Metrics_df_recommended_ref_free[["balance"]]$software, levels = unique(Metrics_df_recommended_ref_free[["balance"]][order(Metrics_df_recommended_ref_free[["balance"]]$mcc_mean)]$software))
Metrics_df_recommended_ref_free[["remove_doublets"]]$software <- factor(Metrics_df_recommended_ref_free[["remove_doublets"]]$software, unique(Metrics_df_recommended_ref_free[["remove_doublets"]][order(Metrics_df_recommended_ref_free[["remove_doublets"]]$ppv_mean)]$software))




#### Statistically test between groups ####
ttest_free_dt_list <- list()

for (method in names(Metrics_df_recommended_ref_free)){
	dt2 <- list()
	for (combo in unique(grep("-", Metrics_df_recommended_ref_free[[method]]$software, value = TRUE))){
		print(combo)
		softwares <- strsplit(combo,"-")[[1]]
		dt <- list()
		for (n in unique(Metrics_df_recommended_ref_free[[method]][software == combo]$size)){
			dt[[method]][[as.character(n)]] <- data.table(software1 = rep(combo, length(softwares)), software2 = softwares, size = n)
			dt[[method]][[as.character(n)]]$statistic <- 0
			dt[[method]][[as.character(n)]]$p <- 0

			for (row in 1:nrow(dt[[method]][[as.character(n)]])){
				soft <- dt[[method]][[as.character(n)]]$software2[row]
				if (method == "balance"){
					ttest <- t.test(Metrics_df_recommended_ref_free[[method]][software == combo & size == n]$mcc, Metrics_df_recommended_ref_free[[method]][software == soft & size == n]$mcc, alternative = "greater")
					dt[[method]][[as.character(n)]][row,"statistic"] <- ttest$statistic
					dt[[method]][[as.character(n)]][row,"p"] <- ttest$p.value

				} else if (method == "remove_doublets"){
					ttest <- t.test(Metrics_df_recommended_ref_free[[method]][software == combo & size == n]$ppv, Metrics_df_recommended_ref_free[[method]][software == soft & size == n]$ppv, alternative = "greater")
					dt[[method]][[as.character(n)]][row,"statistic"] <- ttest$statistic
					dt[[method]][[as.character(n)]][row,"p"] <- ttest$p.value
				}
			}
		}
	dt2[[method]][[combo]] <- rbindlist(dt[[method]])
	}
	ttest_free_dt_list[[method]] <- rbindlist(dt2[[method]])
}



##### Ftest for all softwares vs combos #####
summary_recommended_ref_free_n <- list()
summary_recommended_ref_free_n[["remove_doublets"]] <- unique(summary_w_single_ref_free[, .SD[ppv_mean %in% max(ppv_mean, na.rm = T)], by=size][,c("software", "method", "size")])
summary_recommended_ref_free_n[["balance"]] <- unique(summary_w_single_balance_ref_free[, .SD[mcc_mean %in% max(mcc_mean, na.rm = T)], by=size][,c("software", "method", "size")])

Metrics_df_ftest_free <- list()
Metrics_df_ftest_free[["remove_doublets"]] <- rbind(Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & method == "single_software" & software %in% c(demultiplexing_list[2:5], doublet_detection_list)], Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" ][summary_recommended_ref_free_n[["remove_doublets"]], on = c("size", "software", "method")])
Metrics_df_ftest_free[["balance"]] <- rbind(Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & method == "single_software" & software %in% c(demultiplexing_list[2:5], doublet_detection_list)], Metrics_df_combined[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" ][summary_recommended_ref_free_n[["balance"]], on = c("size", "software", "method")])



ttest_dt_list_all_free <- list()
ftest_dt_list_all_free <- list()

dtf_free <- list()
dt_free <- list()

dt2_free <- list()
dt2f_free <- list()



for (method in names(Metrics_df_ftest_free)){
	for (combo in unique(grep("-", Metrics_df_ftest_free[[method]]$software, value = TRUE))){
		print(combo)
		for (n in unique(Metrics_df_ftest_free[[method]][software == combo]$size)){
			softwares <- unique(Metrics_df_ftest_free[[method]][method == "single_software"]$software)
			dtf_free[[method]][[as.character(n)]] <- data.table(software1 = rep(combo, length(softwares)), software2 = softwares, size = n)
			dtf_free[[method]][[as.character(n)]]$statistic <- as.numeric(NA)
			dtf_free[[method]][[as.character(n)]]$p <- as.numeric(NA)

			dt_free[[method]][[as.character(n)]] <- dtf_free[[method]][[as.character(n)]]

			for (row in 1:nrow(dtf_free[[method]][[as.character(n)]])){
				soft <- dtf_free[[method]][[as.character(n)]]$software2[row]
				print(soft)
				if (method == "balance" & length(which(!is.na(Metrics_df_ftest_free[[method]][software == soft & size == n]$mcc))) > 2){
					ttest <- t.test(Metrics_df_ftest_free[[method]][software == combo & size == n]$mcc, Metrics_df_ftest_free[[method]][software == soft & size == n]$mcc, alternative = "greater")
					dt_free[[method]][[as.character(n)]][row,"statistic"] <- ttest$statistic
					dt_free[[method]][[as.character(n)]][row,"p"] <- ttest$p.value

					ftest <- Bonett.Seier.test(x = Metrics_df_ftest_free[[method]][software == combo & size == n]$mcc, y = Metrics_df_ftest_free[[method]][software == soft & size == n]$mcc, alternative = c("two.sided"))
					dtf_free[[method]][[as.character(n)]][row,"statistic"] <- ftest$Statistic
					dtf_free[[method]][[as.character(n)]][row,"p"] <- ftest$p.value

				} else if (method == "remove_doublets"& length(which(!is.na(Metrics_df_ftest_free[[method]][software == soft & size == n]$ppv))) > 2){
					ttest <- t.test(Metrics_df_ftest_free[[method]][software == combo & size == n]$ppv, Metrics_df_ftest_free[[method]][software == soft & size == n]$ppv, alternative = "greater")
					dt_free[[method]][[as.character(n)]][row,"statistic"] <- ttest$statistic
					dt_free[[method]][[as.character(n)]][row,"p"] <- ttest$p.value

					ftest <- Bonett.Seier.test(Metrics_df_ftest_free[[method]][software == combo & size == n]$ppv, Metrics_df_ftest_free[[method]][software == soft & size == n]$ppv, alternative = c("two.sided"))
					dtf_free[[method]][[as.character(n)]][row,"statistic"] <- ftest$Statistic
					dtf_free[[method]][[as.character(n)]][row,"p"] <- ftest$p.value
				}
			}
			dt_free[[method]][[as.character(n)]]$fdr_p <- p.adjust(dt_free[[method]][[as.character(n)]]$p, method = "fdr")
			dtf_free[[method]][[as.character(n)]]$fdr_p <- p.adjust(dtf_free[[method]][[as.character(n)]]$p, method = "fdr")
		}
	}
	ttest_dt_list_all_free[[method]] <- rbindlist(dt_free[[method]])
	ftest_dt_list_all_free[[method]] <- rbindlist(dtf_free[[method]])
}


lapply(ttest_dt_list_all_free, function(x){x[statistic < 0]})
lapply(ttest_dt_list_all_free, as.data.frame)
lapply(ttest_dt_list_all_free, function(x){length((which(x$fdr_p < 0.05)))})
lapply(ttest_dt_list_all_free, function(x){x[fdr_p > 0.05]})
lapply(ttest_dt_list_all_free, function(x){length((which(!is.na(x$statistic))))})
lapply(ttest_dt_list_all_free, function(x){length((which(x$fdr_p < 0.05)))/length((which(!is.na(x$statistic))))})



pBalanced_ref_free <- ggplot() +
	geom_bar(data = unique(Metrics_df_recommended_ref_free[["balance"]][mt_percent == 0 & ambient_percent == 0 & subsampled == "normal", c("size", "software", "mcc_mean")]), stat = "identity", aes(x = as.factor(size), y = mcc_mean, color = software, fill = software), position = position_dodge(0.75), width=0.4) +
	geom_point(data = Metrics_df_recommended_ref_free[["balance"]][mt_percent == 0 & ambient_percent == 0 & subsampled == "normal"], aes(x = size, y = mcc, group = software), size=0.5, position=position_dodge(0.75)) +
	theme_classic() +
	ylab("MCC") +
	xlab("Number Multiplexed Individuals") +
	scale_color_manual(values = c(software_colors)) +
	scale_fill_manual(values = c(software_colors)) +
	theme(legend.position = "none") +
	scale_y_continuous(limits=c(0.25,0.95), breaks = seq(0.25, 0.95, 0.1), oob = rescale_none)


save_figs(pBalanced_ref_free,  paste0(outdir, "recommendations_balanced_ref_free"), width =12, height = 7)


