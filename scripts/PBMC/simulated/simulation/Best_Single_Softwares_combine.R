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
library(ggpubr)



##### Set up Functions #####
save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}

##### Set up colors for plotting #####
# downsample_colors <- c("#C7EAE5", "#35978F")
downsample_colors <- c("white", "#35978F")
mt_percent_colors <- c("0" = "white", "5" = "#D0E0D3", "10" = "#A3C4A6", "25" = "#1E6D2D")
ambient_percent_colors <- c("0" = "white", "10" = "#FDEDD2", "20" = "#FBD491", "50" = "#E6AB00")
uneven_n_colors <- c("0" = "white", "50" = "#E986AE", "75" = "#C74F7D", "95" = "#A11C44")


##### Set up universal 
demultiplexing_list <- c("demuxalot", "demuxalot_refined","demuxlet", "dropulation", "freemuxlet","scSplit","souporcell","vireo")
names(demultiplexing_list) <- c("Demuxalot", "Demuxalot (refined)","Demuxlet", "Dropulation", "Freemuxlet","ScSplit","Souporcell","Vireo")
doublet_detection_list <- c("DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder", "scDblFinder_known_doublets", "scds",  "scrublet", "solo")
names(doublet_detection_list) <- c("DoubletDecon", "DoubletDetection","DoubletFinder",  "ScDblFinder", "ScDblFinder (known doublets)", "Scds", "Scrublet", "Solo")
softwares <- c(demultiplexing_list, doublet_detection_list)

softwares_dt <- data.table(software = softwares, Software = names(softwares))


software_colors <- c("#575E57", "#008751", "#5DBB4D", "#EFE305", "#F8BE33", "#F97919", "#D21F09", "#8D132B", "#F980BF", "#EAA1D1", "#BAAAFF", "#9B64E0", "#6B82E3", "#63D4FE", "#39BCBC", "#C4E4DF")
names(software_colors) <- names(softwares)


# ###### Set up Directories #####
dir <- "/path/to/output/PBMC/simulation/SimulatedOverlap/"


##### Get Pools #####
pools <- dir(dir, pattern = "size")
pools_sub <- c()
for (pool in pools){
    if (file.exists(paste0(dir,pool, "/percent_correct_per_barcode_single_soft.rds"))){
        pools_sub <- c(pools_sub, pool)
    }
}


##### Read in Data #####
Metrics <- lapply(pools_sub, function(pool){
	fread(paste0(dir, pool,"/single_software_metrics.tsv"), sep = "\t")
})

Metrics_df <- do.call(rbind, Metrics)


fwrite(Metrics_df, paste0(dir, "single_software_metrics.tsv"))


### normal
Metrics_df_nodiff <- Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]
cols <- colnames(Metrics_df_nodiff)[!(colnames(Metrics_df_nodiff) %in% c("mt_percent", "ambient_percent", "subsampled", "uneven"))]
Metrics_df_nodiff <- Metrics_df_nodiff[,..cols]

fwrite(Metrics_df_nodiff, paste0(dir, "single_software_metrics_no_diff.tsv"), sep = "\t")


### ambient
Metrics_df_ambient <- Metrics_df[(mt_percent == 0 & subsampled == "normal" & uneven == 0) & (grepl("size\\d+_SimulatedPool[1,2,3]", pool) | ambient_percent > 0 )]
cols_amb <- colnames(Metrics_df_ambient)[!(colnames(Metrics_df_ambient) %in% c("mt_percent", "subsampled", "uneven"))]
Metrics_df_ambient <- Metrics_df_ambient[,..cols_amb]

fwrite(Metrics_df_ambient, paste0(dir, "single_software_metrics_ambient.tsv"), sep = "\t")


### mt
Metrics_df_mt <- Metrics_df[(ambient_percent == 0 & subsampled == "normal" & uneven == 0) & (grepl("size\\d+_SimulatedPool[1,2,3]", pool) | mt_percent > 0 )]
cols_mt <- colnames(Metrics_df_mt)[!(colnames(Metrics_df_mt) %in% c("ambient_percent", "subsampled", "uneven"))]
Metrics_df_mt <- Metrics_df_mt[,..cols_mt]

fwrite(Metrics_df_mt, paste0(dir, "single_software_metrics_mt.tsv"), sep = "\t")


### subsampled
Metrics_df_subsample <- Metrics_df[(ambient_percent == 0 & mt_percent == 0 & uneven == 0) & (grepl("size\\d+_SimulatedPool[1,2,3]", pool) | subsampled == "normal" )]
cols_subsample <- colnames(Metrics_df_subsample)[!(colnames(Metrics_df_subsample) %in% c("ambient_percent", "mt_percent", "uneven"))]
Metrics_df_subsample <- Metrics_df_subsample[,..cols_subsample]

fwrite(Metrics_df_subsample, paste0(dir, "single_software_metrics_subsample.tsv"), sep = "\t")


### uneven
Metrics_df_uneven <- Metrics_df[(ambient_percent == 0 & mt_percent == 0 & subsampled == "normal") & (grepl("size\\d+_SimulatedPool[1,2,3]", pool) | uneven > 0 )]
cols_uneven <- colnames(Metrics_df_uneven)[!(colnames(Metrics_df_uneven) %in% c("ambient_percent", "mt_percent", "subsampled"))]
Metrics_df_uneven <- Metrics_df_uneven[,..cols_uneven]

fwrite(Metrics_df_uneven, paste0(dir, "single_software_metrics_uneven.tsv"), sep = "\t")





setkey(Metrics_df, mcc)

summary <- Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0,.(MCC_Mean=mean(mcc, na.rm = TRUE), Balanced_Accuracy_Mean = mean(balanced_accuracy, na.rm = TRUE)),.(size,software)];

data.frame(summary[order(size, MCC_Mean)])

summary2 <- Metrics_df[,.(MCC_Mean=mean(mcc, na.rm = TRUE), Balanced_Accuracy_Mean = mean(balanced_accuracy, na.rm = TRUE)),.(size,software,mt_percent,ambient_percent,subsampled,uneven)];


data.frame(summary2[order(size,mt_percent,ambient_percent,subsampled,uneven,MCC_Mean)][size == 2])


Metrics_df[, max(mcc, na.rm = T), by = size]
Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0, .SD[mcc %in% mean(mcc, na.rm = T)], by=size]


##### Compare MCC to balanced accuracy #####
cor(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]$mcc, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]$balanced_accuracy, use = "complete", method = "spearman")
mcc_balanced_acc_cor <- cor.test(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]$mcc, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]$balanced_accuracy, use = "complete")
pt(q = as.numeric(mcc_balanced_acc_cor$statistic), df = 1075, lower.tail=FALSE)*2\


mcc_balanced_acc_spearman <- cor.test(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]$mcc, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]$balanced_accuracy, use = "complete", method = "spearman", ,exact=FALSE)
pt(q = as.numeric(mcc_balanced_acc_spearman$statistic), df = 1075, lower.tail=FALSE)*2


##### Make of MCC and Balanced Accuracy for all softwares across individuals #####
Metrics_df$group <- ifelse(Metrics_df$software %in% demultiplexing_list, "Demultiplexing\nMethods", "Doublet Detecting\nMethods")

pMCC_facet <- ggplot(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal"  & uneven == 0], aes(as.factor(size), mcc, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_grid(~ group) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("MCC") +
	xlab("Number Individuals Multiplexed") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6)

save_figs(pMCC_facet,  paste0(dir, "mcc_single_software_facet"), width =15, height = 6.3)

pBalanced_Accuracy_facet <- ggplot(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal"  & uneven == 0], aes(as.factor(size), balanced_accuracy, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_grid(~ group) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("Balanced Accuracy") +
	xlab("Number Individuals Multiplexed") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6)

save_figs(pBalanced_Accuracy_facet,  paste0(dir, "balanced_accuracy_single_software_facet"), width =15, height = 6.3)


### Statistical tests to include in the paper ###
cortest <- cor.test(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]$mcc, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0]$balanced_accuracy, method = "spearman", na.rm = TRUE)
cortest$p.value

for (soft in unique(Metrics_df$software)){
	print(soft)
	print(cor.test(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & software == soft & uneven == 0]$mcc, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & software == soft & uneven == 0]$balanced_accuracy, method = "spearman", na.rm = TRUE))
}

### t ttest at each time comparing demultiplexign to doublet detecting softwares ###
ttest <- data.table(unique(Metrics_df[,size]))
colnames(ttest) <- c("size")
ttest$mcc_p <- as.numeric(NA)
ttest$balanced_accuracy_p <- as.numeric(NA)

for (row in 1:nrow(ttest)){
	n <- ttest$size[row]
	ttest1 <- t.test(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal"& uneven == 0  & software %in% demultiplexing_list & size == n & software != "scSplit"]$mcc, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0 & size == n & software %in% doublet_detection_list]$mcc, na.rm = TRUE)
	ttest$mcc_p[row] <- ttest1$p.value
	ttest2 <- t.test(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0 & software %in% demultiplexing_list & size == n & software != "scSplit"]$balanced_accuracy, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0 & size == n & software %in% doublet_detection_list]$balanced_accuracy, na.rm = TRUE)
	ttest$balanced_accuracy_p[row] <- ttest2$p.value
}

### Test that souporcell goes down for 64 and 128 individuals compared to 32 ###
ttest_souporcell <- data.table(size1 = c(32, 32), size2 = c(64, 128))
ttest_souporcell$mcc_p <- as.numeric(NA)
ttest_souporcell$balanced_accuracy_p <- as.numeric(NA)

for (row in 1:nrow(ttest_souporcell)){
	n1 <- ttest_souporcell$size1[row]
	n2 <- ttest_souporcell$size2[row]
	ttest1 <- t.test(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0  & size == n1 & software == "souporcell"]$mcc, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & size == n2 & uneven == 0 & software == "souporcell"]$mcc, na.rm = TRUE)
	ttest_souporcell$mcc_p[row] <- ttest1$p.value
	ttest2 <- t.test(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0 & size == n1 & software == "souporcell"]$balanced_accuracy, Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & size == n2 & uneven == 0 & software == "souporcell"]$balanced_accuracy, na.rm = TRUE)
	ttest_souporcell$balanced_accuracy_p[row] <- ttest2$p.value
}


### Test for difference in variance for smaller vs larger pool sizes ###
var_df <- data.table(unique(Metrics_df[,size]))
colnames(var_df) <- c("size")
var_df$var_mcc <- as.numeric(NA)
var_df$var_balanced_accuracy <- as.numeric(NA)

for (row in 1:nrow(var_df)){
	n <- var_df$size[row]
	var_df$var_mcc[row] <- var(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & size == n & uneven == 0 & software %in% c("scds", "solo", "scDblFinder", "DoubletDetection")]$mcc, na.rm = TRUE)
	var_df$var_balanced_accuracy[row] <- var(Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & size == n & uneven == 0 & software %in% c("scds", "solo", "scDblFinder", "DoubletDetection")]$balanced_accuracy, na.rm = TRUE)
}

cor.test(var_df$size, var_df$var_mcc, method = "spearman")
cor.test(var_df$size, var_df$var_balanced_accuracy, method = "spearman")


### test for scdbl finder consistency with and without known doublets ###
dt4scdblfinder_cor <- Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0 & software  == "scDblFinder", c("size","pool","software", "mcc")][Metrics_df[mt_percent == 0 & ambient_percent == 0 & subsampled == "normal" & uneven == 0 & software  == "scDblFinder_known_doublets", c("size", "pool","software", "mcc")], on = c("size", "pool")]

scdblfinder_cor <- cor.test(dt4scdblfinder_cor$mcc, dt4scdblfinder_cor$i.mcc, use = "complete.obs")

pt(q = as.numeric(scdblfinder_cor$statistic), df = 66, lower.tail=FALSE)*2

scdblfinder_ttest_df <- data.table(unique(Metrics_df[,size]))
colnames(scdblfinder_ttest_df) <- "size"
scdblfinder_ttest_df$p <- as.numeric(NA)
scdblfinder_ttest_df$estimate_w_knowndoublets <- as.numeric(NA)
scdblfinder_ttest_df$estimate_NO_knowndoublets <- as.numeric(NA)

for (row in 1:nrow(scdblfinder_ttest_df)){
	n <- scdblfinder_ttest_df$size[row]
	ttest <- t.test(dt4scdblfinder_cor[size == n]$mcc, dt4scdblfinder_cor[size == n]$i.mcc, na.rm = TRUE)
	scdblfinder_ttest_df$p[row] <- ttest$p.value
	scdblfinder_ttest_df$estimate_w_knowndoublets[row] <- ttest$estimate[1]
	scdblfinder_ttest_df$estimate_w_knowndoublets[row] <- ttest$estimate[2]
}


##### Make of Balanced Accuracy for all softwares with difficult pools across individuals with demultiplexing and doublet detecting faceted #####
pMCC_amb_facet <- ggplot(Metrics_df[ambient_percent > 0 | (grepl("Pool[1-3]$", pool))], aes(ambient_percent*100, mcc, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25) +
	facet_grid(group ~ size) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("MCC") +
	xlab("Additional Ambient Percent") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	scale_x_continuous(breaks=c(0,25,50)) +
	theme(legend.position = "none")

save_figs(pMCC_amb_facet,  paste0(dir, "mcc_single_software_ambient_facet"), width =11, height = 7)



pMCC_sub_facet <- ggplot(Metrics_df[subsampled != "normal" | (grepl("Pool[1-3]$", pool))], aes(subsampled, mcc, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25) +
	facet_grid(group ~ size) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("MCC") +
	xlab("Reads Downsampled") +
	geom_smooth(aes(group = software), se = FALSE, method = lm, size = 0.6) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) +
	theme(legend.position = "none")

save_figs(pMCC_sub_facet,  paste0(dir, "mcc_single_software_downsampled_facet"), width =11, height = 8)


pMCC_mt_facet <- ggplot(Metrics_df[mt_percent > 0 | (grepl("Pool[1-3]$", pool))], aes(mt_percent*100, mcc, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25) +
	facet_grid(group ~ size) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("MCC") +
	xlab("Additional Mitochondiral Percent") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
	scale_x_continuous(breaks=c(0,25,50)) +
	theme(legend.position = "none")

save_figs(pMCC_mt_facet,  paste0(dir, "mcc_single_software_mt_percent_facet"), width =11, height = 7)



pMCC_uneven_facet <- ggplot(Metrics_df[(uneven > 0 & uneven < 0.8) | (grepl("Pool[1-3]$", pool))], aes(uneven, mcc, color = Software)) +
	geom_jitter(size = 0.25, width = 0.025) +
	facet_grid(group ~ size) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("MCC") +
	xlab("Spiked Donor - Proportion of Pool") +
	scale_x_continuous(breaks=c(0,0.25,0.5, 0.75)) +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
	theme(legend.position = "none")

save_figs(pMCC_uneven_facet,  paste0(dir, "mcc_single_software_uneven_facet"), width =11, height = 7)


##### Make of Balanced Accuracy for all softwares with difficult pools across individuals #####
pBalanced_Accuracy_amb_facet <- ggplot(Metrics_df[ambient_percent > 0 | (grepl("Pool[1-3]$", pool))], aes(ambient_percent*100, balanced_accuracy, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25) +
	facet_grid(group ~ size) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("Balanced Accuracy") +
	xlab("Additional Ambient Percent") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	scale_x_continuous(breaks=c(0,25,50)) +
	theme(legend.position = "none")

save_figs(pBalanced_Accuracy_amb_facet,  paste0(dir, "balanced_accuracy_single_software_ambient_facet"), width =11, height = 7)



pBalanced_Accuracy_sub_facet <- ggplot(Metrics_df[subsampled != "normal" | (grepl("Pool[1-3]$", pool))], aes(subsampled, balanced_accuracy, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25) +
	facet_grid(group ~ size) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("Balanced Accuracy") +
	xlab("Reads Downsampled") +
	geom_smooth(aes(group = software), se = FALSE, method = lm, size = 0.6) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) +
	theme(legend.position = "none")

save_figs(pBalanced_Accuracy_sub_facet,  paste0(dir, "balanced_accuracy_single_software_downsampled_facet"), width =11, height = 8)



pBalanced_Accuracy_mt_facet <- ggplot(Metrics_df[mt_percent > 0 | (grepl("Pool[1-3]$", pool))], aes(mt_percent*100, balanced_accuracy, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25) +
	facet_grid(group ~ size) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("Balanced Accuracy") +
	xlab("Additional Mitochondiral Percent") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
	scale_x_continuous(breaks=c(0,25,50)) +
	theme(legend.position = "none")

save_figs(pBalanced_Accuracy_mt_facet,  paste0(dir, "balanced_accuracy_single_software_mt_percent_facet"), width =11, height = 7)


pBalanced_Accuracy_uneven_facet <- ggplot(Metrics_df[ (uneven > 0 & uneven < 0.8) | (grepl("Pool[1-3]$", pool))], aes(uneven, balanced_accuracy, color = Software)) +
	geom_jitter(size = 0.25, width = 0.025) +
	facet_grid(group ~ size) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("Balanced Accuracy") +
	xlab("Spiked Donor - Proportion of Pool") +
	scale_x_continuous(breaks=c(0,0.25,0.5, 0.75)) +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
	theme(legend.position = "none")

save_figs(pBalanced_Accuracy_uneven_facet,  paste0(dir, "balanced_accuracy_single_software_uneven_donors_facet"), width =11, height = 7)


