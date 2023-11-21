library("tidyr")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("ComplexHeatmap")
library(awtools)
library(circlize)
library(ggpubr)
library(viridis)
library(ggnewscale)
library(ggforce)
library(ComplexUpset)
library(RColorBrewer)
library(plyr)
library(data.table)
library(facetscales)
library(ggh4x)
library(colorspace)
library(future.apply)

set.seed(79)

##### Set directories #####
dir <- "/path/to/output/PBMC/" ### the outdir used in the Snakefile
out <- "/path/to/output/PBMC/Round1Overlap/" ### path to write files to, including files to use for pool simulation
outdir <- "/path/to/output/PBMC/Round1Overlap/Figures/"
meta <- fread("/path/to/Demuxafy_manuscript/files/PBMC/PBMC_sample_meta.tsv", sep = "\t")
dir.create(out)
sample_list <- dir(path = dir, pattern = "scFibroblast_EQTL_Sample")

softwares <- c("demuxalot", "demuxalot_refined","demuxlet", "dropulation", "freemuxlet","scSplit","souporcell","vireo", "DoubletDecon", "DoubletDetection","DoubletFinder",  "scDblFinder", "scds", "scrublet", "solo")
names(softwares) <- c("Demuxalot", "Demuxalot (refined)","Demuxlet", "Dropulation", "Freemuxlet","ScSplit","Souporcell","Vireo", "DoubletDecon", "DoubletDetection","DoubletFinder",  "ScDblFinder", "Scds", "Scrublet", "Solo")
demultiplexing_softwares <- c("demuxalot", "demuxalot_refined","demuxlet", "dropulation", "freemuxlet","scSplit","souporcell","vireo")
doublet_detection_softwares <- c("DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder","scds",  "scrublet", "solo")

software_names <- c("Demuxalot", "Demuxalot (refined)","Demuxlet", "Dropulation", "Freemuxlet","ScSplit","Souporcell","Vireo", "DoubletDecon", "DoubletDetection","DoubletFinder",  "ScDblFinder", "Scds", "Scrublet", "Solo")
demultiplexing_names <- c("Demuxalot", "Demuxalot (refined)","Demuxlet", "Dropulation", "Freemuxlet","ScSplit","Souporcell","Vireo")
doubletdetecting_names <- c("DoubletDecon", "DoubletDetection","DoubletFinder",  "ScDblFinder", "Scds", "Scrublet", "Solo")


colors <- c("#575E57", "#008751", "#5DBB4D", "#EFE305", "#F8BE33", "#F97919", "#D21F09", "#8D132B", "#F980BF", "#EAA1D1", "#BAAAFF", "#9B64E0", "#6B82E3", "#63D4FE", "#39BCBC", "#C4E4DF", "#E2E2E2")
names(colors) <- c("Demuxalot","Demuxalot (refined)","Demuxlet","Dropulation","Freemuxlet","ScSplit","Souporcell","Vireo","DoubletDecon", "DoubletDetection","DoubletFinder","ScDblFinder","ScDblFinder (known doublets)","Scds","Scrublet","Solo","Solo (known doublets)")


result_files <- paste0(dir,sample_list,"/CombinedResults/CombinedDropletAssignments.tsv")
names(result_files) <- sample_list


##### Functions #####
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





##### Read in Files #####
results_list <- lapply(result_files, function(x){
  fread(x, sep = "\t")
})
names(results_list) <- sample_list

##### Read in the key tables for freemuxlet, scSplit and souporcell #####
key <- readRDS(paste0(out,"PoolKeys.rds"))
key <- lapply(key, function(x){
  x$Cluster_ID <- gsub("CLUST","",x$Cluster_ID)
  x$Software <- paste0(x$Software, "_Assignment")
  return(data.table(x))
})

##### Pivot the dataframes longer for joining#####
results_list_long <- lapply(results_list, function(x){
  melt(x, measure.vars = c("demuxalot_Assignment", "demuxalot_refined_Assignment","demuxlet_Assignment", "dropulation_Assignment", "freemuxlet_Assignment","scSplit_Assignment","souporcell_Assignment","vireo_Assignment"), variable.name = "Software", value.name = "Assignment")
})

##### Left_join the common assignments to the dataframe #####
results_list_long <- lapply(names(key), function(x){
  temp <- key[[x]][results_list_long[[x]], on = c("Cluster_ID" = "Assignment", "Software" = "Software")]
  colnames(temp) <- gsub("Cluster_ID", "Assignment", colnames(temp))
  # temp <- left_join(results_list_long[[x]],key[[x]], by = c("Assignment" = "Cluster_ID", "Software" = "Software"))
  temp$Genotype_ID <- ifelse(temp$Software == "demuxalot_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Software == "demuxalot_refined_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Software == "demuxlet_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Software == "dropulation_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Software == "vireo_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Assignment == "doublet", "doublet", temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(is.na(temp$Genotype_ID), "unassigned", temp$Genotype_ID)
  colnames(temp) <- gsub("V42", "solo_DropletScore", colnames(temp))
  temp$Pool <- x
  return(temp)
})
names(results_list_long) <- names(key)

any(lapply(results_list_long, function(x) any(is.na(x$Genotype_ID))))


##### Pivot wider to get the counts of singlet, doublets... and intersection across softwares #####
results_list_wide <- lapply(results_list_long, function(x){
  x$Assignment <- NULL
  x$Correlation <- NULL
  x <- pivot_wider(x,names_from = "Software", values_from = "Genotype_ID")
  return(data.table(x))
})
names(results_list_wide) <- names(key)


##### Add information on number of singlet and doublet calls #####
results_list_wide <- lapply(results_list_wide, function(x){
  columns <- grep("_DropletType", colnames(x), value = TRUE)
  demultiplexing_columns <- paste0(demultiplexing_softwares, "_DropletType")
  doublet_detecting_columns <- paste0(doublet_detection_softwares, "_DropletType")
  x[,"Singlets"] <- rowSums(x[,..columns] == "singlet", na.rm = TRUE)
  x[,"Doublets"] <- rowSums(x[,..columns] == "doublet", na.rm = TRUE)
  x[,"Demultiplexing_Singlets"]<- rowSums(x[, ..demultiplexing_columns] == "singlet", na.rm = TRUE)
  x[,"Demultiplexing_Doublets"]<- rowSums(x[, ..demultiplexing_columns] == "doublet", na.rm = TRUE)  
  x[,"DoubletDetecting_Singlets"]<- rowSums(x[, ..doublet_detecting_columns] == "singlet", na.rm = TRUE)
  x[,"DoubletDetecting_Doublets"]<- rowSums(x[, ..doublet_detecting_columns] == "doublet", na.rm = TRUE)
  return(x)
})


results_wide <- do.call(rbind, results_list_wide)



##### Bring in metadata - produced by MakeSeuratMetaData.R #####
QC_df <- fread(paste0(out, "Onek1k_metadata.tsv"), sep = "\t")
QC_df_updated <- QC_df[,c("Pool", "Barcode", "percent.mt", "percent.rb", "nCount_RNA", "nFeature_RNA")]

QC_df_updated <- QC_df_updated[results_wide, on = c("Barcode", "Pool")]



QC_metrics_long <- pivot_longer(QC_df_updated[,c("Pool", "percent.mt", "nCount_RNA", "nFeature_RNA", "Singlets", "Doublets", "Demultiplexing_Singlets", "Demultiplexing_Doublets", "DoubletDetecting_Singlets", "DoubletDetecting_Doublets")],
                      cols = c("percent.mt", "nCount_RNA", "nFeature_RNA"), names_to = "QC_Metric", values_to = "QC_Metric_Value")

N_dt <- data.table(table(QC_metrics_long$Singlets))
colnames(N_dt) <- c("Singlets", "N")
N_dt$QC_Metric <- "nCount_RNA"


scale_y <- list(
  percent.mt = scale_y_continuous(limits = c(0, 100)),
  nCount_RNA = scale_y_continuous(limits = c(0, 110000)),
  nFeature_RNA = scale_y_continuous(limits = c(0, 7000))
)

p_singlets_qc <- ggplot(QC_metrics_long, aes(factor(Singlets, levels = c(seq(0,15), "all")), QC_Metric_Value)) +
  geom_violin() +
  stat_summary(, fun = "mean", size = 0.5,
                 geom="point", color="firebrick3") +
  theme_bw() +
  facet_grid_sc(QC_Metric ~ ., scales = list(y=scale_y), switch = "y", labeller = as_labeller(c(percent.mt = "Mitochondrial\nPercent", nCount_RNA = "Number\nUMIs", nFeature_RNA = "Number\nGenes")))  +
  ylab(NULL) +
  xlab("Number Softwares Identified\nDroplet as Singlet") +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  force_panelsizes(rows = c(1.5, 1, 1)) +
  geom_text(data = N_dt, aes(y = 85000, label=N, hjust=0,vjust=0,angle=45))

ggsave(p_singlets_qc, filename = paste0(outdir,"Singlets_v_QC.pdf"), width = 10, height = 8, units = "cm")
ggsave(p_singlets_qc, filename = paste0(outdir,"Singlets_v_QC.png"), width = 10, height = 8, units = "cm")


N_doub_dt <- data.table(table(QC_metrics_long$Doublets))
colnames(N_doub_dt) <- c("Doublets", "N")
N_doub_dt$QC_Metric <- "nCount_RNA"

p_doublets_qc <- ggplot(QC_metrics_long, aes(factor(Doublets, levels = c(seq(0,15), "all")), QC_Metric_Value)) +
  geom_violin() +
  stat_summary(fun = "mean", size = 0.5,
                 geom="point", color="firebrick3") +
  theme_bw() +
  facet_grid_sc(QC_Metric ~ ., scales = list(y=scale_y), switch = "y", labeller = as_labeller(c(percent.mt = "Mitochondrial Percent", nCount_RNA = "Number UMIs", nFeature_RNA = "Number Genes")))  +
  ylab(NULL) +
  xlab("Number of Softwares Identified Droplet as Doublets") +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  force_panelsizes(rows = c(1.5, 1, 1)) +
  geom_text(data = N_doub_dt, aes(y = 85000, label=N, hjust=0,vjust=0,angle=45))

ggsave(p_doublets_qc, filename = paste0(outdir,"Doublets_v_QC.pdf"), width = 17, height = 14, units = "cm")
ggsave(p_doublets_qc, filename = paste0(outdir,"Doublets_v_QC.png"), width = 17, height = 14, units = "cm")



N_demux_dt <- data.table(table(QC_metrics_long$Demultiplexing_Singlets))
colnames(N_demux_dt) <- c("Demultiplexing_Singlets", "N")
N_demux_dt$QC_Metric <- "nCount_RNA"

p_demultiplexing_singlets_qc <- ggplot(QC_metrics_long, aes(factor(Demultiplexing_Singlets, c(seq(0,8), "all")), QC_Metric_Value)) +
  geom_violin() +
  stat_summary(fun = "mean", size = 0.5,
                 geom="point", color="firebrick3") +
  theme_bw() +
  facet_grid_sc(QC_Metric ~ ., scales = list(y=scale_y), switch = "y", labeller = as_labeller(c(percent.mt = "Mitochondrial Percent", nCount_RNA = "Number UMIs", nFeature_RNA = "Number Genes")))  +
  ylab(NULL) +
  xlab("Number of Demultiplexing Softwares\nIdentified Droplet as Singlet") +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  force_panelsizes(rows = c(1.5, 1, 1)) +
  geom_text(data = N_demux_dt, aes(y = 85000, label=N, hjust=0,vjust=0,angle=45))

ggsave(p_demultiplexing_singlets_qc, filename = paste0(outdir,"Demultiplexing_Singlets_v_QC.pdf"), width = 10, height = 14, units = "cm")
ggsave(p_demultiplexing_singlets_qc, filename = paste0(outdir,"Demultiplexing_Singlets_v_QC.png"), width = 10, height = 14, units = "cm")



N_demux_doub_dt <- data.table(table(QC_metrics_long$Demultiplexing_Doublets))
colnames(N_demux_doub_dt) <- c("Demultiplexing_Doublets", "N")
N_demux_doub_dt$QC_Metric <- "nCount_RNA"


p_demultiplexing_doublets_qc <- ggplot(QC_metrics_long, aes(factor(Demultiplexing_Doublets, seq(0,8)), QC_Metric_Value)) +
  geom_violin() +
  stat_summary(fun = "mean", size = 0.5,
                 geom="point", color="firebrick3") +
  theme_bw() +
  facet_grid_sc(QC_Metric ~ ., scales = list(y=scale_y), switch = "y", labeller = as_labeller(c(percent.mt = "Mitochondrial Percent", nCount_RNA = "Number UMIs", nFeature_RNA = "Number Genes")))  +
  ylab(NULL) +
  xlab("Number of Demultiplexing Softwares Identified Droplet as Doublet") +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  force_panelsizes(rows = c(1.5, 1, 1)) +
  geom_text(data = N_demux_doub_dt, aes(y = 85000, label=N, hjust=0,vjust=0,angle=45))

ggsave(p_demultiplexing_doublets_qc, filename = paste0(outdir,"Demultiplexing_Doublets_v_QC.pdf"), width = 10, height = 14, units = "cm")
ggsave(p_demultiplexing_doublets_qc, filename = paste0(outdir,"Demultiplexing_Doublets_v_QC.png"), width = 10, height = 14, units = "cm")



N_dd_dt <- data.table(table(QC_metrics_long$DoubletDetecting_Singlets))
colnames(N_dd_dt) <- c("DoubletDetecting_Singlets", "N")
N_dd_dt$QC_Metric <- "nCount_RNA"

p_doublet_detecting_singlets_qc <- ggplot(QC_metrics_long, aes(factor(DoubletDetecting_Singlets, c(seq(0,7), "all")), QC_Metric_Value)) +
  geom_violin() +
  stat_summary(fun = "mean", size = 0.5,
                 geom="point", color="firebrick3") +
  theme_bw() +
  facet_grid_sc(QC_Metric ~ ., scales = list(y=scale_y), switch = "y", labeller = as_labeller(c(percent.mt = "Mitochondrial Percent", nCount_RNA = "Number UMIs", nFeature_RNA = "Number Genes")))  +
  ylab(NULL) +
  xlab("Number of Doublet Detecting Softwares\nIdentified Droplet as Singlet") +
  theme(strip.background = element_blank(),
        strip.placement = "outside")  +
  force_panelsizes(rows = c(1.5, 1, 1)) +
  geom_text(data = N_dd_dt, aes(y = 85000, label=N, hjust=0,vjust=0,angle=45))

ggsave(p_doublet_detecting_singlets_qc, filename = paste0(outdir,"DoubletDetecting_Singlets_v_QC.pdf"), width = 10, height = 14, units = "cm")
ggsave(p_doublet_detecting_singlets_qc, filename = paste0(outdir,"DoubletDetecting_Singlets_v_QC.png"), width = 10, height = 14, units = "cm")


N_dd_doub_dt <- data.table(table(QC_metrics_long$DoubletDetecting_Doublets))
colnames(N_dd_doub_dt) <- c("DoubletDetecting_Doublets", "N")
N_dd_doub_dt$QC_Metric <- "nCount_RNA"


p_doublet_detecting_doublets_qc <- ggplot(QC_metrics_long, aes(factor(DoubletDetecting_Doublets, seq(0,7)), QC_Metric_Value)) +
  geom_violin() +
  stat_summary(fun = "mean", size = 0.5,
                 geom="point", color="firebrick3") +
  theme_bw() +
  facet_grid_sc(QC_Metric ~ ., scales = list(y=scale_y), switch = "y", labeller = as_labeller(c(percent.mt = "Mitochondrial Percent", nCount_RNA = "Number UMIs", nFeature_RNA = "Number Genes")))  +
  ylab(NULL) +
  xlab("Number of Doublet Detecting Softwares\nIdentified Droplet as Doublet") +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  force_panelsizes(rows = c(1.5, 1, 1)) +
  geom_text(data = N_dd_doub_dt, aes(y = 85000, label=N, hjust=0,vjust=0,angle=45))

ggsave(p_doublet_detecting_doublets_qc, filename = paste0(outdir,"DoubletDetecting_Doublets_v_QC.pdf"), width = 10, height = 14, units = "cm")
ggsave(p_doublet_detecting_doublets_qc, filename = paste0(outdir,"DoubletDetecting_Doublets_v_QC.png"), width = 10, height = 14, units = "cm")



##### Singlets for demultiplexing and doublet detecting faceted #####
QC_metrics_long_singlets <- data.table(melt(QC_metrics_long, id.vars = c("Pool", "QC_Metric", "QC_Metric_Value"),
                measure.vars = c("Demultiplexing_Singlets", "DoubletDetecting_Singlets"),
								variable.name = "Singlets", value.name = "Software_N"))

N_demux_dt_updated <- cbind(N_demux_dt, data.table(Singlets = "Demultiplexing_Singlets"))
colnames(N_demux_dt_updated) <- gsub("Demultiplexing_Singlets", "Software_N", colnames(N_demux_dt_updated))
N_dd_dt_updated <- cbind(N_dd_dt, data.table(Singlets = "DoubletDetecting_Singlets"))
colnames(N_dd_dt_updated) <- gsub("DoubletDetecting_Singlets", "Software_N", colnames(N_dd_dt_updated))

N_sing_demux_dd_dt <- rbind(N_demux_dt_updated, N_dd_dt_updated)

p_demultiplexing_doubletdetecting_singlets_qc <- ggplot(QC_metrics_long_singlets, aes(factor(Software_N, c(seq(0,8), "all")), QC_Metric_Value)) +
  geom_violin() +
  stat_summary(fun = "mean", size = 0.5,
                 geom="point", color="firebrick3") +
  theme_bw() +
  facet_grid_sc(QC_Metric ~ Singlets, scales = "free", switch = "y", labeller = as_labeller(c(percent.mt = "Mitochondrial\nPercent", nCount_RNA = "Number\nUMIs", nFeature_RNA = "Number\nGenes", Demultiplexing_Singlets = "Demultiplexing\nSoftwares", DoubletDetecting_Singlets = "Doublet Detecting\nSoftwares")))  +
  ylab(NULL) +
  xlab("Number Softwares\nIdentified Droplet as Singlet") +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  force_panelsizes(rows = c(1.5, 1, 1)) +
  geom_text(data = N_sing_demux_dd_dt, aes(y = 40000, label=N, hjust=0,vjust=0,angle=45))

ggsave(p_demultiplexing_doubletdetecting_singlets_qc, filename = paste0(outdir,"Demultiplexing_DoubletDetecting_Singlets_v_QC.pdf"), width = 10, height = 8, units = "cm")
ggsave(p_demultiplexing_doubletdetecting_singlets_qc, filename = paste0(outdir,"Demultiplexing_DoubletDetecting_Singlets_v_QC.png"), width = 10, height = 8, units = "cm")




##### Plot % agreement per software with facet for number #####
singlet_doublet_long_list <- lapply(results_list_wide, function(x){
    columns <- c("Barcode", grep("_DropletType", colnames(x), value = TRUE))
    temp <- data.table(pivot_longer(x[,..columns], cols = grep("_DropletType", colnames(x), value = TRUE), values_to = "DropletType", names_to = "Software"))
    temp$Software <- gsub("_DropletType", "", temp$Software) %>%
        gsub("demux", "Demux", .) %>%
        gsub("dropulation", "Dropulation", .) %>%
        gsub("Demuxalot_refined", "Demuxalot (refined)", .) %>%
        gsub("freemuxlet", "Freemuxlet", .) %>%
        gsub("sc", "Sc", .) %>%
        gsub("vireo", "Vireo", .) %>%
        gsub("solo", "Solo", .) %>%
        gsub("souporcell", "Souporcell", .)
    return(temp)
})


proportion_list <- lapply(names(singlet_doublet_long_list), function(x){
  temp <- data.table(prop.table(table(singlet_doublet_long_list[[x]]$Software, singlet_doublet_long_list[[x]]$DropletType), margin = 1))
  temp$Pool <- x
  colnames(temp) <- gsub("V1", "Software", colnames(temp)) %>%
                      gsub("V2", "DropletType", .)
  return(temp)
})


proportion <- do.call(rbind, proportion_list)
proportion$DropletType <- factor(proportion$DropletType, levels = rev(c("singlet", "doublet", "unassigned")))
proportion$N <- as.numeric(proportion$N)*100
proportion$Software <- factor(proportion$Software, levels = software_names)

p_proportion <- ggbarplot(proportion, 
                                  x = "Software", 
                                  y = "N", 
                                  add = "mean_se",
                                  fill = "DropletType")  +
                          ylab("Percent") +
                          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=9)) +
                          scale_fill_manual(values = c("grey85", "#5C3E9E", "#47A448")) +
                          xlab(NULL) +
                          scale_y_continuous(expand = c(0,0))

ggsave(p_proportion, filename = paste0(outdir, "singlets_doublets_proportions.png"), width = 2.1, height = 3.75)
ggsave(p_proportion, filename = paste0(outdir, "singlets_doublets_proportions.pdf"), width = 2.1, height = 3.75)


##### Demultiplexing #####
singlets_long_list <- lapply(results_list_wide, function(x){
    columns <- c("Barcode", grep("_DropletType", colnames(x), value = TRUE), "Singlets", "Demultiplexing_Singlets", "DoubletDetecting_Singlets")
    temp <- data.table(pivot_longer(x[,..columns], cols = grep("_DropletType", colnames(x), value = TRUE), values_to = "DropletType", names_to = "Software"))
    temp <- temp[DropletType == "singlet"]
    temp$Software <- gsub("_DropletType", "", temp$Software) %>%
        gsub("demux", "Demux", .) %>%
        gsub("dropulation", "Dropulation", .) %>%
        gsub("Demuxalot_refined", "Demuxalot (refined)", .) %>%
        gsub("freemuxlet", "Freemuxlet", .) %>%
        gsub("sc", "Sc", .) %>%
        gsub("vireo", "Vireo", .) %>%
        gsub("solo", "Solo", .) %>%
        gsub("souporcell", "Souporcell", .)
    return(temp)
})

singlets_proportion_list <- lapply(names(singlets_long_list), function(x){
  temp <- data.table(prop.table(table(singlets_long_list[[x]]$Software, singlets_long_list[[x]]$Singlets), margin = 1))
  temp$Pool <- x
  colnames(temp) <- gsub("V1", "Software", colnames(temp)) %>%
                      gsub("V2", "Singlets", .)
  return(temp)
})

singlets_proportion <- do.call(rbind, singlets_proportion_list)
singlets_proportion$Singlets <- factor((as.numeric(as.character(singlets_proportion$Singlets)) - 1), levels = 0:14)
singlets_proportion$N <- as.numeric(singlets_proportion$N)*100
singlets_proportion$Software <- factor(singlets_proportion$Software, levels = software_names)

singlets_proportions <- ggbarplot(singlets_proportion, 
                                  x = "Software", 
                                  y = "N", 
                                  add = "mean_se",
                                  fill = "Singlets")  +
                          ylab("Percent") +
                          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                          xlab(NULL) +
                          scale_y_continuous(expand = c(0,0)) +
                          scale_fill_discrete_sequential(palette = "Greens") +
                          labs(fill="Number other Methods\nClassified as Singlet")

ggsave(singlets_proportions, filename = paste0(outdir, "singlets_agreement_proportions.png"), width = 3.5, height = 6)
ggsave(singlets_proportions, filename = paste0(outdir, "singlets_agreement_proportions.pdf"), width = 3.5, height = 6)





## Demultiplexing Singlets
singlets_demux_long_list <- lapply(singlets_long_list, function(x){
   x[Software %in% demultiplexing_names]
})

singlets_demux_proportion_list <- lapply(names(singlets_demux_long_list), function(x){
  temp <- data.table(prop.table(table(singlets_demux_long_list[[x]]$Software, singlets_demux_long_list[[x]]$Demultiplexing_Singlets), margin = 1))
  temp$Pool <- x
  colnames(temp) <- gsub("V1", "Software", colnames(temp)) %>%
                      gsub("V2", "Demultiplexing_Singlets", .)
  return(temp)
})

singlets_demux_proportion <- do.call(rbind, singlets_demux_proportion_list)
singlets_demux_proportion$Demultiplexing_Singlets <- factor(as.numeric(as.character(singlets_demux_proportion$Demultiplexing_Singlets)) - 1, levels = 0:8)
singlets_demux_proportion$N <- as.numeric(singlets_demux_proportion$N)*100
singlets_demux_proportion$Software <- factor(singlets_demux_proportion$Software, levels = demultiplexing_names)

p_singlets_demux_proportions <- ggbarplot(singlets_demux_proportion, 
                                  x = "Software", 
                                  y = "N", 
                                  add = "mean_se",
                                  fill = "Demultiplexing_Singlets")  +
                          ylab("Percent") +
                          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                          xlab(NULL) +
                          scale_y_continuous(expand = c(0,0)) +
                          scale_fill_discrete_sequential(palette = "Greens") +
                          labs(fill="Number other\nDemultiplexing Methods\nClassified as Singlet")

ggsave(p_singlets_demux_proportions, filename = paste0(outdir, "singlets_demux_agreement_proportions.png"), width = 2.5, height = 6)
ggsave(p_singlets_demux_proportions, filename = paste0(outdir, "singlets_demux_agreement_proportions.pdf"), width = 2.5, height = 6)

## Doublet Detecting Singlets
singlets_dd_long_list <- lapply(singlets_long_list, function(x){
   x[Software %in% doubletdetecting_names]
})

singlets_dd_proportion_list <- lapply(names(singlets_dd_long_list), function(x){
  temp <- data.table(prop.table(table(singlets_dd_long_list[[x]]$Software, singlets_dd_long_list[[x]]$DoubletDetecting_Singlets), margin = 1))
  temp$Pool <- x
  colnames(temp) <- gsub("V1", "Software", colnames(temp)) %>%
                      gsub("V2", "DoubletDetecting_Singlets", .)
  return(temp)
})

singlets_dd_proportion <- do.call(rbind, singlets_dd_proportion_list)
singlets_dd_proportion$DoubletDetecting_Singlets <- factor(as.numeric(singlets_dd_proportion$DoubletDetecting_Singlets) - 1, levels = 0:7)
singlets_dd_proportion$N <- as.numeric(singlets_dd_proportion$N)*100
singlets_dd_proportion$Software <- factor(singlets_dd_proportion$Software, levels = doubletdetecting_names)

p_singlets_dd_proportions <- ggbarplot(singlets_dd_proportion, 
                                  x = "Software", 
                                  y = "N", 
                                  add = "mean_se",
                                  fill = "DoubletDetecting_Singlets")  +
                          ylab("Percent") +
                          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                          xlab(NULL) +
                          scale_y_continuous(expand = c(0,0)) +
                          scale_fill_discrete_sequential(palette = "Greens") +
                          labs(fill="Number other\nDoublet Detecting Methods\nClassified as Singlet")

ggsave(p_singlets_dd_proportions, filename = paste0(outdir, "singlets_dd_agreement_proportions.png"), width = 2.5, height = 6)
ggsave(p_singlets_dd_proportions, filename = paste0(outdir, "singlets_dd_agreement_proportions.pdf"), width = 2.5, height = 6)



### Doublets ###
doublets_long_list <- lapply(results_list_wide, function(x){
  columns <- c("Barcode", grep("_DropletType", colnames(x), value = TRUE), "Doublets", "Demultiplexing_Doublets", "DoubletDetecting_Doublets")
    temp <- data.table(pivot_longer(x[,..columns], cols = grep("_DropletType", colnames(x), value = TRUE), values_to = "DropletType", names_to = "Software"))
    temp <- temp[DropletType == "doublet"]
    temp$Software <- gsub("_DropletType", "", temp$Software)
    temp$Software <- gsub("_DropletType", "", temp$Software) %>%
        gsub("demux", "Demux", .) %>%
        gsub("dropulation", "Dropulation", .) %>%
        gsub("Demuxalot_refined", "Demuxalot (refined)", .) %>%
        gsub("freemuxlet", "Freemuxlet", .) %>%
        gsub("sc", "Sc", .) %>%
        gsub("vireo", "Vireo", .) %>%
        gsub("solo", "Solo", .) %>%
        gsub("souporcell", "Souporcell", .)
    return(temp)
})


doublets_proportion_list <- lapply(names(doublets_long_list), function(x){
  temp <- data.table(prop.table(table(doublets_long_list[[x]]$Software, doublets_long_list[[x]]$Doublets), margin = 1))
  temp$Pool <- x
  colnames(temp) <- gsub("V1", "Software", colnames(temp)) %>%
                      gsub("V2", "Doublets", .)
  return(temp)
})

doublets_proportion <- do.call(rbind, doublets_proportion_list)
doublets_proportion$Doublets <- factor(as.numeric(as.character(doublets_proportion$Doublets)) -1, levels = 0:14)
doublets_proportion$N <- as.numeric(doublets_proportion$N)*100
doublets_proportion$Software <- factor(doublets_proportion$Software, levels = software_names)

p_doublets_proportion <- ggbarplot(doublets_proportion, 
                                  x = "Software", 
                                  y = "N", 
                                  add = "mean_se",
                                  fill = "Doublets")  +
                          ylab("Percent") +
                          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                          xlab(NULL) +
                          scale_y_continuous(expand = c(0,0)) +
                          scale_fill_discrete_sequential(palette = "Purples")

ggsave(p_doublets_proportion, filename = paste0(outdir, "doublets_agreement_proportions.png"), width = 3.5, height = 6)
ggsave(p_doublets_proportion, filename = paste0(outdir, "doublets_agreement_proportions.pdf"), width = 3.5, height = 6)


## Demultiplexing softwares ##
doublets_demux_long_list <- lapply(doublets_long_list, function(x){
    x[Software %in% demultiplexing_names]
})


doublets_demux_proportion_list <- lapply(names(doublets_demux_long_list), function(x){
  temp <- data.table(prop.table(table(doublets_demux_long_list[[x]]$Software, doublets_demux_long_list[[x]]$Demultiplexing_Doublets), margin = 1))
  temp$Pool <- x
  colnames(temp) <- gsub("V1", "Software", colnames(temp)) %>%
                      gsub("V2", "Demultiplexing_Doublets", .)
  return(temp)
})

doublets_demux_proportion <- do.call(rbind, doublets_demux_proportion_list)
doublets_demux_proportion$Demultiplexing_Doublets <- factor(as.numeric(as.character(doublets_demux_proportion$Demultiplexing_Doublets)) - 1, levels = 0:7)
doublets_demux_proportion$N <- as.numeric(doublets_demux_proportion$N)*100
doublets_demux_proportion$Software <- factor(doublets_demux_proportion$Software, levels = demultiplexing_names)

p_doublets_demux_proportion <- ggbarplot(doublets_demux_proportion, 
                                  x = "Software", 
                                  y = "N", 
                                  add = "mean_se",
                                  fill = "Demultiplexing_Doublets")  +
                          ylab("Percent") +
                          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                          xlab(NULL) +
                          scale_y_continuous(expand = c(0,0)) +
                          scale_fill_discrete_sequential(palette = "Purples")

ggsave(p_doublets_demux_proportion, filename = paste0(outdir, "doublets_demultiplexing_agreement_proportions.png"), width = 3.5, height = 6)
ggsave(p_doublets_demux_proportion, filename = paste0(outdir, "doublets_demultiplexing_agreement_proportions.pdf"), width = 3.5, height = 6)


## Doublet Detecting softwares ##
doublets_dd_long_list <- lapply(doublets_long_list, function(x){
    x[Software %in% doubletdetecting_names]
})


doublets_dd_proportion_list <- lapply(names(doublets_dd_long_list), function(x){
  temp <- data.table(prop.table(table(doublets_dd_long_list[[x]]$Software, doublets_dd_long_list[[x]]$DoubletDetecting_Doublets), margin = 1))
  temp$Pool <- x
  colnames(temp) <- gsub("V1", "Software", colnames(temp)) %>%
                      gsub("V2", "DoubletDetecting_Doublets", .)
  return(temp)
})

doublets_dd_proportion <- do.call(rbind, doublets_dd_proportion_list)
doublets_dd_proportion$DoubletDetecting_Doublets <- factor(as.numeric(as.character(doublets_dd_proportion$DoubletDetecting_Doublets)) - 1, levels = 0:7)
doublets_dd_proportion$N <- as.numeric(doublets_dd_proportion$N)*100
doublets_dd_proportion$Software <- factor(doublets_dd_proportion$Software, levels = doubletdetecting_names)

p_doublets_dd_proportion <- ggbarplot(doublets_dd_proportion, 
                                  x = "Software", 
                                  y = "N", 
                                  add = "mean_se",
                                  fill = "DoubletDetecting_Doublets")  +
                          ylab("Percent") +
                          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                          xlab(NULL) +
                          scale_y_continuous(expand = c(0,0)) +
                          scale_fill_discrete_sequential(palette = "Purples")

ggsave(p_doublets_dd_proportion, filename = paste0(outdir, "doublets_dd_agreement_proportions.png"), width = 3.5, height = 6)
ggsave(p_doublets_dd_proportion, filename = paste0(outdir, "doublets_dd_agreement_proportions.pdf"), width = 3.5, height = 6)



##### Facet singlets and doublets to one figure #####
p_singlets_doublets_proportions <- cowplot::plot_grid(p_doublets_proportion + theme(text = element_text(size=9), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank()), singlets_proportions + theme(text = element_text(size=9), legend.position = "none"),  align = "v", ncol = 1, rel_heights = c(1,1.675))

ggsave(p_singlets_doublets_proportions, filename = paste0(outdir, "stacked_proportions_singlets_doublets.png"), width = 2.1, height = 4.2)
ggsave(p_singlets_doublets_proportions, filename = paste0(outdir, "stacked_proportions_singlets_doublets.pdf"), width = 2.1, height = 4.2)




### Consider the majority call correct and see how each software compares ###
individual_assignment_list <- future_apply(results_wide[,.SD, .SDcols = grep("_Assignment", colnames(results_wide), value = TRUE)], 1, function(y) table(y)[!(rownames(table(y)) %in% c("unassigned","doublet"))])

individual_assignment <- data.table(MajorityAssignment_ID = unlist(lapply(individual_assignment_list, function(y) ifelse(is.null(names(which.max.simple(y[names(y) != "doublet"]))), NA, names(which.max.simple(y[names(y) != "doublet"]))))),
                  MajorityAssignment_N =  unlist(lapply(individual_assignment_list, function(y) ifelse(is.null(names(which.max.simple(y[names(y) != "doublet"]))), NA, max(y[names(y) != "doublet"]))))) ## Need to remove doublet count from table

results_wide <- cbind(results_wide, individual_assignment)

results_wide$MajorityDropletType <- ifelse(results_wide$Singlets >= 8 & results_wide$MajorityAssignment_N >= 4, "singlet", 
                                      ifelse(results_wide$Doublets >=8, "doublet", "unassigned"))




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


for (pool in unique(results_wide$Pool)){
	# Metrics[[pool]] <- data.table(matrix(nrow = length(softwares), ncol = 13))
	Metrics[[pool]] <- data.table("software" = softwares, "true_singlet" = as.numeric(NA), "false_singlet" = as.numeric(NA), "true_doublet" = as.numeric(NA), "false_doublet" = as.numeric(NA), "accuracy" = as.numeric(NA), "balanced_accuracy" = as.numeric(NA), "sensitivity" = as.numeric(NA), "specificity" = as.numeric(NA), "npv" = as.numeric(NA), "ppv" = as.numeric(NA), "precision" = as.numeric(NA), "mcc" = as.numeric(NA))
  Metrics[[pool]]$pool <- pool


	for (soft in softwares){
		## True singlet (TP)
    column <- paste0(soft, "_DropletType")
		if (soft %in% demultiplexing_softwares){
      demux_column <- paste0(soft, "_Assignment")
			Metrics[[pool]][software == soft,]$true_singlet <- as.numeric(length(which(results_wide[Pool == pool, ..column] == "singlet" & results_wide[Pool == pool, ..column] == results_wide[Pool == pool]$MajorityDropletType & results_wide[Pool == pool, ..demux_column] == results_wide[Pool == pool]$MajorityAssignment_ID)))
		} else {
			Metrics[[pool]][software == soft,]$true_singlet <- as.numeric(length(which(results_wide[Pool == pool, ..column] == "singlet" & results_wide[Pool == pool, ..column] == results_wide[Pool == pool]$MajorityDropletType)))
		}
		# Metrics[[pool]][software == soft,]$true_singlet[Metrics[[pool]][software == soft,]$true_singlet == 0] <- NA

		## False singlet (FP)
		Metrics[[pool]][software == soft,]$false_singlet <- as.numeric(length(which(results_wide[Pool == pool, ..column] == "singlet")) - Metrics[[pool]][software == soft,]$true_singlet)
		# Metrics[[pool]][software == soft,]$false_singlet[Metrics[[pool]][software == soft,]$false_singlet == 0] <- NA

		### Update NA for pools where could not run (ie DoubletDecon)
		if (Metrics[[pool]][software == soft,]$true_singlet == 0 & Metrics[[pool]][software == soft,]$false_singlet == 0){
			Metrics[[pool]][software == soft,]$true_doublet <- NA
			Metrics[[pool]][software == soft,]$false_doublet <- NA
			Metrics[[pool]][software == soft,]$accuracy <- NA
			Metrics[[pool]][software == soft,]$sensitivity <- NA
			Metrics[[pool]][software == soft,]$specificity <- NA
			Metrics[[pool]][software == soft,]$npv <- NA
			Metrics[[pool]][software == soft,]$ppv <- NA
			Metrics[[pool]][software == soft,]$balanced_accuracy <- NA
			Metrics[[pool]][software == soft,]$precision <- NA
			Metrics[[pool]][software == soft,]$mcc <- NA
		} else{

			## True doublet (TN)
			Metrics[[pool]][software == soft,]$true_doublet <- as.numeric(length(which(results_wide[Pool == pool, ..column] == "doublet" & results_wide[Pool == pool, ..column] == results_wide[Pool == pool]$MajorityDropletType)))


			## False doublet (FN)
			Metrics[[pool]][software == soft,]$false_doublet <- as.numeric(length(which(results_wide[Pool == pool, ..column] == "doublet")) - Metrics[[pool]][software == soft,]$true_doublet)


			## Accuracy
			Metrics[[pool]][software == soft,]$accuracy <- as.numeric((Metrics[[pool]][software == soft,]$true_singlet + Metrics[[pool]][software == soft,]$true_doublet)/nrow(results_wide[Pool == pool]))


			## Sensitivity (TPR)
			Metrics[[pool]][software == soft,]$sensitivity <- as.numeric(Metrics[[pool]][software == soft,]$true_singlet/length(which(results_wide[Pool == pool]$MajorityDropletType == "singlet")))


			## Specificity (TNR)
			Metrics[[pool]][software == soft,]$specificity <- as.numeric(Metrics[[pool]][software == soft,]$true_doublet/length(which(results_wide[Pool == pool]$MajorityDropletType == "doublet")))


			## Negative Predictive Value
			Metrics[[pool]][software == soft,]$npv <- as.numeric(Metrics[[pool]][software == soft,]$true_doublet/length(which(results_wide[Pool == pool, ..column] == "doublet")))


			## Positive Predictive Value
			Metrics[[pool]][software == soft,]$ppv <- as.numeric(Metrics[[pool]][software == soft,]$true_singlet/length(which(results_wide[Pool == pool, ..column] == "singlet")))


			## Balanced Accuracy
			Metrics[[pool]][software == soft,]$balanced_accuracy <- as.numeric((Metrics[[pool]][software == soft,]$sensitivity + Metrics[[pool]][software == soft,]$specificity)/2)


			## Precision
			Metrics[[pool]][software == soft,]$precision <- as.numeric(Metrics[[pool]][software == soft,]$true_singlet/(Metrics[[pool]][software == soft,]$true_singlet + Metrics[[pool]][software == soft,]$false_singlet))

			
			## MCC
			Metrics[[pool]][software == soft,]$mcc <- as.numeric(((Metrics[[pool]][software == soft,]$true_singlet * Metrics[[pool]][software == soft,]$true_doublet) - (Metrics[[pool]][software == soft,]$false_singlet * Metrics[[pool]][software == soft,]$false_doublet))/sqrt((Metrics[[pool]][software == soft,]$true_singlet + Metrics[[pool]][software == soft,]$false_singlet) * (Metrics[[pool]][software == soft,]$true_singlet + Metrics[[pool]][software == soft,]$false_doublet) * (Metrics[[pool]][software == soft,]$true_doublet + Metrics[[pool]][software == soft,]$false_singlet) * (Metrics[[pool]][software == soft,]$true_doublet + Metrics[[pool]][software == soft,]$false_doublet)))
		}
	}
}

Metrics_dt <- do.call(rbind, Metrics)


Metrics_dt_long <- melt(Metrics_dt,
                measure.vars =c("accuracy", "balanced_accuracy", "sensitivity", "specificity", "npv", "ppv", "precision", "mcc"),
								variable.name = "Metric", value.name = "Value")

Metrics_dt_long$software <-  gsub("demux", "Demux", Metrics_dt_long$software) %>%
        gsub("dropulation", "Dropulation", .) %>%
        gsub("Demuxalot_refined", "Demuxalot (refined)", .) %>%
        gsub("freemuxlet", "Freemuxlet", .) %>%
        gsub("sc", "Sc", .) %>%
        gsub("vireo", "Vireo", .) %>%
        gsub("solo", "Solo", .) %>%
        gsub("souporcell", "Souporcell", .)

Metrics_dt_long$software <- factor(Metrics_dt_long$software, levels = software_names)

Metrics_dt_long$Metric <- gsub("mcc", "MCC", Metrics_dt_long$Metric) %>%
                            gsub("sensitivity", "Sensitivity", .) %>%
                            gsub("specificity", "Sepcificity", .) %>%
                            gsub("precision", "Precision", .) %>%
                            gsub("npv", "NPV", .) %>%
                            gsub("ppv", "PPV", .)


Pmajority_metrics <- ggplot(Metrics_dt_long[Metric %in% c("MCC", "Sensitivity", "Specificity", "Precision", "NPV", "PPV")], aes(software, Value, color = software)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5, width = 0.1) +
  theme_classic() +
  facet_grid(Metric ~ ., scales = "free_y") +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab(NULL) +
  theme(legend.position = "none")+
  ylab("Metric Value")


ggsave(Pmajority_metrics, filename = paste0(outdir, "facet_metrics_majority.png"), width = 3, height = 5.5)
ggsave(Pmajority_metrics, filename = paste0(outdir, "facet_metrics_majority.pdf"), width = 3, height = 5.5)



Pmajority_MCC <- ggplot(Metrics_dt_long[Metric %in% c("MCC")], aes(software, Value, color = software)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.75, width = 0.1) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("MCC") +
  xlab(NULL) +
  theme(legend.position = "none", text = element_text(size=11))

ggsave(Pmajority_MCC, filename = paste0(outdir, "facet_MCC_majority.png"), height = 2.5, width = 4)
ggsave(Pmajority_MCC, filename = paste0(outdir, "facet_MCC_majority.pdf"), height = 2.5, width = 4)





##### Make a dataframe that has the number of singlets, doublets, unassigned in each pool for each software #####
singlets_temp <- data.table(table(results_wide[Singlets == 15 & MajorityAssignment_N == 8 & !(Pool %in% c("OneK1K_scRNA_Sample40", "OneK1K_scRNA_Sample48", "OneK1K_scRNA_Sample66"))][,c("Pool", "MajorityAssignment_ID")]))[N > 0]
colnames(singlets_temp) <- c("Pool","Individual", "Count")


doublets_temp <- data.table(table(results_wide[Doublets == 15 & !(Pool %in% c("OneK1K_scRNA_Sample40", "OneK1K_scRNA_Sample48", "OneK1K_scRNA_Sample66"))][,c("Pool", "MajorityDropletType")]))
doublets_temp <- rbind(doublets_temp, data.table(Pool = unique(results_wide$Pool[!(results_wide$Pool %in% doublets_temp$Pool)]), MajorityDropletType = "doublet", N = 0))
colnames(doublets_temp) <- c("Pool","Individual", "Count")

results_singlets_4df <- rbind(singlets_temp, doublets_temp)



pSingletsBox <- ggplot(singlets_temp, aes(x=factor(Pool), y = Count)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(filename = paste0(out, "SingletsBoxPlot.png"), plot = pSingletsBox, width = 16, height = 9)



### Add N individuals in each pool ###
## Read in individuals files ##
individuals <- lapply(unique(singlets_temp$Pool), function(x){
  temp <- fread(paste0(dir, x, "/popscle/Individuals.txt"), col.names = c("Individual"), header = FALSE)
  temp2 <- data.table(Pool = x, N_genotyped = nrow(temp))
  return(temp2)
})
individual_n <- do.call(rbind, individuals)

singlets_temp <- singlets_temp[individual_n, on = c("Pool")]

N_individuals <- data.table(table(singlets_temp$Pool))
colnames(N_individuals) <- c("Pool", "N_identified")

singlets_temp <- singlets_temp[N_individuals, on = "Pool"]

singlets_temp <- singlets_temp[order(N_genotyped, N_identified)]
singlets_temp$Pool <-  factor(singlets_temp$Pool, levels = unique(singlets_temp$Pool))



doublets_temp <-doublets_temp[data.table(Pool = levels(singlets_temp$Pool)), on = "Pool"]
doublets_temp$Pool <-  factor(doublets_temp$Pool, levels = unique(singlets_temp$Pool))

colors <- list("Number Individuals in Pool" = brewer.pal(length(unique(singlets_temp$N_genotyped)),"Blues"))
colors[["Number Individuals in Pool"]] <- factor(colors[["Number Individuals in Pool"]][1:length(unique(singlets_temp$N_genotyped))], levels = colors[["Number Individuals in Pool"]][1:length(unique(singlets_temp$N_genotyped))])
names(colors[["Number Individuals in Pool"]]) <- unique(singlets_temp$N_genotyped)

pSingletsBox_doubletBar <- ggplot() +
  geom_bar(data = doublets_temp, aes(x=Pool, y = Count, fill = ""), stat="identity") +
  scale_fill_manual("Doublets", values = "grey37") +
  new_scale("fill") +
  geom_boxplot(data = singlets_temp, aes(x=Pool, y = Count, fill = as.factor(N_genotyped), color = ""), outlier.size = 0.5) +
  scale_color_manual("Singlets Per\nIndividual", values=1) +
  scale_fill_manual(values = c("#EFF3FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"), name = "Number\nIndividuals\nin Pool") +
  theme_classic() +
  theme(text = element_text(size=18),
      legend.position="right", 
      legend.box = "verticle",
      plot.title = element_text(hjust = 0.5),
      legend.background = element_blank(),
      legend.box.background = element_rect(colour = "black"),
      axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # ggtitle("Number of Singlets and Doublets\nCommon Across All Demultiplexing Softwares") +
  ylab("Number Droplets")  +
  guides(color = guide_legend(order = 1, nrow = 1),
         fill = guide_legend(order = 2, nrow = 3)) +
  scale_y_continuous(expand = c(0, 0))
  
ggsave(filename = paste0(outdir, "SingletsBoxPlot_DoubletBar_genotyped_test.png"), plot = pSingletsBox_doubletBar, width = 8.25, height = 2.3)
ggsave(filename = paste0(outdir, "SingletsBoxPlot_DoubletBar_genotyped_test.pdf"), plot = pSingletsBox_doubletBar, width = 8.25, height = 2.3)






##### Make a dataframe that has the number of unique inidividuals identified per software per pool #####
Indiv_identified <- data.table(Pool = unique(results_wide[!(Pool %in% c("OneK1K_scRNA_Sample40", "OneK1K_scRNA_Sample48", "OneK1K_scRNA_Sample66"))]$Pool),"demuxalot" = as.numeric(NA), "demuxalot_refined" = as.numeric(NA), "demuxlet" = as.numeric(NA), "dropulation" = as.numeric(NA), "freemuxlet" = as.numeric(NA),"scSplit" = as.numeric(NA),"souporcell" = as.numeric(NA),"vireo" = as.numeric(NA))


for (pool in Indiv_identified$Pool){
  Indiv_identified[which(Indiv_identified$Pool == pool),"demuxalot"] <- length(unique(results_wide[Pool == pool]$demuxalot_Assignment)[!(unique(results_wide[Pool == pool]$demuxalot_Assignment) %in% c("doublet","unsure","unassigned"))])
  Indiv_identified[which(Indiv_identified$Pool == pool),"demuxalot_refined"] <- length(unique(results_wide[Pool == pool]$demuxalot_refined_Assignment)[!(unique(results_wide[Pool == pool]$demuxalot_refined_Assignment) %in% c("doublet","unsure","unassigned"))])
  Indiv_identified[which(Indiv_identified$Pool == pool),"demuxlet"] <- length(unique(results_wide[Pool == pool]$demuxlet_Assignment)[!(unique(results_wide[Pool == pool]$demuxlet_Assignment) %in% c("doublet","unsure","unassigned"))])
  Indiv_identified[which(Indiv_identified$Pool == pool),"dropulation"] <- length(unique(results_wide[Pool == pool]$dropulation_Assignment)[!(unique(results_wide[Pool == pool]$dropulation_Assignment) %in% c("doublet","unsure","unassigned"))])
  Indiv_identified[which(Indiv_identified$Pool == pool),"freemuxlet"] <- length(unique(results_wide[Pool == pool]$freemuxlet_Assignment)[!(unique(results_wide[Pool == pool]$freemuxlet_Assignment) %in% c("doublet","unsure","unassigned"))])
  Indiv_identified[which(Indiv_identified$Pool == pool),"scSplit"] <- length(unique(results_wide[Pool == pool]$scSplit_Assignment)[!(unique(results_wide[Pool == pool]$scSplit_Assignment) %in% c("doublet","unsure","unassigned"))])
  Indiv_identified[which(Indiv_identified$Pool == pool),"souporcell"] <- length(unique(results_wide[Pool == pool]$souporcell_Assignment)[!(unique(results_wide[Pool == pool]$souporcell_Assignment) %in% c("doublet","unsure","unassigned"))])
  Indiv_identified[which(Indiv_identified$Pool == pool),"vireo"] <- length(unique(results_wide[Pool == pool]$vireo_Assignment)[!(unique(results_wide[Pool == pool]$vireo_Assignment) %in% c("doublet","unsure","unassigned"))][!str_detect(unique(results_list_wide[[pool]]$vireo_Assignment)[!(unique(results_list_wide[[pool]]$vireo_Assignment) %in% c("doublet","unsure","unassigned"))],pattern="donor")])
}



Indiv_identified <- left_join(Indiv_identified, meta, by = "Pool")
Indiv_identified$Genotyped <- NA
for (row in 1:nrow(Indiv_identified)){
  Indiv_identified$Genotyped[row] <- (str_count(Indiv_identified$Individuals[row], pattern = ",") + 1)
}
Indiv_identified$Individuals <- NULL
Indiv_identified$Barcodes <- NULL


Indiv_missed <- Indiv_identified
Indiv_missed$Demuxalot_percentage <- Indiv_identified$Genotyped- Indiv_identified$demuxalot
Indiv_missed$`Demuxalot (refined)_percentage` <- Indiv_identified$Genotyped- Indiv_identified$`demuxalot_refined`
Indiv_missed$Demuxlet_percentage <- Indiv_identified$Genotyped- Indiv_identified$demuxlet
Indiv_missed$Dropulation_percentage <- Indiv_identified$Genotyped- Indiv_identified$dropulation
Indiv_missed$Freemuxlet_percentage <- Indiv_identified$Genotyped- Indiv_identified$freemuxlet
Indiv_missed$ScSplit_percentage <- Indiv_identified$Genotyped- Indiv_identified$scSplit
Indiv_missed$Souporcell_percentage <- Indiv_identified$Genotyped- Indiv_identified$souporcell
Indiv_missed$Vireo_percentage <- Indiv_identified$Genotyped- Indiv_identified$vireo
Indiv_missed$TotalMissing <- Indiv_missed$Demuxalot_percentage + Indiv_missed$`Demuxalot (refined)_percentage` + Indiv_missed$Demuxlet_percentage + Indiv_missed$Dropulation_percentage + Indiv_missed$Freemuxlet_percentage + Indiv_missed$ScSplit_percentage + Indiv_missed$Souporcell_percentage + Indiv_missed$Vireo_percentage
Indiv_missed$demuxalot <- NULL
Indiv_missed$demuxalot_refined <- NULL
Indiv_missed$demuxlet <- NULL
Indiv_missed$dropulation <- NULL
Indiv_missed$freemuxlet <- NULL
Indiv_missed$scSplit <- NULL
Indiv_missed$souporcell <- NULL
Indiv_missed$vireo <- NULL
colnames(Indiv_missed) <- gsub("_percentage","",colnames(Indiv_missed))
Indiv_missed <- Indiv_missed[order( Genotyped, TotalMissing ),]


ht_opt(legend_title_position = "topleft", 
  legend_labels_gp = gpar(fontsize = 18),
  legend_title_gp = gpar(fontsize = 18, fontface = "bold"),
  heatmap_row_names_gp = gpar(fontsize = 14),
  heatmap_column_title_gp = gpar(fontsize = 18))
colors <- list("Number Individuals in Pool" = brewer.pal(length(unique(Indiv_missed$Genotyped)),"Blues"))
colors[["Number Individuals in Pool"]] <- factor(colors[["Number Individuals in Pool"]][1:length(unique(Indiv_missed$Genotyped))], levels = colors[["Number Individuals in Pool"]][1:length(unique(Indiv_missed$Genotyped))])
names(colors[["Number Individuals in Pool"]]) <- unique(Indiv_missed$Genotyped)
column_ha = HeatmapAnnotation("Number Individuals in Pool" = Indiv_missed$Genotyped, col = colors,
  annotation_name_side = "left",
  annotation_legend_param = list(
        "Number Individuals in Pool" = list(
                title = "Number Individuals in Pool",
                at = names(colors[["Number Individuals in Pool"]]),
                nrow = 1
            )))

matrix <- as.matrix(t(Indiv_missed[,c("Demuxalot","Demuxalot (refined)","Demuxlet","Dropulation","Freemuxlet","ScSplit","Souporcell","Vireo")]))
colnames(matrix) <- NULL


pdf(paste0(outdir, "NumberIndividualsMissed_heatmap.pdf"), width = 14, height = 2.5)
draw(Heatmap(matrix, name = "Number Individuals Missed", 
  col=c("#BABABA", "#B56A73", "#B21F2C"),
  row_order = c("Demuxalot","Demuxalot (refined)","Demuxlet","Dropulation","Freemuxlet","ScSplit","Souporcell","Vireo"),
  column_order = 1:ncol(matrix),
  bottom_annotation = column_ha,
  column_split = Indiv_missed$Genotyped,
  rect_gp = gpar(col = "grey", lwd = 1),
  column_title = "Pools", 
  column_title_side = "bottom",
  row_names_side = "left",
  heatmap_legend_param = list(nrow = 1)
  ), merge_legend = TRUE,annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()


##### Make figures for % Agreement #####
### Make list of dataframes that have each pair of softwares and the % agreement between them ###
pct_agree_list <- lapply(results_list_wide, function(x){
	combn_df <- as.data.frame(t(combn(softwares, 2, simplify = TRUE)))
	colnames(combn_df) <- c("Software_1", "Software_2")
	combn_df$Percent_Agreement <- NA
	return(combn_df)
})

pct_agree_list <- lapply(names(pct_agree_list), function(x){
	for (row in 1:nrow(pct_agree_list[[x]])){
		if ((pct_agree_list[[x]][row,"Software_1"] %in% demultiplexing_softwares & pct_agree_list[[x]][row,"Software_2"] %in% demultiplexing_softwares)){
      var1 = paste0(pct_agree_list[[x]][row,"Software_1"], "_Assignment")
      var2 = paste0(pct_agree_list[[x]][row,"Software_2"], "_Assignment")
			pct_agree_list[[x]][row, "Percent_Agreement"] <- length(which(results_list_wide[[x]][, ..var1] == results_list_wide[[x]][, ..var2]))/nrow(results_list_wide[[x]])
		} else {
      var1 = paste0(pct_agree_list[[x]][row,"Software_1"], "_DropletType")
      var2 = paste0(pct_agree_list[[x]][row,"Software_2"], "_DropletType")
			pct_agree_list[[x]][row, "Percent_Agreement"] <- length(which(results_list_wide[[x]][, ..var1] == results_list_wide[[x]][, ..var2]))/nrow(results_list_wide[[x]])
		}
	}
	pct_agree_list[[x]]$Pool <- x
	return(pct_agree_list[[x]])
})
names(pct_agree_list) <- names(results_list_wide)


pct_agree_long_df <- do.call(rbind, pct_agree_list)

pct_agree_long_df %>% group_by(., Software_1, Software_2) %>% summarize(., mean = mean(Percent_Agreement))

fwrite(pct_agree_long_df, paste0(outdir, "Percent_agreement_long.tsv"), sep = "\t")

pct_agree_summary <- ddply(pct_agree_long_df, c("Software_1", "Software_2"), summarise,
               N    = length(Percent_Agreement),
               mean = mean(Percent_Agreement, na.rm=TRUE),
               sd   = sd(Percent_Agreement, na.rm=TRUE),
               se   = sd / sqrt(N)
)


pct_agree_summary$Software_1 <- factor(pct_agree_summary$Software_1, levels = softwares)
pct_agree_summary$Software_2 <- factor(pct_agree_summary$Software_2, levels = rev(softwares))

pct_agree_summary_opposite <- pct_agree_summary[,c("Software_2", "Software_1", "N", "mean", "sd", "se")]
colnames(pct_agree_summary_opposite) <- c("Software_1", "Software_2", "N", "mean", "sd", "se")

pct_agree_summary_combined <- rbind(pct_agree_summary,pct_agree_summary_opposite)

for (software in softwares){
	row_n <- nrow(pct_agree_summary_combined) + 1
	pct_agree_summary_combined[row_n,"Software_1"] <- software
	pct_agree_summary_combined[row_n,"Software_2"] <- software
	pct_agree_summary_combined[row_n,"N"] <- 74
	pct_agree_summary_combined[row_n,"mean"] <- 1
	pct_agree_summary_combined[row_n,"sd"] <- 0
	pct_agree_summary_combined[row_n,"se"] <- 0
}

for (software in softwares){
  pct_agree_summary_combined$Software_1 <- gsub(paste0("^",software, "$"), names(softwares[softwares == software]), pct_agree_summary_combined$Software_1)
  pct_agree_summary_combined$Software_2 <- gsub(paste0("^",software, "$"), names(softwares[softwares == software]), pct_agree_summary_combined$Software_2)
}



fwrite(pct_agree_summary_combined, paste0(outdir, "Percent_agreement.tsv"), sep = "\t")
# pct_agree_summary_combined <- fread(paste0(outdir, "Percent_agreement.tsv"), sep = "\t")


pct_agree_summary_combined$Software_1 <- factor(pct_agree_summary_combined$Software_1, levels = names(softwares))
pct_agree_summary_combined$Software_2 <- factor(pct_agree_summary_combined$Software_2, levels = rev(names(softwares)))



heatmap_plot <- ggplot(pct_agree_summary_combined, aes(Software_1, Software_2, fill=mean)) +
	geom_tile() + 
	labs(x = NULL, y = NULL, fill = "Proportion\nAgreement", title="Percent Droplet Assignment\nAgreement Between Softwares") + 
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
		plot.title = element_text(hjust = 0.5)) +
	scale_fill_gradient(low="#FBFEF9",high="#A63446", limits=c(0.69,1)) 

ggsave(heatmap_plot, filename = paste0(outdir, "Percent_agreement.png"), width = 4.7, height = 4)
ggsave(heatmap_plot, filename = paste0(outdir, "Percent_agreement.eps"), width = 4.7, height = 4)
ggsave(heatmap_plot, filename = paste0(outdir, "Percent_agreement.pdf"), width = 4.7, height = 4)



##### Make a dataframe of singlets #####
## Barcode, genotype ID ##
results_singlets <- lapply(results_list_wide, function(x){
  temp <- x[Singlets == 15 &
              demuxlet_Assignment == freemuxlet_Assignment & 
              demuxlet_Assignment == scSplit_Assignment &
              demuxlet_Assignment == souporcell_Assignment &
              demuxlet_Assignment == demuxalot_Assignment &
              demuxlet_Assignment == demuxalot_refined_Assignment &
              demuxlet_Assignment == dropulation_Assignment,c("Barcode","demuxlet_Assignment")]
  colnames(temp) <- c("Barcode","Individual")
  return(temp)
})

## Make a list of individuals in list of pools that has just barcode IDs as DF ##
results_singlets_individuals_list <- lapply(results_singlets, function(x){
  temp <- list()
  for (ind in unique(x$Individual)){
    temp[[ind]] <- x[which(x$Individual == ind),"Barcode"]
  }
  return(temp)
})

lapply(names(results_singlets_individuals_list), function(x){
  lapply(names(results_singlets_individuals_list[[x]]), function(y){
    write_delim(results_singlets_individuals_list[[x]][[y]], delim = "\t", path = paste0(out,x,"/",x,"_",y,"_droplet_barcodes.tsv"), col_names = FALSE)
  })
})




lapply(names(results_singlets), function(x){
  fwrite(results_singlets[[x]], paste0(out,x,"/",x,"individuals.tsv"), col.names = FALSE, sep = "\t")
})

indiv_pool_list <- lapply(names(results_singlets), function(x){
    data.table(Pool = x, Individuals = unique(results_singlets[[x]]$Individual))
})

indiv_pool <- do.call(rbind, indiv_pool_list)
indiv_pool <- indiv_pool[!is.na(Individuals)]

fwrite(indiv_pool, paste0(out,"/Pool_individuals_meta.tsv"), sep = "\t")




##### Make a boxplot figure of the number of singlets per individual per pool #####
results_singlets_4df <- lapply(names(results_singlets), function(x){
  if (nrow(results_singlets[[x]]) > 0){
    results_singlets[[x]] <- as.data.frame(table(results_singlets[[x]][,c("Individual")]))
    colnames(results_singlets[[x]]) <- c("Individual","Count")
    results_singlets[[x]]$Pool <- x
  } else {
    results_singlets[[x]] <-NULL
  }
  return(results_singlets[[x]])
})

results_singlets_df <- do.call(rbind,results_singlets_4df)
results_singlets_df <- left_join(results_singlets_df, meta, by = c("Pool"))
results_singlets_df$N_genotyped <- NA
for (row in 1:nrow(results_singlets_df)){
  results_singlets_df$N_genotyped[row] <- (str_count(results_singlets_df$Individuals[row], pattern = ",") + 1)
}

results_singlets_df$Individuals <- NULL
results_singlets_df$Barcodes <- NULL
results_singlets_df <- results_singlets_df[order(pull(results_singlets_df,N_genotyped), pull(results_singlets_df,Pool)),]
results_singlets_df$Pool  <- gsub("OneK1K_scRNA_Sample","Pool", results_singlets_df$Pool)
results_singlets_df$Pool <- factor(results_singlets_df$Pool, levels = unique(results_singlets_df$Pool))

# ### Get the number of singlets per individual average ###
mean(results_singlets_df$Count) # 720.2254
mean(results_singlets_df[!(results_singlets_df$Pool %in% c("Pool40", "Pool48" ,"Pool66")),]$Count) # 731.96



##### Make a dataframe with doublets for each set of singlets #####
Ndoublets_dt <- data.table("Ntotal" = seq(1:150000), "Ndoublets" = as.numeric(NA), "Nsinglets" = as.numeric(NA))
Ndoublets_dt$Ndoublets <- apply(Ndoublets_dt, 1, function(x) {
  print(x['Ntotal'])
  return(round(min((0.2*x['Ntotal']),(((x['Ntotal']^2)*0.008)/1000))))
})
Ndoublets_dt$Nsinglets <- Ndoublets_dt$Ntotal - Ndoublets_dt$Ndoublets


##### Make new pools #####
##### Make new pools for estimating large numbers #####
if (!file.exists(paste0(out,"/SimulatedPoolsbams_increasingSizes.tsv"))){
  results_singlets_df$Pool  <- gsub("Pool","OneK1K_scRNA_Sample", results_singlets_df$Pool)
  results_singlets_df <- results_singlets_df[which(results_singlets_df$Pool != "OneK1K_scRNA_Sample40"),]
  results_singlets_df <- results_singlets_df[which(results_singlets_df$Pool != "OneK1K_scRNA_Sample48"),]
  results_singlets_df <- results_singlets_df[which(results_singlets_df$Pool != "OneK1K_scRNA_Sample66"),]


  results_singlets[["OneK1K_scRNA_Sample40"]] <- NULL
  results_singlets[["OneK1K_scRNA_Sample48"]] <- NULL
  results_singlets[["OneK1K_scRNA_Sample77"]] <- NULL
  results_singlets_barcodes_df <- do.call(rbind, results_singlets)

  results_singlets_barcodes_df$Barcode <- paste0(gsub("1", "", results_singlets_barcodes_df$Barcode), results_singlets_barcodes_df$Individual)
  results_singlets_barcodes_df <- data.table(results_singlets_barcodes_df)

 ### Get barcode list for each combination and add individual ids ###


  sizes <- list()
  pools <- list()
  bams <- list()
  barcodes <- list()
  Ndoublets <- list()
  Ntotal <- list()


  for (number in 1:7){
    i <- 2^number
    sizes[[paste0("size",i)]] <- as.data.frame(matrix(nrow=10, ncol=i))
    pools[[paste0("pool",i)]] <- as.data.frame(matrix(nrow=10, ncol=i))
    bams[[paste0("bam",i)]] <- as.data.frame(matrix(nrow=10, ncol=i))
    Ndoublets[[paste0("barcode",i)]] <- as.data.frame(matrix(nrow=10, ncol=1))
    colnames(Ndoublets[[paste0("barcode",i)]]) <- "Ndoublets"
    Ntotal[[paste0("barcode",i)]] <- as.data.frame(matrix(nrow=10, ncol=1))
    colnames(Ntotal[[paste0("barcode",i)]]) <- "Ntotal"

    for (rep in 1:10){
      sizes[[paste0("size",i)]][rep,] <- results_singlets_df[sample(1:nrow(results_singlets_df), i),"Individual"]
      for (ind in 1:ncol(sizes[[paste0("size",i)]])){
        pools[[paste0("pool",i)]][rep,ind] <- results_singlets_df[which(results_singlets_df$Individual == sizes[[paste0("size",i)]][rep,ind]),"Pool"]
        bams[[paste0("bam",i)]][rep,ind] <- paste0(out,pools[[paste0("pool",i)]][rep,ind],"/", sizes[[paste0("size",i)]][rep,ind],"_updated.bam")
      }
      Ntotal[[paste0("barcode",i)]][rep,"Ntotal"] <- sum(results_singlets_df[which(results_singlets_df$Individual %in% sizes[[paste0("size",i)]][rep,]),]$Count)
      Ndoublets[[paste0("barcode",i)]][rep,"Ndoublets"] <- round(min(0.2*(sum(results_singlets_df[which(results_singlets_df$Individual %in% sizes[[paste0("size",i)]][rep,]),]$Count)/0.8), (Ndoublets_dt[Nsinglets == sum(results_singlets_df[which(results_singlets_df$Individual %in% sizes[[paste0("size",i)]][rep,]),]$Count)]$Ndoublets[1]) ))
      barcodes[[paste0("size", i, "_SimulatedPool",rep)]] <- results_singlets_barcodes_df[Individual %in% unlist(sizes[[paste0("size",i)]][rep,]), "Barcode"]
    }
  }

  ### check doulbet calculations
  lapply(names(Ntotal), function(x) cbind(Ntotal[[x]], Ndoublets[[x]]))


  ##### Combine into a comma-separated list, the individuals, bams, barcodes and pools for reads from each individual to be pulled from
  sizes_combined <- lapply(sizes, function(x){
    unite(x, sep = ",", col = "Individuals")
  })

  pools_combined <- lapply(pools, function(x){
    unite(x, sep = ",", col = "Pools")
  })

  bams_combined <- lapply(bams, function(x){
    unite(x, sep = ",", col = "Bams")
  })

  # barcodes_combined <- lapply(barcodes, function(x){
  #   unite(x, sep = ",", col = "Barcodes")
  # })


  combined <- lapply(gsub("size","", names(sizes_combined)), function(x){
    temp <- data.table(
        Pool =  paste0("size", x, "_SimulatedPool",1:nrow(sizes_combined[[paste0("size",x)]])),
        size = paste0("size", x),
        N = as.numeric(as.character(x)),
        Ndoublets[[paste0("barcode",x)]],
        sizes_combined[[paste0("size",x)]],
        bams_combined[[paste0("bam",x)]]
    )
    # temp <- cbind(sizes_combined[[paste0("size",x)]], pools_combined[[paste0("pool",x)]], bams_combined[[paste0("bam",x)]], barcodes_combined[[paste0("barcode",x)]])
    # temp$
    # temp$dir <- paste0(temp$size, "_", temp$Pool)
    return(temp)
  })


  sizeCOMBINED <- do.call(rbind, combined)

  ### Add column with N of doublets to add  (Ndoublets) ###
  write_delim(sizeCOMBINED, paste0(out,"SimulatedPoolsbams_increasingSizes.tsv"), delim = "\t")

  ### Also write out barcodes4simulation for each pool (remember to replace -1 with -indiv) ###
  lapply(names(barcodes), function(pool){
    path=paste0(out, pool, "/")
    dir.create(path)
    fwrite(barcodes[[pool]], paste0(path, "barcodes4simulation.tsv"), col.names = FALSE)
  })

}