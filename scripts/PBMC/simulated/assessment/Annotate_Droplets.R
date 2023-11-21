library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(data.table)
library(ggpubr)
library(colorspace)
library(ComplexUpset)
library(Seurat)



save_figs <- function(plot, basename, width = 17, height = 17, units = "cm"){
    ggsave(plot, filename = paste0(basename,".png"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".pdf"), height = height, width = width, units = units)
    ggsave(plot, filename = paste0(basename,".eps"), height = height, width = width, units = units)
}

##### Set up directories #####
demuxafy_dir <- "/path/to/Demuxafy_manuscript/" ## path that has the files from Demuxafy_manuscript from zenodo (recommended) or github
sim_dir <- "/path/to/output/PBMC/simulation/"
dir <- paste0(sim_dir,"/SimulatedOverlap/")
anno_dir <- paste0(sim_dir,"/../Round1Overlap/CellClassification_noRefDoublets/")
outdir <- paste0(dir,"/DropletAnnotation/")
dir.create(outdir)



##### Set up universal 
demultiplexing_list <- c("demuxalot", "demuxalot_refined","demuxlet", "dropulation", "freemuxlet","scSplit","souporcell","vireo")
names(demultiplexing_list) <- c("Demuxalot", "Demuxalot (refined)","Demuxlet", "Dropulation", "Freemuxlet","ScSplit","Souporcell","Vireo")
doublet_detection_list <- c("DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder", "scDblFinder_known_doublets", "scds",  "scrublet", "solo")
names(doublet_detection_list) <- c("DoubletDecon", "DoubletDetection","DoubletFinder",  "ScDblFinder", "ScDblFinder (known doublets)", "Scds", "Scrublet", "Solo")
softwares <- c(demultiplexing_list, doublet_detection_list)

softwares_dt <- data.table(software = softwares, Software = names(softwares))


software_colors <- c("#575E57", "#008751", "#5DBB4D", "#EFE305", "#F8BE33", "#F97919", "#D21F09", "#8D132B", "#F980BF", "#EAA1D1", "#BAAAFF", "#9B64E0", "#6B82E3", "#63D4FE", "#39BCBC", "#C4E4DF")
names(software_colors) <- names(softwares)


##### Read in metadata file #####
sim_pool_meta <- read_delim(paste0(sim_dir,"/SimulatedPoolsbams_increasingSizes.tsv"), delim = "\t")

##### Pull list of pools from metadata files #####
pools <- sim_pool_meta$Pool



##### Read in simulated barcode files #####
sim_barcode_list <- lapply(pools, function(x){
    tmp <- fread(paste0(sim_dir,x,"/matrix_out/barcodes.tsv.gz"), sep = "\t", col.names = "Simulated_Barcode", header = FALSE)
	tmp[, c("Droplet1", "Droplet2") := tstrsplit(Simulated_Barcode, ":", fixed=TRUE)]
	tmp[, c("Barcode1", "Individual1") := tstrsplit(Droplet1, "-", fixed=TRUE)]
	tmp[, c("Barcode2", "Individual2") := tstrsplit(Droplet2, "-", fixed=TRUE)]
	return(tmp)
})
names(sim_barcode_list) <- pools



##### Make list of dataframes that has the the individusal id, the origianl pool name and the individual number as columns #####
individual_conversion_list <- lapply(pools, function(x){
    indivs <- unlist(strsplit(sim_pool_meta[which(sim_pool_meta$Pool == x),]$Individuals, split = ","))
    orig_pools <- gsub("/.+", "", gsub(paste0(sim_dir,"/Round1Overlap/"),"",unlist(strsplit(sim_pool_meta[which(sim_pool_meta$Pool == x),]$Bams, split = ","))))
    data.table("Individual" = indivs, "Original_Pool" = orig_pools)
})

names(individual_conversion_list) <- pools


##### Left join the sim barcode files with the individual_conversion_list
sim_barcode_list <- lapply(pools, function(x){
    tmp <- individual_conversion_list[[x]][sim_barcode_list[[x]], on = c("Individual" = "Individual2")]
    tmp <- individual_conversion_list[[x]][tmp, on = c("Individual" = "Individual1")]
	colnames(tmp) <- gsub("i.Original_Pool", "Original_Pool2", gsub("i.Individual", "Individual2", gsub("^Original_Pool$", "Original_Pool1", gsub("^Individual$", "Individual1", colnames(tmp)))))
	return(tmp)
})
names(sim_barcode_list) <- pools


##### Read in annotate cell info #####
original_pools_meta <- list.files(anno_dir, pattern = "_metadata_no_doublet.RDS")
original_pools <- gsub("_metadata_no_doublet.RDS", "",original_pools_meta)

cell_annotations <- lapply(original_pools_meta, function(x){
    readRDS(paste0(anno_dir, x))
})

names(cell_annotations) <- original_pools


cell_annotations_sub <- lapply(cell_annotations, function(x){
    df <- x[,c("Barcode","predicted.celltype.l2")]
    df$Pool <- gsub("_[ATCG]+-1","",rownames(df))
    return(df)
})

cell_annotations_sub_df <- data.table(do.call(rbind, cell_annotations_sub))
cell_annotations_sub_df$Barcode <- gsub("-1", "", cell_annotations_sub_df$Barcode)



##### Add annotations to sim_barcode_list #####
sim_barcode_list_annos <- lapply(sim_barcode_list, function(x){
    tmp <- cell_annotations_sub_df[x, on = c("Barcode" = "Barcode1", "Pool" = "Original_Pool1")]
	tmp <- cell_annotations_sub_df[tmp, on = c("Barcode" = "Barcode2", "Pool" = "Original_Pool2")]
	colnames(tmp) <- gsub("^Barcode$", "Barcode2", gsub("^Pool$", "Pool2",gsub("i.Pool", "Pool1", gsub("i.Barcode", "Barcode1", gsub("i.CellType1", "CellType1", gsub("^CellType1$", "CellType2", colnames(tmp)))))))
	return(tmp)
})
names(sim_barcode_list_annos) <- pools


lapply(names(sim_barcode_list), function(x){
    nrow(sim_barcode_list[[x]]) - nrow(sim_barcode_list_annos[[x]])
})



##### Add if homotypic or heterotypic doublet and if from same or different individuals ###### 
sim_barcode_list_annos_doublets <- lapply(sim_barcode_list_annos, function(x){
    x$doublet_singlet <- ifelse(is.na(x$Barcode2), "singlet","doublet")
    x$heterotypic_homotypic <- ifelse((x$doublet_singlet == "doublet" & x$CellType1 == x$CellType2), "Homotypic",
                                        ifelse((x$doublet_singlet == "doublet" & x$CellType1 != x$CellType2), "Heterotypic", NA))
    x$individuals <- ifelse((x$doublet_singlet == "doublet" & x$Individual1 == x$Individual2), "Homogenic", 
                            ifelse((x$doublet_singlet == "doublet" & x$Individual1 != x$Individual2), "Heterogenic", NA))
    return(x)
})



### Identified an issue with some doublets not having annotations for one of the cells
### This turned out to be an issue with the droplets used for doublet making and were pools that were rerun with the droplet only in the original pool not used in the simulated pool.
### I can't figure out HOW this happened but it's clear that those doublets are not true doublets because the second cell didn't exist in the bam so it's just the singlet of the first cell but labeled as a doublet
### Make a list of the droplets that this happened to in each simulated pool => use to remove before metric calculations in best methods scripts
### Also remove here

barcodes2remove <- lapply(sim_barcode_list_annos_doublets, function(x){
	barcodes <- x[is.na(CellType1) & doublet_singlet == "doublet"]$Simulated_Barcode
	barcodes <- c(barcodes,x[is.na(CellType2) & doublet_singlet == "doublet"]$Simulated_Barcode)
	return(barcodes)
})
names(barcodes2remove) <- names(sim_barcode_list_annos_doublets)

saveRDS(barcodes2remove, paste0(outdir,"barcodes2remove.rds"))
barcodes2remove <- readRDS(paste0(outdir,"barcodes2remove.rds"))


sim_barcode_list_annos_doublets <- lapply(names(sim_barcode_list_annos_doublets), function(pool){
	sim_barcode_list_annos_doublets[[pool]][!(Simulated_Barcode %in% barcodes2remove[[pool]])]
})
names(sim_barcode_list_annos_doublets) <- names(barcodes2remove)
lapply(sim_barcode_list_annos_doublets, dim)


saveRDS(sim_barcode_list_annos_doublets, paste0(outdir,"simulated_droplet_types.rds"))


sim_barcode_list_annos_doublets <- readRDS(paste0(outdir,"simulated_droplet_types.rds"))




##### Make Singlets and Doublets per N figure #####
prop_df_sing_doub_list <- lapply(names(sim_barcode_list_annos_doublets), function(x){
    temp <- as.data.frame(prop.table(table(sim_barcode_list_annos_doublets[[x]]$doublet_singlet)))
    colnames(temp) <- c("doublet_singlet", "Proportion")
    temp$Pool_Size <- as.numeric(as.character(gsub("size","",x) %>% gsub("_SimulatedPool\\d+", "", .)))
	temp$N <- nrow(sim_barcode_list_annos_doublets[[x]])
	temp$N_doublets <- length(which(!is.na(sim_barcode_list_annos_doublets[[x]]$Individual2)))
	temp$N_singlets <- length(which(is.na(sim_barcode_list_annos_doublets[[x]]$Individual2)))
    return(temp)
})

prop_df_sing_doub <- data.table(do.call(rbind, prop_df_sing_doub_list))
prop_df_sing_doub$Pool_Size <- factor(prop_df_sing_doub$Pool_Size)

prop_df_sing_doub$Percentage <- prop_df_sing_doub$Proportion * 100

mean_N <- data.table(Pool_Size = unique(prop_df_sing_doub$Pool_Size))
mean_N$mean <- as.character(NA)

for (n in mean_N$Pool_Size){
	mean_N[Pool_Size == n]$mean <- formatC(round(mean(unique(prop_df_sing_doub[Pool_Size == n,c("Pool_Size", "N", "N_doublets", "N_singlets")])$N)), format="d", big.mark=",")
}


p_doub_sing <- ggbarplot(prop_df_sing_doub, 
							x = "Pool_Size", 
							y = "Percentage", 
							add = "mean_se", 
							fill = "doublet_singlet",
							color = "black",
							palette = c("#6a38b3", "#81B338"),
							xlab = "Number Individuals\nMultiplexed",
							ylab = "Percentage of Droplets",
							legend = "right"
							) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							geom_text(data = mean_N, aes(x = Pool_Size, y = 103, label = mean), angle = 90, vjust = 0.5, hjust=0) +
							scale_y_continuous(expand = c(0, 0), breaks = c(0, 25, 50, 75, 100), limits = c(0, 130))

ggsave(p_doub_sing, filename = paste0(outdir, "Singlet_Doublet_Proportions.png"), width = 3.4, height = 4)
ggsave(p_doub_sing, filename = paste0(outdir, "Singlet_Doublet_Proportions.pdf"), width = 3.4, height = 4)



##### Visualize doublets separately for each  #####
lapply(sim_barcode_list_annos_doublets, function(x){
    as.data.frame(prop.table(table(x$heterotypic_homotypic, x$individuals)))
})

celltypes <- unique(unlist(lapply(sim_barcode_list_annos_doublets, function(x){
	unique(x$CellType1)
})))

sim_barcode_list_annos_table <- lapply(names(sim_barcode_list_annos_doublets), function(x){
	sim_barcode_list_annos_doublets[[x]]$CellType_Doublet <- ifelse(sim_barcode_list_annos_doublets[[x]]$doublet_singlet == "doublet", paste0("Doublet: ", sim_barcode_list_annos_doublets[[x]]$heterotypic_homotypic, " & ", sim_barcode_list_annos_doublets[[x]]$individuals), sim_barcode_list_annos_doublets[[x]]$CellType1)
	tbl <- data.table(table(sim_barcode_list_annos_doublets[[x]]$CellType_Doublet))
	tmp1 <- t(tbl$N)
	colnames(tmp1) <- tbl$V1
	temp2 <- data.table(Pool = x, "Number Multiplexed" = gsub("size", "", gsub("_SimulatedPool\\d+", "", x)), "Number Droplets" = sum(as.numeric(tbl$N)), tmp1)
	for (droplet in c("Doublet: Heterotypic & Heterogenic", "Doublet: Heterotypic & Homogenic", "Doublet: Homotypic & Homogenic", "Doublet: Homotypic & Heterogenic",celltypes[!(celltypes %in% colnames(temp2))])){
		if (!(droplet %in% colnames(temp2))){
			print(droplet)
			dt <- data.table(celltype = 0)
			colnames(dt) <- droplet
			temp2 <- cbind(temp2,dt)
		}
	}
	no_doublet_cols <- c("Pool", "Number Multiplexed", "Number Droplets", celltypes, "Doublet: Heterotypic & Heterogenic", "Doublet: Heterotypic & Homogenic", "Doublet: Homotypic & Homogenic", "Doublet: Homotypic & Heterogenic") 
	setcolorder(temp2,no_doublet_cols)
    return(temp2)
})

sim_barcode_list_annos_table_dt <- do.call(rbind, sim_barcode_list_annos_table)


fwrite(sim_barcode_list_annos_table_dt, paste0(outdir, "sim_pool_characterisation.tsv"), sep = "\t")



##### Make a bar graph of the hetertypic vs homotypic (colored by same or different individual) and split by size of pool #####
prop_df_list <- lapply(names(sim_barcode_list_annos_doublets), function(x){
    temp <- as.data.frame(prop.table(table(sim_barcode_list_annos_doublets[[x]]$heterotypic_homotypic, sim_barcode_list_annos_doublets[[x]]$individuals)))
    colnames(temp) <- c("heterotypic_homotypic", "heterogenic_homogenic", "Proportion")
    temp$Pool_Size <- as.numeric(as.character(gsub("size","",x) %>% gsub("_SimulatedPool\\d+", "", .)))
	temp$N_doublets <- length(which(!is.na(sim_barcode_list_annos_doublets[[x]]$Individual2)))
	temp$N_singlets <- length(which(is.na(sim_barcode_list_annos_doublets[[x]]$Individual2)))
	temp$N_doublet_het_homo <- temp$N_doublets * temp$Proportion
	temp$N <- round(data.table(table(sim_barcode_list_annos_doublets[[x]]$heterotypic_homotypic, sim_barcode_list_annos_doublets[[x]]$individuals))$N)
    return(temp)
})


prop_df <- data.table(do.call(rbind, prop_df_list))
prop_df$Pool_Size <- factor(prop_df$Pool_Size)

prop_df$Hetero_combined <- paste0(prop_df$heterotypic_homotypic, "-", prop_df$heterogenic_homogenic)
prop_df$Hetero_combined <- factor(prop_df$Hetero_combined, levels = rev(c("Heterotypic-Heterogenic", "Heterotypic-Homogenic", "Homotypic-Heterogenic", "Homotypic-Homogenic")))
mean_prop_df <- unique(prop_df[,.(Mean = round(mean(N), 0)),.(heterotypic_homotypic,heterogenic_homogenic, Pool_Size)])


doublet_types <- ggplot(prop_df, aes(Pool_Size, Proportion, fill = heterogenic_homogenic)) + 
  geom_bar(position = "stack", stat = "summary") +
  facet_wrap(vars(heterotypic_homotypic), nrow = 1) +
  theme_classic() +
  scale_fill_manual(values = c("#6A51A3", "#238B45")) +
  labs(fill = "Individuals Doublets\nCreated From") +
  xlab("Number Individuals per Pool")

ggsave(doublet_types, filename = paste0(outdir, "doublet_type_proportions.png"))


prop_df$Percentage <- prop_df$Proportion * 100

doublet_types_stacked <- ggbarplot(prop_df, 
									x = "Pool_Size", 
									y = "Percentage", 
									add = "mean_se", 
									fill = "Hetero_combined",
									color = "black",
									palette = c("grey", "#d3212d", "#3c50b1", "#6a38b3"),
									xlab = "Number Individuals\nMultiplexed",
									ylab = "Percentage of Doublets",
									legend = "right"
									) +
									theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
									scale_y_continuous(expand = c(0, 0))


ggsave(doublet_types_stacked, filename = paste0(outdir, "stacked_hetero_doublets.png"), width = 4.1, height = 3.3)
ggsave(doublet_types_stacked, filename = paste0(outdir, "stacked_hetero_doublets.pdf"), width = 4.1, height = 3.3)



prop_df$Hetero_combined <- factor(prop_df$Hetero_combined, levels = rev(levels(prop_df$Hetero_combined)))

doublet_types_numbers <- ggbarplot(prop_df, 
									x = "Hetero_combined", 
									y = "N", 
									facet.by = "Pool_Size",
									nrow = 1,
									label = TRUE,
									label.pos = "out",
									scales = "free_x",
									add = "mean_se", 
									fill = "Hetero_combined",
									color = "black",
									palette = rev(c("grey", "#d3212d", "#3c50b1", "#6a38b3")),
									xlab = NULL,
									ylab = "Number of Droplets",
									legend = "right"
									) + coord_flip() +
									rotate_x_text(90) 


ggsave(doublet_types_numbers, filename = paste0(outdir, "hetero_doublets_numbers.png"), width = 9, height = 3)
ggsave(doublet_types_numbers, filename = paste0(outdir, "hetero_doublets_numbers.pdf"), width = 8.3, height = 3)




### Make a figure using the average percentages to get 20,000 cells for different # individuals multiplexed ###
summarized_prop_df <- prop_df[,.(Mean=mean(Proportion)),.(Pool_Size,heterotypic_homotypic,heterogenic_homogenic,Hetero_combined)]
summarized_prop_df$droplets_20k <- round(20000 * 0.16 * summarized_prop_df$Mean,0)
summarized_prop_df$Hetero_combined <- factor(summarized_prop_df$Hetero_combined, levels = c("Heterotypic-Heterogenic", "Heterotypic-Homogenic", "Homotypic-Heterogenic", "Homotypic-Homogenic"))


doublet_types_numbers_20k <- ggbarplot(summarized_prop_df[Pool_Size %in% c(2,4,8,16,32)], 
									x = "Pool_Size", 
									y = "droplets_20k", 
									position = position_dodge2(), 
									# facet.by = "Pool_Size",
									nrow = 1,
									label = TRUE,
									label.pos = "out",
									scales = "free_x",
									add = "mean_se", 
									fill = "Hetero_combined",
									color = "black",
									palette = rev(c("grey", "#d3212d", "#3c50b1", "#6a38b3")),
									xlab = "Number Multiplexed Individuals",
									ylab = "Number of Droplets",
									legend = "right"
									)

ggsave(doublet_types_numbers_20k, filename = paste0(outdir, "hetero_doublets_numbers_20k.png"), width = 7, height = 3)
ggsave(doublet_types_numbers_20k, filename = paste0(outdir, "hetero_doublets_numbers_20k.pdf"), width = 7, height = 3)


##### Try an upset plot #####
### Prepare dataframe ###
upset_df <- prop_df
colnames(upset_df) <- c("Heterotypic", "Heterogenic", "Proportion","Pool_Size")
upset_df$Heterotypic <- as.numeric(gsub("Heterotypic", 1, upset_df$Heterotypic) %>% gsub("Homotypic", 0, .))
upset_df$Heterogenic <- as.numeric(gsub("different", 1, upset_df$Heterogenic) %>% gsub("same", 0, .))

upset_data(upset_df,c("Heterogenic", "Heterotypic"))


pUPSET <- upset(upset_df,
	c("Heterogenic", "Heterotypic"), 
	set_sizes=FALSE, 
    base_annotations=list(
        "Proportion" = (
                    ggplot(upset_df, aes( y = Proportion, fill = Pool_Size)) +
                      geom_bar(stat="identity", position=position_dodge()) +
                      theme_classic() +
					  theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        plot.title = element_text(hjust = 0.5)) +
                      ylab(NULL)
                  )
    ),
    width_ratio=0.1
)

pdf(file= paste0(outdir,"Heterotypic_Heterogenic_upset.pdf"), height = 4, width = 6) # or other device
pUPSET
dev.off()

### Want the proportion of each N in each category, make ggplot to replace the counts in the upset plot

##### Make figures comparing the cell type to singlet doublet classifications ######
singlets_list <- lapply(pools, function(x){
	print(x)
	readRDS(paste0(dir, x, "/singlet_counts_per_barcode_single_soft.rds"))
})
names(singlets_list) <- pools

singlets_sim_barcode_list_annos_doublets_list <- lapply(names(singlets_list), function(x){
	tmp <- left_join(sim_barcode_list_annos_doublets[[x]], singlets_list[[x]], by = c("Simulated_Barcode" = "Barcode"))
	tmp[,names(which(colSums(tmp[,..softwares]) == 0))] <- NA
	return(tmp)
})
names(singlets_sim_barcode_list_annos_doublets_list) <- names(singlets_list)


soft_classification_annotation_list <- lapply(names(singlets_sim_barcode_list_annos_doublets_list), function(x){
	temp2 <- lapply(colnames(singlets_list[[1]])[2:ncol(singlets_list[[1]])], function(y){
		temp <- data.table(matrix(nrow = nrow(singlets_sim_barcode_list_annos_doublets_list[[x]]), ncol = 3))
		colnames(temp) <- c("Barcode" ,"Reference", "Predicted")
		temp$Barcode <- singlets_sim_barcode_list_annos_doublets_list[[x]][,"Simulated_Barcode"]
		temp$Reference <- ifelse(singlets_sim_barcode_list_annos_doublets_list[[x]]$doublet_singlet == "singlet", singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1, paste0("doublet_",singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals))
		temp$Predicted <- ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 1 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "singlet"), singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1,
								ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 1 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "doublet"), "doublet > singlet",
									ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 0 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "doublet"), paste0("doublet_", singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals), "singlet > doublet")))
		temp$software <- y
		temp$status <- factor(ifelse(temp$Reference == temp$Predicted, "Correct", temp$Predicted), levels = c("Correct", "doublet > singlet", "singlet > doublet"))


		temp3 <- data.table(prop.table(table(paste0(temp$Reference, "-", temp$software), temp$status), 1))
		colnames(temp3) <- c("Reference_software", "status", "N")
		temp3 <- temp3[, c("Reference", "software") := tstrsplit(Reference_software, "-", fixed = TRUE)]
		temp3$Reference_software <- NULL


		return(temp3)
	})
	temp3 <- do.call(rbind, temp2)
	temp3$pool <- x
	temp3$Number_Individuals <- as.numeric(gsub("size", "", temp3$pool) %>% gsub("_SimulatedPool.+","",.))
	temp3$Additional_Mt_Percent <- ifelse(grepl("mt", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_mt_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	temp3$Additional_Ambient_Percent <- ifelse(grepl("ambient", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_ambient_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	temp3$Downsampled <- ifelse(grepl("subsampled", temp3$pool), TRUE,FALSE)
	temp3$Uneven <- ifelse(grepl("uneven", temp3$pool), TRUE,FALSE)
	return(temp3)
})


soft_classification_annotation <- do.call(rbind, soft_classification_annotation_list)

soft_classification_annotation$software <- factor(soft_classification_annotation$software, levels = softwares)
soft_classification_annotation$Number_Individuals <- factor(soft_classification_annotation$Number_Individuals)


soft_classification_annotation$Reference <- gsub("doublet_", "Doublet:\n", soft_classification_annotation$Reference) %>% 
									gsub("_Heterogenic", " &\nHeterogenic", .) %>%
									gsub("singlet", "Singlet", .) %>%
									gsub("_Homogenic", " &\nHomogenic", .) %>%
									gsub(" Proliferating", "\nProliferating", .) %>%
									gsub("NK_", "NK\n", .) %>%
									gsub(" intermediate", "\nintermediate", .)



soft_classification_annotation$Reference <- factor(soft_classification_annotation$Reference, levels = c(unique(soft_classification_annotation$Reference)[!grepl("Doublet", unique(soft_classification_annotation$Reference))], "Doublet:\nHeterotypic &\nHeterogenic", "Doublet:\nHomotypic &\nHeterogenic", "Doublet:\nHeterotypic &\nHomogenic", "Doublet:\nHomotypic &\nHomogenic"))

soft_classification_annotation$software_type <- ifelse(soft_classification_annotation$software %in% demultiplexing_list, "Demultiplexing\nSoftwares", "Doublet Detecting\nSoftwares")

### Save classification annotations for use in other scripts ###
saveRDS(soft_classification_annotation, paste0(outdir,"single_software_classifications_annotations_df.rds"))
saveRDS(soft_classification_annotation_list, paste0(outdir,"single_software_classifications_annotations_list.rds"))

soft_classification_annotation <- readRDS(paste0(outdir,"single_software_classifications_annotations_df.rds"))
soft_classification_annotation_list <- readRDS(paste0(outdir,"single_software_classifications_annotations_list.rds"))


soft_classification_annotation_known <- soft_classification_annotation[software %in% c("scDblFinder","scDblFinder_known_doublets") & Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & Uneven == FALSE]
soft_classification_annotation_known_wide <- dcast(soft_classification_annotation_known, status + Reference + pool + Number_Individuals + software_type ~ software, value.var = "N")


ttest_known_df <- data.table(unique(soft_classification_annotation[software %in% c("scDblFinder") & Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE,
																			c("status", "Reference", "software","Number_Individuals")]))
ttest_known_df$software2 <- ifelse(ttest_known_df$software == "scDblFinder", "scDblFinder_known_doublets")

ttest_known_df <- ttest_known_df[status == "Correct" ]

ttest_known_df$statistic <- as.numeric(NA)
ttest_known_df$p <- as.numeric(NA)
ttest_known_df$mean1 <- as.numeric(NA)
ttest_known_df$mean2 <- as.numeric(NA)

for (row in 1:nrow(ttest_known_df)){
	soft1 <- ttest_known_df$software[row]
	soft2 <- ttest_known_df$software2[row]
	celltype <- ttest_known_df$Reference[row]
	group <- ttest_known_df$status[row]
	n <- ttest_known_df$Number_Individuals[row]
	if ((all(is.na(soft_classification_annotation_known_wide[Reference == celltype & status == group & Number_Individuals == n, ..soft1][[1]])) & all(is.na(soft_classification_annotation_known_wide[Reference == celltype & status == status & Number_Individuals == n, ..soft2][[1]]))) | nrow(soft_classification_annotation_known_wide[Reference == celltype & status == group & Number_Individuals == n]) < 3){
		ttest_known_df$p[row] <- NA
		ttest_known_df$statistic[row] <- NA
	} else {
		tryCatch({
			ttest <- t.test(soft_classification_annotation_known_wide[Reference == celltype & status == group & Number_Individuals == n, ..soft1][[1]], soft_classification_annotation_known_wide[Reference == celltype & status == group & Number_Individuals == n,..soft2][[1]])
			ttest_known_df$p[row] <-  ttest$p.value
			ttest_known_df$statistic[row] <- ttest$statistic
			ttest_known_df$mean1[row] <- ttest$estimate[1]
			ttest_known_df$mean2[row] <- ttest$estimate[2]
		},
		error = function(e) {
			ttest_known_df$p[row] <-  NA
			ttest_known_df$statistic[row] <- NA
		})
	}
}


min(ttest_known_df[software == "scDblFinder"]$p, na.rm = TRUE)
max(ttest_known_df[software == "scDblFinder"]$p, na.rm = TRUE)
ttest_known_df[software == "scDblFinder" & p < 0.05]
ttest_known_df[software == "scDblFinder" & !is.na(mean1)]



##### Make figure of stacked bar plot showing proportion of cells in each cell cycle group #####
pNormal_stacked <- ggbarplot(soft_classification_annotation[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE],"Number_Individuals", "N", add = c("mean_se"), fill = "status", facet.by = c("software", "Reference")) +
	rotate_x_text(90) +
	# coord_flip() +
	scale_fill_manual(values = c("#05734A", "#FF9801", "#FDFD96"))

save_figs(pNormal_stacked,  paste0(outdir, "stacked_software_annotations"), width =80, height = 40)


### Make plots with % correct (facet by cell type but put all softwares on same plot)
### Demultiplexing Softwares
pNormal_demultiplexing_doublet_cell_type <- ggplot(soft_classification_annotation[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares"], aes(Number_Individuals, N*100, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_wrap(vars(Reference), nrow = 5) +
	rotate_x_text(90) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	scale_fill_manual(values = software_colors) +
	ylab("Percent Correct") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
		legend.position = "top")

save_figs(pNormal_demultiplexing_doublet_cell_type,  paste0(outdir, "demultiplex_doublet_correct_cell_type"), width =17, height = 24)


### Doublet Detecting Softwares
pNormal_txn_doublet_cell_type <- ggplot(soft_classification_annotation[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & status == "Correct" & software_type == "Doublet Detecting\nSoftwares"], aes(Number_Individuals, N*100, color = software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_wrap(vars(Reference), nrow = 5) +
	rotate_x_text(90) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	scale_fill_manual(values = software_colors) +
	ylab("Percent Correct") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
		legend.position = "top")

save_figs(pNormal_txn_doublet_cell_type,  paste0(outdir, "txn_doublet_correct_cell_type"), width =17, height = 22.5)



### Demultiplexing Software heterogenic doublets homotypic vs heterotypic ###
pNormal_demultiplexing_doublet_typic <- ggplot(soft_classification_annotation[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares"], aes(Number_Individuals, N*100, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_wrap(vars(Reference), nrow = 5) +
	rotate_x_text(90) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	scale_fill_manual(values = software_colors) +
	ylab("Percent Correct") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
		legend.position = "top")

save_figs(pNormal_demultiplexing_doublet_typic,  paste0(outdir, "demultiplex_doublet_correct_cell_type"), width =17, height = 25)



lapply(unique(soft_classification_annotation$software), function(soft){
	pNormal_stacked_list <- ggbarplot(soft_classification_annotation[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & software == soft],"Number_Individuals", "N", add = c("mean_se"), fill = "status", facet.by = c("Reference", "software")) +
		rotate_x_text(90) +
		scale_fill_manual(values = c("#05734A", "#FF9801", "#FDFD96"))

	save_figs(pNormal_stacked_list,  paste0(outdir, soft, "_stacked_software_annotations"), width =5, height = 60)
})





soft_annotation_list <- lapply(names(singlets_sim_barcode_list_annos_doublets_list), function(x){
	temp2 <- lapply(colnames(singlets_list[[1]])[2:ncol(singlets_list[[1]])], function(y){
		temp <- data.table(matrix(nrow = nrow(singlets_sim_barcode_list_annos_doublets_list[[x]]), ncol = 3))
		colnames(temp) <- c("Barcode" ,"Reference", "Predicted")
		temp$Barcode <- singlets_sim_barcode_list_annos_doublets_list[[x]][,"Simulated_Barcode"]
		temp$Reference <- ifelse(singlets_sim_barcode_list_annos_doublets_list[[x]]$doublet_singlet == "singlet", singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1, paste0("doublet_",singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals))
		temp$Predicted <- ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 1 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "singlet"), singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1,
								ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 1 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "doublet"), "doublet > singlet",
									ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 0 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "doublet"), paste0("doublet_", singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals), "singlet > doublet")))
		temp$software <- y
		temp$status <- factor(ifelse(temp$Reference == temp$Predicted, "Correct", temp$Predicted), levels = c("Correct", "doublet > singlet", "singlet > doublet"))
		temp$Reference_status <- ifelse(grepl("doublet", temp$Reference), temp$Reference, "singlet")


		temp3 <- data.table(prop.table(table(paste0(temp$Reference_status, "-", temp$software), temp$status), 1))
		colnames(temp3) <- c("Reference_software", "status", "N")
		temp3 <- temp3[, c("Reference", "software") := tstrsplit(Reference_software, "-", fixed = TRUE)]
		temp3$Reference_software <- NULL

		return(temp3)
	})
	temp3 <- do.call(rbind, temp2)
	temp3$pool <- x
	temp3$Number_Individuals <- as.numeric(gsub("size", "", temp3$pool) %>% gsub("_SimulatedPool.+","",.))
	temp3$Additional_Mt_Percent <- ifelse(grepl("mt", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_mt_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	temp3$Additional_Ambient_Percent <- ifelse(grepl("ambient", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_ambient_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	temp3$Downsampled <- ifelse(grepl("subsampled", temp3$pool), TRUE,FALSE)
	temp3$Uneven <- ifelse(grepl("uneven", temp3$pool), TRUE,FALSE)
	return(temp3)
})



soft_classification <- do.call(rbind, soft_annotation_list)

soft_classification$software <- factor(soft_classification$software, levels = softwares)

soft_classification <- data.table(soft_classification)

soft_classification$software_type <- ifelse(soft_classification$software %in% demultiplexing_list, "Demultiplexing\nSoftwares", "Doublet Detecting\nSoftwares")
soft_classification$Number_Individuals <- factor(soft_classification$Number_Individuals)

soft_classification$Reference <- gsub("doublet_", "Doublet:\n", soft_classification$Reference) %>% 
									gsub("_Heterogenic", " &\nHeterogenic", .) %>%
									gsub("singlet", "Singlet", .) %>%
									gsub("_Homogenic", " &\nHomogenic", .)


fwrite(soft_classification, paste0(outdir, "software_classifications.tsv"), sep = "\t")
soft_classification <- fread(paste0(outdir, "software_classifications.tsv"), sep = "\t")

soft_classification <- softwares_dt[soft_classification, on = "software"]

soft_classification$Reference <- factor(soft_classification$Reference, levels = c("Singlet", "Doublet:\nHeterotypic &\nHeterogenic", "Doublet:\nHomotypic &\nHeterogenic", "Doublet:\nHeterotypic &\nHomogenic", "Doublet:\nHomotypic &\nHomogenic"))


##### Make figure of stacked bar plot showing proportion of cells in each cell cycle group #####
pNormal_stacked <- ggbarplot(soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & Uneven == FALSE],"Number_Individuals", "N", add = c("mean_se"), fill = "status", facet.by = c("software", "Reference")) +
	rotate_x_text(90) +
	scale_fill_manual(values = c("#05734A", "#FF9801", "#FDFD96"))

save_figs(pNormal_stacked,  paste0(outdir, "stacked_software"), width =20, height = 30)


##### Make figure of stacked bar plot showing proportion of cells in each cell cycle group #####
pNormal_demultiplexing_doublet <- ggplot(soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & Uneven == FALSE & status == "Correct"], aes(factor(Number_Individuals), N*100, color = Software)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_grid(software_type ~ Reference) +
	rotate_x_text(90) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	ylab("Percent Correct") +
	geom_smooth(aes(group = Software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

save_figs(pNormal_demultiplexing_doublet,  paste0(outdir, "demultiplex_doublet_correct"), width =25, height = 10)





##### Make figure of showing diffrerence of heterotypic vs homotypic % correct for demultiplexing softwares #####
pNormal_demultiplexing_doublet <- ggplot(soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares" & grepl("Heterogenic", Reference)], aes(Number_Individuals, N*100, color = Reference)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_wrap(vars(software), scales = "free_y", nrow = 1) +
	rotate_x_text(90) +
	theme_classic() +
	scale_color_manual(values = c("#d3212d", "#6a38b3")) +
	scale_fill_manual(values = c("#d3212d", "#6a38b3")) +
	ylab("Percent Correct") +
	geom_smooth(aes(group = Reference), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

save_figs(pNormal_demultiplexing_doublet,  paste0(outdir, "demultiplex_doublet_correct_type_comparison"), width =25, height = 5)


pNormal_demultiplexing_doublet_grid <- ggplot(soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares" & grepl("Heterogenic", Reference)], aes(Number_Individuals, N*100, color = Reference)) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	facet_grid(~software, scales = "free_y") +
	rotate_x_text(90) +
	theme_classic() +
	scale_color_manual(values = c("#d3212d", "#6a38b3")) +
	scale_fill_manual(values = c("#d3212d", "#6a38b3")) +
	ylab("Percent Correct") +
	geom_smooth(aes(group = Reference), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

save_figs(pNormal_demultiplexing_doublet_grid,  paste0(outdir, "demultiplex_doublet_correct_type_comparison_grid"), width =25, height = 5)




### Test for significant differences ###
ttest <- data.table(unique(soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & Uneven == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares" &  Reference == "Doublet:\nHeterotypic &\nHeterogenic", c("Number_Individuals", "software")]))
ttest$pval <- as.numeric(NA)

for(number in unique(soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & Uneven == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares" & grepl("Heterogenic", Reference)]$Number_Individuals)){
	for (soft in unique(soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & Uneven == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares" & grepl("Heterogenic", Reference)]$software)){
		tryCatch({
		ttest[Number_Individuals == number & software == soft]$pval <- t.test(soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & Uneven == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares" & Number_Individuals == number & Reference == "Doublet:\nHeterotypic &\nHeterogenic" & software == soft]$N,
			soft_classification[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & Uneven == FALSE & status == "Correct" & software_type == "Demultiplexing\nSoftwares" & Number_Individuals == number & Reference == "Doublet:\nHomotypic &\nHeterogenic" & software == soft]$N)$p.value
		},
		error = function(e) {
			ttest[Number_Individuals == number & software == soft]$pval <-  NA
		})
	}
}


### Test for difference in variance for demultiplexing vs doublet detecting ###
var_dt <- data.table(unique(soft_classification[,c("Number_Individuals", "Reference")]))
var_dt$p <- as.numeric(NA)

for (row in 1:nrow(var_dt)){
	ref <- var_dt$Reference[row]
	n <- var_dt$Number_Individuals[row]
	ftest <- var.test(soft_classification[status == "Correct" & Reference == ref & software %in% demultiplexing_list & Number_Individuals == n]$N, soft_classification[status == "Correct" & Reference == ref & software %in% doublet_detection_list & Number_Individuals == n]$N, 
         alternative = "two.sided")
	var_dt$p[row] <- ftest$p.value
}

### Get proportions of each cell type in each pool ###
soft_cell_type_props <- lapply(names(singlets_sim_barcode_list_annos_doublets_list), function(x){
	temp <- data.table(matrix(nrow = nrow(singlets_sim_barcode_list_annos_doublets_list[[x]]), ncol = 2))
	colnames(temp) <- c("Barcode" ,"Reference")
	temp$Barcode <- singlets_sim_barcode_list_annos_doublets_list[[x]][,"Simulated_Barcode"]
	temp$Reference <- ifelse(singlets_sim_barcode_list_annos_doublets_list[[x]]$doublet_singlet == "singlet", singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1, paste0("doublet_",singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals))


	temp3 <- data.table(prop.table(table(temp$Reference)))
	colnames(temp3) <- c("Reference", "N")


	temp3$pool <- x
	temp3$Number_Individuals <- as.numeric(gsub("size", "", temp3$pool) %>% gsub("_SimulatedPool.+","",.))
	temp3$Additional_Mt_Percent <- ifelse(grepl("mt", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_mt_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	temp3$Additional_Ambient_Percent <- ifelse(grepl("ambient", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_ambient_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	temp3$Downsampled <- ifelse(grepl("subsampled", temp3$pool), TRUE,FALSE)
	return(temp3)
})

soft_cell_type_props_df <- data.table(do.call(rbind, soft_cell_type_props))

soft_cell_type_props_df_summary <- soft_cell_type_props_df[,.(Mean=mean(N), SD = sd(N)),.(Reference,Number_Individuals)]

soft_cell_type_props_df_summary_summary <- soft_cell_type_props_df[,.(Mean = mean(N), Min = min(N), Max = max(N)),.(Reference)]



### Look at characteristics of false doublets and false singlets ###
soft_correct_incorrcet_list <- lapply(names(singlets_sim_barcode_list_annos_doublets_list), function(x){
	print(x)
	temp2 <- lapply(colnames(singlets_list[[1]])[2:ncol(singlets_list[[1]])], function(y){
		# print(x)
		print(y)
		temp <- data.table(matrix(nrow = nrow(singlets_sim_barcode_list_annos_doublets_list[[x]]), ncol = 3))
		colnames(temp) <- c("Barcode" ,"Reference", "Predicted")
		temp$Barcode <- singlets_sim_barcode_list_annos_doublets_list[[x]][,"Simulated_Barcode"]
		temp$Reference <- ifelse(singlets_sim_barcode_list_annos_doublets_list[[x]]$doublet_singlet == "singlet", singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1, paste0("doublet_",singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals))
		temp$Predicted <- ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 1 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "singlet"), singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1,
								ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 1 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "doublet"), "doublet > singlet",
									ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 0 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "doublet"), paste0("doublet_", singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals), "singlet > doublet")))
		temp$software <- y
		temp$status <- factor(ifelse(temp$Reference == temp$Predicted, "Correct", temp$Predicted), levels = c("Correct", "doublet > singlet", "singlet > doublet"))
		temp$pool <- x
		temp$Number_Individuals <- as.numeric(gsub("size", "", temp$pool) %>% gsub("_SimulatedPool.+","",.))
		temp$Additional_Mt_Percent <- ifelse(grepl("mt", temp$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_mt_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
		temp$Additional_Ambient_Percent <- ifelse(grepl("ambient", temp$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_ambient_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
		temp$Downsampled <- ifelse(grepl("subsampled", temp$pool), TRUE,FALSE)

		# print(temp)
		correct <- temp[status == "Correct"]
		incorrect1 <- temp[status == "doublet > singlet"]
		incorrect1$status2 <- "doublet > singlet"
		incorrect2 <- temp[status == "singlet > doublet"]
		incorrect2$status2 <- "singlet > doublet"


		## Pull random number correct to match number of incorrect
		if (nrow(correct) < (nrow(incorrect1) + nrow(incorrect2))){
			if (nrow(incorrect1) < nrow(incorrect2)){
				incorrect2 <- incorrect2[sample(1:nrow(incorrect2), nrow(correct) - nrow(incorrect1), replace = FALSE),]
			} else if (nrow(incorrect1) > nrow(incorrect2)){
				incorrect1 <- incorrect1[sample(1:nrow(incorrect1), nrow(correct) - nrow(incorrect2), replace = FALSE),]
			}
		}
		correct_updated1 <- correct[sample(1:nrow(correct), nrow(incorrect1) + nrow(incorrect2), replace = FALSE),]
		correct_updated1$status2 <- c(rep("doublet > singlet", nrow(incorrect1)), rep("singlet > doublet", nrow(incorrect2)))

		if(nrow(incorrect1) > 0 & nrow(incorrect2) > 0){
			updated <- rbind(correct_updated1, rbind(incorrect1, incorrect2))
		} else if (nrow(incorrect1) > 0){
			updated <- rbind(correct_updated1, incorrect1)
		} else if (nrow(incorrect2) > 0){
			updated <- rbind(correct_updated1, incorrect2)
		} else {
			updated <- data.table(temp[0,], "status2" = NA)[0,]
			# updated <- NULL
		}

		# print(updated)
		if(!is.null(updated)){
			return(updated)
		} else {
			return(NULL)
		}
	})
	# print(temp2)
	temp3 <- do.call(rbind, temp2)
	temp3$status <- as.character(temp3$status)
	temp3$status2 <- as.character(temp3$status2)

	temp4 <- as.data.frame(dcast(temp3, pool + Number_Individuals + Barcode ~ software, value.var = c("status", "status2")))
	rownames(temp4) <- temp4$Barcode

	return(temp4)
})

names(soft_correct_incorrcet_list) <- names(singlets_sim_barcode_list_annos_doublets_list)



### Read in data as seurat objects ###
seurat_list <- lapply(pools, function(x){
	print(x)
	Read10X(paste0(sim_dir,x,"/matrix_out/"), gene.column = 1)
})

seurat_list <- lapply(seurat_list, function(x){
	rownames(x) <- gsub("\\.\\d+", "", rownames(x))
	return(x)
})

seurat_obj_list <- lapply(seurat_list, function(x){
	CreateSeuratObject(counts = x)
})

rb_genes <- fread(paste0(demuxafy_dir,"/RibosomalGeneList_GeneID_ENSG.txt"))
mt_genes <- fread(paste0(demuxafy_dir,"/MtGeneList_GeneID_ENSG.txt"))


seurat_obj_list <- lapply(seurat_obj_list, function(x){
	print(x)
	rb_list <- rb_genes$ENSG[rb_genes$ENSG %in% rownames(x)]
	mt_list <- mt_genes$ENSG[mt_genes$ENSG %in% rownames(x)]
	x[["rb.percent"]] <- PercentageFeatureSet(x, features = rb_list)
	x[["mt.percent"]] <- PercentageFeatureSet(x, features = mt_list)
	return(x)
})

names(seurat_obj_list) <- pools


seurat_obj_list_meta <- lapply(names(seurat_obj_list), function(x){
	AddMetaData(seurat_obj_list[[x]], soft_correct_incorrcet_list[[x]])
})

metadata <- lapply(seurat_obj_list_meta, function(x){
	tmp <- data.table(x@meta.data)
	# print(tmp)
	if (!("status_DoubletDecon" %in% colnames(tmp))){
		tmp$status_DoubletDecon <- NA
		tmp$status2_DoubletDecon <- NA
	}
	if (!("status_scDblFinder" %in% colnames(tmp))){
		tmp$status_scDblFinder <- NA
		tmp$status2_scDblFinder <- NA
	}
	if (!("status_scDblFinder_known_doublets" %in% colnames(tmp))){
		tmp$status_scDblFinder_known_doublets <- NA
		tmp$status2_scDblFinder_known_doublets <- NA
	}
	if (!("status_solo" %in% colnames(tmp))){
		tmp$status_solo <- NA
		tmp$status2_solo <- NA
	}
	if (!("status_scSplit" %in% colnames(tmp))){
		tmp$status_scSplit <- NA
		tmp$status2_scSplit <- NA
	}
	if (!("status_souporcell" %in% colnames(tmp))){
		tmp$status_souporcell <- NA
		tmp$status2_souporcell <- NA
	}
	tmp <- tmp[,c("orig.ident", "nCount_RNA","nFeature_RNA", "rb.percent","mt.percent", "pool","Number_Individuals", "Barcode","status_demuxalot", "status_demuxalot_refined", "status_dropulation", "status_DoubletDecon", "status_DoubletDetection", "status_DoubletFinder", "status_demuxlet", "status_freemuxlet", "status_scDblFinder", "status_scDblFinder_known_doublets", "status_scSplit", "status_scds", "status_scrublet", "status_solo", "status_souporcell", "status_vireo", "status2_demuxalot", "status2_demuxalot_refined", "status2_dropulation", "status2_DoubletDecon", "status2_DoubletDetection", "status2_DoubletFinder", "status2_demuxlet", "status2_freemuxlet", "status2_scDblFinder", "status2_scDblFinder_known_doublets", "status2_scSplit", "status2_scds", "status2_scrublet", "status2_solo", "status2_souporcell", "status2_vireo")]
	return(tmp)
})

lapply(metadata, ncol)


metadata_df <- do.call(rbind, metadata)

metadata_df_long <- melt(metadata_df, measure.vars = c(paste0("status_",softwares), paste0("status2_",softwares)), variable.name = "software", value.name = c("status"))
metadata_df_long <- metadata_df_long[!is.na(metadata_df_long$status)]
metadata_df_long$type <- ifelse(grepl("status_", metadata_df_long$software), "status", "status2")
metadata_df_long$software <- gsub("status_", "", metadata_df_long$software) %>% gsub("status2_", "", .)

metadata_df_long_wide <- dcast(metadata_df_long, orig.ident + nCount_RNA + nFeature_RNA + rb.percent + mt.percent + pool + Number_Individuals + Barcode + software~ type, value.var = "status")
metadata_df_long_wide$Number_Individuals <- factor(metadata_df_long_wide$Number_Individuals, levels = c(2,4,8,16,32,64,128))

fwrite(metadata_df_long_wide, paste0(outdir,"metadata_df.tsv"), sep = "\t")

metadata_df_long_wide <- fread(paste0(outdir,"metadata_df.tsv"), sep = "\t")



### Make barplots + SE ###
metadata_df_long_wide_avg <- metadata_df_long_wide[,.(mt_mean = mean(mt.percent, na.rm = TRUE), rb_mean = mean(rb.percent, na.rm = TRUE), mean_UMI = mean(nCount_RNA), mean_features = mean(nFeature_RNA)),.(pool,software,status,status2,Number_Individuals)]

metadata_df_long_wide_avg <- softwares_dt[metadata_df_long_wide_avg, on = "software"]
metadata_df_long_wide_avg$Software <- gsub(" ", "\n", metadata_df_long_wide_avg$Software)


p_mt_false_sing_demux <- ggbarplot(metadata_df_long_wide_avg[software %in% demultiplexing_list & status2 == "doublet > singlet"], 
								x = "Number_Individuals", y = "mt_mean", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							rotate_x_text() +
							scale_fill_manual(values = c("#3dbeb0","#cd8d5a")) +
							ylab("Mt %") +
							theme_classic() +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(p_mt_false_sing_demux, filename = paste0(outdir, "bar_false_singlet_demultiplexing_mt.png"), height = 2, width = 9)
ggsave(p_mt_false_sing_demux, filename = paste0(outdir, "bar_false_singlet_demultiplexing_mt.pdf"), height = 2, width = 9)



p_mt_false_sing_doublet_detecting <- ggbarplot(metadata_df_long_wide_avg[software %in% doublet_detection_list & status2 == "doublet > singlet"], 
								x = "Number_Individuals", y = "mt_mean", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#cd8d5a")) +
							ylab("Mt %") +
							theme_classic()
ggsave(p_mt_false_sing_doublet_detecting, filename = paste0(outdir, "bar_false_singlet_doublet_detecting_mt.png"), height = 2, width = 9)
ggsave(p_mt_false_sing_doublet_detecting, filename = paste0(outdir, "bar_false_singlet_doublet_detecting_mt.pdf"), height = 2, width = 9)

p_mt_false_doub_demux <- ggbarplot(metadata_df_long_wide_avg[software %in% demultiplexing_list & status2 == "singlet > doublet"], 
								x = "Number_Individuals", y = "mt_mean", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#dea9cc")) +
							ylab("Mt %") +
							theme_classic()
ggsave(p_mt_false_doub_demux, filename = paste0(outdir, "bar_false_doublet_demultiplexing_mt.png"), height = 2, width = 9)
ggsave(p_mt_false_doub_demux, filename = paste0(outdir, "bar_false_doublet_demultiplexing_mt.pdf"), height = 2, width = 9)

p_mt_false_doub_doublet_detecting <- ggbarplot(metadata_df_long_wide_avg[software %in% doublet_detection_list & status2 == "singlet > doublet"], 
								x = "Number_Individuals", y = "mt_mean", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#dea9cc")) +
							ylab("Mt %") +
							theme_classic()
ggsave(p_mt_false_doub_doublet_detecting, filename = paste0(outdir, "bar_false_doublet_doublet_detecting_mt.png"), height = 2, width = 9)
ggsave(p_mt_false_doub_doublet_detecting, filename = paste0(outdir, "bar_false_doublet_doublet_detecting_mt.pdf"), height = 2, width = 9)


mt_stats <- data.table(unique(metadata_df_long_wide_avg[,c("software", "Number_Individuals","status2")]))
mt_stats$p <- 0
mt_stats$statistic <- 0

for (row in 1:nrow(mt_stats)){
	soft <- mt_stats[row, software]
	n <- mt_stats[row, Number_Individuals]
	type <- mt_stats[row, status2]
	if (nrow(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]) >= 3){
		# print(nrow(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]))
		test <- t.test(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == type]$mt_mean,
				metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]$mt_mean)
		mt_stats[status2 == type & Number_Individuals == n & software == soft]$p <- as.numeric(test$p.value)
		mt_stats[status2 == type & Number_Individuals == n & software == soft]$statistic <- test$statistic
	} else {
		# mt_stats <- mt_stats[!(status2 == type & Number_Individuals == n & demultiplexing_list == soft)]
		mt_stats[status2 == type & Number_Individuals == n & software == soft]$p <- NA
		mt_stats[status2 == type & Number_Individuals == n & software == soft]$statistic <- NA
	}
}

mt_stats$bonferroni_p <- p.adjust(mt_stats$p, method = "bonferroni")


fwrite(mt_stats, paste0(outdir,"mt_percent_stats.tsv"), sep = "\t")
mt_stats <- fread(paste0(outdir,"mt_percent_stats.tsv"), sep = "\t")


min(mt_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$p)
max(mt_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$p)
min(mt_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$statistic)
max(mt_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$statistic)

min(mt_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$p)
max(mt_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$p)
min(mt_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$statistic)
max(mt_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$statistic)




# ## Rb percent ##
p_rb_false_sing_demux <- ggbarplot(metadata_df_long_wide_avg[software %in% demultiplexing_list & status2 == "doublet > singlet"], 
								x = "Number_Individuals", y = "rb_mean", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#cd8d5a")) +
							ylab("Rb %") +
							theme_classic()
ggsave(p_rb_false_sing_demux, filename = paste0(outdir, "bar_false_singlet_demultiplexing_rb.png"), height = 2, width = 9)
ggsave(p_rb_false_sing_demux, filename = paste0(outdir, "bar_false_singlet_demultiplexing_rb.pdf"), height = 2, width = 9)

p_rb_false_sing_doublet_detecting <- ggbarplot(metadata_df_long_wide_avg[software %in% doublet_detection_list & status2 == "doublet > singlet"], 
								x = "Number_Individuals", y = "rb_mean", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#cd8d5a")) +
							ylab("Rb %") +
							theme_classic()
ggsave(p_rb_false_sing_doublet_detecting, filename = paste0(outdir, "bar_false_singlet_doublet_detecting_rb.png"), height = 2, width = 9)
ggsave(p_rb_false_sing_doublet_detecting, filename = paste0(outdir, "bar_false_singlet_doublet_detecting_rb.pdf"), height = 2, width = 9)

p_rb_false_doub_demux <- ggbarplot(metadata_df_long_wide_avg[software %in% demultiplexing_list & status2 == "singlet > doublet"], 
								x = "Number_Individuals", y = "rb_mean", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#dea9cc")) +
							ylab("Rb %") +
							theme_classic()
ggsave(p_rb_false_doub_demux, filename = paste0(outdir, "bar_false_doublet_demultiplexing_rb.png"), height = 2, width = 9)
ggsave(p_rb_false_doub_demux, filename = paste0(outdir, "bar_false_doublet_demultiplexing_rb.pdf"), height = 2, width = 9)

p_rb_false_doub_doublet_detecting <- ggbarplot(metadata_df_long_wide_avg[software %in% doublet_detection_list & status2 == "singlet > doublet"], 
								x = "Number_Individuals", y = "rb_mean", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#dea9cc")) +
							ylab("Rb %") +
							theme_classic()
ggsave(p_rb_false_doub_doublet_detecting, filename = paste0(outdir, "bar_false_doublet_doublet_detecting_rb.png"), height = 2, width = 9)
ggsave(p_rb_false_doub_doublet_detecting, filename = paste0(outdir, "bar_false_doublet_doublet_detecting_rb.pdf"), height = 2, width = 9)


rb_stats <- data.table(unique(metadata_df_long_wide_avg[,c("software", "Number_Individuals","status2")]))
rb_stats$p <- 0
rb_stats$statistic <- 0

for (row in 1:nrow(rb_stats)){
	print(row)
	print(rb_stats)
	soft <- rb_stats[row, software]
	n <- rb_stats[row, Number_Individuals]
	type <- rb_stats[row, status2]
	if (nrow(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]) >= 3){
		test <- t.test(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == type]$rb_mean,
				metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]$rb_mean)
		rb_stats[status2 == type & Number_Individuals == n & software == soft]$p <- as.numeric(test$p.value)
		rb_stats[status2 == type & Number_Individuals == n & software == soft]$statistic <- test$statistic
	} else {
		print("not enough")
		rb_stats[status2 == type & Number_Individuals == n & software == soft]$p <- NA
		rb_stats[status2 == type & Number_Individuals == n & software == soft]$statistic <- NA
		# rb_stats <- rb_stats[!(status2 == type & Number_Individuals == n & software == soft)]
	}
}

rb_stats <- rb_stats[!is.na(statistic)]

rb_stats$bonferroni_p <- p.adjust(rb_stats$p, method = "bonferroni")

fwrite(rb_stats, paste0(outdir,"rb_percent_stats.tsv"), sep = "\t")
rb_stats <- fread(paste0(outdir,"rb_percent_stats.tsv"), sep = "\t")


min(rb_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "singlet > doublet"]$p)
max(rb_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "singlet > doublet"]$p)
min(rb_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "singlet > doublet"]$statistic)
max(rb_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "singlet > doublet"]$statistic)


rb_stats[bonferroni_p < 0.05][order(status2, software, Number_Individuals)]



## UMI ##
p_umi_false_sing_demux <- ggbarplot(metadata_df_long_wide_avg[software %in% demultiplexing_list & status2 == "doublet > singlet"], 
								x = "Number_Individuals", y = "mean_UMI", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#cd8d5a")) +
							ylab("N UMI") +
							theme_classic()
ggsave(p_umi_false_sing_demux, filename = paste0(outdir, "bar_false_singlet_demultiplexing_umi.png"), height = 2, width = 9)
ggsave(p_umi_false_sing_demux, filename = paste0(outdir, "bar_false_singlet_demultiplexing_umi.pdf"), height = 2, width = 9)

p_umi_false_sing_doublet_detecting <- ggbarplot(metadata_df_long_wide_avg[software %in% doublet_detection_list & status2 == "doublet > singlet"], 
								x = "Number_Individuals", y = "mean_UMI", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#cd8d5a")) +
							ylab("N UMI") +
							theme_classic()
ggsave(p_umi_false_sing_doublet_detecting, filename = paste0(outdir, "bar_false_singlet_doublet_detecting_umi.png"), height = 2, width = 9)
ggsave(p_umi_false_sing_doublet_detecting, filename = paste0(outdir, "bar_false_singlet_doublet_detecting_umi.pdf"), height = 2, width = 9)

p_umi_false_doub_demux <- ggbarplot(metadata_df_long_wide_avg[software %in% demultiplexing_list & status2 == "singlet > doublet"], 
								x = "Number_Individuals", y = "mean_UMI", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#dea9cc")) +
							ylab("N UMI") +
							theme_classic()
ggsave(p_umi_false_doub_demux, filename = paste0(outdir, "bar_false_doublet_demultiplexing_umi.png"), height = 2, width = 9)
ggsave(p_umi_false_doub_demux, filename = paste0(outdir, "bar_false_doublet_demultiplexing_umi.pdf"), height = 2, width = 9)

p_umi_false_doub_doublet_detecting <- ggbarplot(metadata_df_long_wide_avg[software %in% doublet_detection_list & status2 == "singlet > doublet"], 
								x = "Number_Individuals", y = "mean_UMI", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#dea9cc")) +
							ylab("N UMI") +
							theme_classic()
ggsave(p_umi_false_doub_doublet_detecting, filename = paste0(outdir, "bar_false_doublet_doublet_detecting_umi.png"), height = 2, width = 9)
ggsave(p_umi_false_doub_doublet_detecting, filename = paste0(outdir, "bar_false_doublet_doublet_detecting_umi.pdf"), height = 2, width = 9)



umi_stats <- data.table(unique(metadata_df_long_wide_avg[,c("software", "Number_Individuals","status2")]))
umi_stats$p <- 0
umi_stats$statistic <- 0


for (row in 1:nrow(umi_stats)){
	soft <- umi_stats[row, software]
	n <- umi_stats[row, Number_Individuals]
	type <- umi_stats[row, status2]
	if (nrow(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]) > 3){
		test <- t.test(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == type]$mean_UMI,
				metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]$mean_UMI)
		umi_stats[status2 == type & Number_Individuals == n & software == soft]$p <- as.numeric(test$p.value)
		umi_stats[status2 == type & Number_Individuals == n & software == soft]$statistic <- test$statistic
	} else {
		umi_stats[status2 == type & Number_Individuals == n & software == soft]$p <- NA
		umi_stats[status2 == type & Number_Individuals == n & software == soft]$statistic <- NA
		# umi_stats <- umi_stats[!(status2 == type & Number_Individuals == n & software == soft)]
	}
}

umi_stats <- umi_stats[!is.na(statistic)]

umi_stats$bonferroni_p <- p.adjust(umi_stats$p, method = "bonferroni")

umi_stats[bonferroni_p < 0.05][order(status2, software, Number_Individuals)]
umi_stats[bonferroni_p < 0.05 & bonferroni_p > 0.01][order(status2, software, Number_Individuals)]
umi_stats[bonferroni_p < 0.01 & bonferroni_p > 0.001][order(status2, software, Number_Individuals)]
umi_stats[bonferroni_p > 0.05][order(status2, software, Number_Individuals)]

fwrite(umi_stats, paste0(outdir,"umi_stats.tsv"), sep = "\t")
umi_stats <- fread(paste0(outdir,"umi_stats.tsv"), sep = "\t")


min(umi_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "singlet > doublet"]$p)
max(umi_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "singlet > doublet"]$p)
min(umi_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "singlet > doublet"]$statistic)
max(umi_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "singlet > doublet"]$statistic)




# ## N Features ##
p_n_features_false_sing_demux <- ggbarplot(metadata_df_long_wide_avg[software %in% demultiplexing_list & status2 == "doublet > singlet"], 
								x = "Number_Individuals", y = "mean_features", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#cd8d5a")) +
							ylab("N UMI") +
							theme_classic()
ggsave(p_n_features_false_sing_demux, filename = paste0(outdir, "bar_false_singlet_demultiplexing_n_features.png"), height = 2, width = 9)
ggsave(p_n_features_false_sing_demux, filename = paste0(outdir, "bar_false_singlet_demultiplexing_n_features.pdf"), height = 2, width = 9)

p_n_features_false_sing_doublet_detecting <- ggbarplot(metadata_df_long_wide_avg[software %in% doublet_detection_list & status2 == "doublet > singlet"], 
								x = "Number_Individuals", y = "mean_features", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#cd8d5a")) +
							ylab("N UMI") +
							theme_classic()
ggsave(p_n_features_false_sing_doublet_detecting, filename = paste0(outdir, "bar_false_singlet_doublet_detecting_n_features.png"), height = 2, width = 9)
ggsave(p_n_features_false_sing_doublet_detecting, filename = paste0(outdir, "bar_false_singlet_doublet_detecting_n_features.pdf"), height = 2, width = 9)

p_n_features_false_doub_demux <- ggbarplot(metadata_df_long_wide_avg[software %in% demultiplexing_list & status2 == "singlet > doublet"], 
								x = "Number_Individuals", y = "mean_features", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#dea9cc")) +
							ylab("N UMI") +
							theme_classic()
ggsave(p_n_features_false_doub_demux, filename = paste0(outdir, "bar_false_doublet_demultiplexing_n_features.png"), height = 2, width = 9)
ggsave(p_n_features_false_doub_demux, filename = paste0(outdir, "bar_false_doublet_demultiplexing_n_features.pdf"), height = 2, width = 9)

p_n_features_false_doub_doublet_detecting <- ggbarplot(metadata_df_long_wide_avg[software %in% doublet_detection_list & status2 == "singlet > doublet"], 
								x = "Number_Individuals", y = "mean_features", fill = "status", 
								add = "mean_se", 
								facet.by = c("Software"), nrow = 1, 
								position = position_dodge(0.7),
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							scale_fill_manual(values = c("#3dbeb0","#dea9cc")) +
							ylab("N UMI") +
							theme_classic()
ggsave(p_n_features_false_doub_doublet_detecting, filename = paste0(outdir, "bar_false_doublet_doublet_detecting_n_features.png"), height = 2, width = 9)
ggsave(p_n_features_false_doub_doublet_detecting, filename = paste0(outdir, "bar_false_doublet_doublet_detecting_n_features.pdf"), height = 2, width = 9)


metadata_df_long_wide_avg$software <- factor(metadata_df_long_wide_avg$software, levels = c(demultiplexing_list, doublet_detection_list))


features_stats <- data.table(unique(metadata_df_long_wide_avg[,c("software", "Number_Individuals","status2")]))
features_stats$p <- 0
features_stats$statistic <- 0

for (row in 1:nrow(features_stats)){
	soft <- features_stats[row, software]
	n <- features_stats[row, Number_Individuals]
	type <- features_stats[row, status2]
	if (nrow(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]) >3){
		test <- t.test(metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == type]$mean_features,
				metadata_df_long_wide_avg[software == soft & status2 == type & Number_Individuals == n & status == "Correct"]$mean_features)
		features_stats[status2 == type & Number_Individuals == n & software == soft]$p <- as.numeric(test$p.value)
		features_stats[status2 == type & Number_Individuals == n & software == soft]$statistic <- test$statistic
	} else {
		features_stats[status2 == type & Number_Individuals == n & software == soft]$p <- NA
		features_stats[status2 == type & Number_Individuals == n & software == soft]$statistic <- NA
	}
}

features_stats <- features_stats[!is.na(statistic)]

features_stats$bonferroni_p <- p.adjust(features_stats$p, method = "bonferroni")

features_stats[bonferroni_p < 0.05][order(status2, software, Number_Individuals)]
features_stats[bonferroni_p < 0.05 & bonferroni_p > 0.01][order(status2, software, Number_Individuals)]
features_stats[bonferroni_p < 0.01 & bonferroni_p > 0.001][order(status2, software, Number_Individuals)]
features_stats[bonferroni_p > 0.05][order(status2, software, Number_Individuals)]

fwrite(features_stats, paste0(outdir,"features_stats.tsv"), sep = "\t")
features_stats <- fread(paste0(outdir,"features_stats.tsv"), sep = "\t")


min(features_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$p)
max(features_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$p)
min(features_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$statistic)
max(features_stats[software %in% demultiplexing_list & p < 0.05 & status2 == "doublet > singlet"]$statistic)


umi_features_cor <- cor.test(features_stats$statistic, umi_stats$statistic, method = "spearman")

pt(q = as.numeric(umi_features_cor$statistic), df = 215, lower.tail=FALSE)*2



##### Get the total percent correct per software in each pool #####
soft_annotation_correct_list <- lapply(names(singlets_sim_barcode_list_annos_doublets_list), function(x){
	temp2 <- lapply(colnames(singlets_list[[1]])[2:ncol(singlets_list[[1]])], function(y){
		temp <- data.table(matrix(nrow = nrow(singlets_sim_barcode_list_annos_doublets_list[[x]]), ncol = 3))
		colnames(temp) <- c("Barcode" ,"Reference", "Predicted")
		temp$Barcode <- singlets_sim_barcode_list_annos_doublets_list[[x]][,"Simulated_Barcode"]
		temp$Reference <- ifelse(singlets_sim_barcode_list_annos_doublets_list[[x]]$doublet_singlet == "singlet", singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1, paste0("doublet_",singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals))
		temp$Predicted <- ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 1 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "singlet"), singlets_sim_barcode_list_annos_doublets_list[[x]]$CellType1,
								ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 1 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "doublet"), "doublet > singlet",
									ifelse((singlets_sim_barcode_list_annos_doublets_list[[x]][,..y] == 0 & singlets_sim_barcode_list_annos_doublets_list[[x]][,"doublet_singlet"] == "doublet"), paste0("doublet_", singlets_sim_barcode_list_annos_doublets_list[[x]]$heterotypic_homotypic, "_", singlets_sim_barcode_list_annos_doublets_list[[x]]$individuals), "singlet > doublet")))
		temp$software <- y
		temp$status <- factor(ifelse(temp$Reference == temp$Predicted, "Correct", "Incorrect"), levels = c("Correct", "Incorrect"))
		temp$Reference_status <- ifelse(grepl("doublet", temp$Reference), temp$Reference, "singlet")


		temp3 <- data.table(prop.table(table(paste0(temp$software), temp$status), 1))
		colnames(temp3) <- c("Reference_software", "status", "N")
		temp3 <- temp3[, c("Reference", "software") := tstrsplit(Reference_software, "-", fixed = TRUE)]
		temp3$Reference_software <- NULL

		return(temp3)
	})
	temp3 <- do.call(rbind, temp2)
	temp3$pool <- x
	temp3$Number_Individuals <- as.numeric(gsub("size", "", temp3$pool) %>% gsub("_SimulatedPool.+","",.))
	temp3$Additional_Mt_Percent <- ifelse(grepl("mt", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_mt_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	temp3$Additional_Ambient_Percent <- ifelse(grepl("ambient", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_ambient_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	temp3$Downsampled <- ifelse(grepl("subsampled", temp3$pool), TRUE,FALSE)
	temp3$Uneven <- ifelse(grepl("uneven", temp3$pool), (100*as.numeric(gsub("size\\d+_SimulatedPool\\d+_mt_","",rownames(pct_correct_avg)) %>% gsub("pctl","", .))), 0)
	return(temp3)
})



soft_correct <- do.call(rbind, soft_annotation_correct_list)

soft_correct$software <- factor(soft_correct$software, levels = softwares)

soft_correct <- data.table(soft_correct)

soft_correct$software_type <- ifelse(soft_correct$software %in% demultiplexing_list, "Demultiplexing\nSoftwares", "Doublet Detecting\nSoftwares")
soft_correct$Number_Individuals <- factor(soft_correct$Number_Individuals)

max(soft_correct[status == "Incorrect"]$N, na.rm = TRUE)
min(soft_correct[status == "Incorrect"]$N, na.rm = TRUE)

soft_correct <- softwares_dt[soft_correct,on = "software"]


##### Make figure showing the % incorrect
pNormal_demultiplexing_doublet_cell_type <- ggplot(soft_correct[Additional_Mt_Percent == 0 & Additional_Ambient_Percent == 0 & Downsampled == FALSE & status == "Incorrect" & Uneven == 0], aes(Number_Individuals, N*100, color = Software)) +
	# geom_boxplot(outlier.size = 0.25, aes(fill = software), alpha = 0.5) +
	geom_jitter(size = 0.25, width = 0.25, height = 0) +
	# facet_wrap(vars(Reference), nrow = 5) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	theme_classic() +
	scale_color_manual(values = software_colors) +
	scale_fill_manual(values = software_colors) +
	ylab("Percent Incorrect") +
	xlab("Number Multiplexed Individuals") +
	geom_smooth(aes(group = software), se = FALSE, method = loess, size = 0.6) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
		legend.position = "top") +
	theme(legend.position = "none")

save_figs(pNormal_demultiplexing_doublet_cell_type,  paste0(outdir, "percent_incorrect"), width = 7, height = 10)

