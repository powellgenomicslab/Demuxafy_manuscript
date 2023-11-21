# Script information ------------------------------------------------------

# title: used to annotate each of the droplets identified as singlets by all the methods (doublet detecting and demultiplexing)
# run after completing Demuxafy_manuscript/scripts/PBMC/empirical/annotate_cells.R


# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library(data.table)

# Secondary
library("Seurat")
library("SeuratDisk")
library("RColorBrewer")
library("viridis") 


# Read data ---------------------------------------------------------------

# Input
dir = "/path/to/output/PBMC/simulation/" ### the directory that contains all the results from cellranger for the fibroblast pools
out = paste0(dir,"/Round1Overlap/CellClassification_noRefDoublets/")
meta_file = "/path/to/Demuxafy_manuscript/files/PBMC/PBMC_sample.meta.txt"
ref_path = paste0("/path/tp/Demuxafy_manuscript/files/PBMC/pbmc_multimodal.h5seurat") ## available on Zenodo

dir.create(out)
sample_list <- dir(path = dir, pattern = "OneK1K_scRNA_Sample")
sample_list <- sample_list[!grepl("V1", sample_list)]

softwares <- c("demuxalot", "demuxalot_refined", "demuxlet","dropulation", "freemuxlet","scSplit","souporcell","vireo", "DoubletDetection", "DoubletFinder", "DoubletDecon", "scds", "scDblFinder", "scrublet", "solo")
demultiplexing_softwares <- c("demuxalot", "demuxalot_refined", "demuxlet","dropulation", "freemuxlet","scSplit","souporcell","vireo")
doublet_detection_softwares <- c("DoubletDetection", "DoubletFinder", "DoubletDecon", "scds", "scDblFinder", "scrublet", "solo")



# Read file
inicio("Reading gene expression data")
##### Make Seurat Objects #####
dirs10x <- paste0("/path/to/10x/data/dir/", sample_list, "_V1/outs/filtered_gene_bc_matrices/Homo_sapiens_GRCh38p10/")
## Read in expression data
counts_list <- lapply(dirs10x, function(x){
    Read10X(x, gene.column = 2)
})
names(counts_list) <- sample_list

## Add poolnames to cell names so can easily match back if there are any cells with same barcoe
counts_list <- lapply(names(counts_list), function(x){
    colnames(counts_list[[x]]) <- paste0(x, "_", colnames(counts_list[[x]]))
    return(counts_list[[x]])
})
names(counts_list) <- sample_list



result_files <- paste0(dir,sample_list,"/CombinedResults/CombinedDropletAssignments.tsv")
names(result_files) <- sample_list
meta <- read_delim(meta_file, delim = "\t")

##### Read in Files #####
results_list <- lapply(result_files, function(x){
  fread(x, sep = "\t")
})
names(results_list) <- sample_list

##### Read in the key tables for freemuxlet, scSplit and souporcell #####
key <- readRDS(paste0(dir,"/Round1Overlap/PoolKeys.rds"))
key <- lapply(key, function(x){
  x$Cluster_ID <- gsub("CLUST","",x$Cluster_ID)
  x$Software <- paste0(x$Software, "_Assignment")
  return(x)
})
key <- key[sample_list]

##### Pivot the dataframes longer for joining#####
results_list_long <- lapply(results_list, function(x){
  tmp <- melt(x,measure.vars = paste0(demultiplexing_softwares, "_Assignment"),
								variable.name = "Software", value.name = "Assignment")
  return(tmp)
})

##### Left_join the common assignments to the dataframe #####
results_list_long <- lapply(names(key), function(x){
  print(x)
  temp <- data.table(key[[x]])[results_list_long[[x]], on = c("Cluster_ID" = "Assignment", "Software")]
  colnames(temp) <- gsub("Cluster_ID", "Assignment", colnames(temp))
  temp$Genotype_ID <- ifelse(temp$Software == "demuxalot_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Software == "demuxalot_refined_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Software == "demuxlet_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Software == "dropulation_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Software == "vireo_Assignment", temp$Assignment, temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(temp$Assignment == "doublet", "doublet", temp$Genotype_ID)
  temp$Genotype_ID <- ifelse(is.na(temp$Genotype_ID), "unassigned", temp$Genotype_ID)
  return(temp)
})
names(results_list_long) <- names(key)

##### Pivot wider to get the counts of singlet, doublets... and intersection across softwares #####
results_list_wide <- lapply(results_list_long, function(x){
  x$Assignment <- NULL
  x$Correlation <- NULL
  # xlong <- dcast(x, colnames(x)[!(colnames(x) %in% c("Barcode"))] ~ Software, value.var = "Genotype_ID")
  xlong <- data.table(pivot_wider(x,names_from = "Software", values_from = "Genotype_ID"))
  return(xlong)
})
names(results_list_wide) <- names(key)


##### Add information on number of singlet and doublet calls #####
results_list_wide <- lapply(results_list_wide, function(x){
    drop_type <- grep("_DropletType", colnames(x))
    x[,"Singlets"] <- rowSums(x[,..drop_type] == "singlet", na.rm = TRUE)
    x[,"Doublets"] <- rowSums(x[,..drop_type] == "doublet", na.rm = TRUE)
    demux_type <- colnames(x)[colnames(x) %in% (paste0(demultiplexing_softwares,"_DropletType"))]
    doub_type <- colnames(x)[colnames(x) %in% (paste0(doublet_detection_softwares,"_DropletType"))]
    x[,"Demultiplexing_Singlets"]<- rowSums(x[, ..demux_type] == "singlet", na.rm = TRUE)
    x[,"Demultiplexing_Doublets"]<- rowSums(x[, ..demux_type] == "doublet", na.rm = TRUE)  
    x[,"DoubletDetecting_Singlets"]<- rowSums(x[, ..doub_type] == "singlet", na.rm = TRUE)
    x[,"DoubletDetecting_Doublets"]<- rowSums(x[, ..doub_type] == "doublet", na.rm = TRUE)
    return(x)
})


results_list_wide <- lapply(names(results_list_wide), function(x){
  results_list_wide[[x]] <- as.data.frame(results_list_wide[[x]])
  rownames(results_list_wide[[x]]) <- paste0(x, "_", results_list_wide[[x]]$Barcode)
  return(results_list_wide[[x]])
})
names(results_list_wide) <- sample_list

seurat_list <- lapply(names(counts_list), function(x){
    CreateSeuratObject(counts = counts_list[[x]], meta.data = results_list_wide[[x]])
})
names(seurat_list) <- sample_list


fin()


pryr::mem_used()
seurat_list <- lapply(seurat_list, function(x){
  DefaultAssay(x) <- "RNA"
  return(x)
})
gc()
pryr::mem_used()



# Read reference data -----------------------------------------------------

reference <- LoadH5Seurat(ref_path) ## available on Zenodo
ref_umap <- DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave(ref_umap, filename = paste0(out,"reference_data_umap.png"))

graphs_wknn_temp <- reference@graphs$wknn[WhichCells(reference, idents = "Doublet", invert = TRUE),WhichCells(reference, idents = "Doublet", invert = TRUE)]
graphs_wsnn_temp <- reference@graphs$wsnn[WhichCells(reference, idents = "Doublet", invert = TRUE),WhichCells(reference, idents = "Doublet", invert = TRUE)]

reference_updated <- subset(reference, subset = celltype.l2 != "Doublet")
reference_updated@graphs$wknn <- graphs_wknn_temp
reference_updated@graphs$wsnn <- graphs_wsnn_temp

ref_umap_updated <- DimPlot(object = reference_updated, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
ggsave(ref_umap_updated, filename = paste0(out,"reference_data_noDoublets_umap.png"))




map_cells <- function(x, y, ref){
  
  cat("Annotating: ", y, "\n", sep = "")
  
  x <- SCTransform(x)
  
  anchors <- FindTransferAnchors(
    reference = ref,
    query = x,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50,
    recompute.residuals = FALSE
  )
  
  x <- MapQuery(
    anchorset = anchors,
    query = x,
    reference = ref,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  
  x
  
}


inicio("Annotate cells")
anno_data <- imap(seurat_list, map_cells, reference)
names(anno_data) <- sample_list

anno_data_updated <- imap(seurat_list, map_cells, reference_updated)
names(anno_data) <- sample_list

fin()



# Save metadata for each pool ---------------------------------------------

inicio("Save metadata")
iwalk(anno_data, ~ {
  saveRDS(.x[[]], here(out, .y %p% "_metadata.RDS"))
})
fin()


inicio("Save no_doublet metadata")
iwalk(anno_data_updated, ~ {
  saveRDS(.x[[]], here(out, .y %p% "_metadata_no_doublet.RDS"))
})
fin()

pryr::mem_used()
rm(seurat_list)


# Save pools --------------------------------------------------------------

inicio("Save objects")
iwalk(anno_data, ~ {
  cat("Saving: ", .y, "\n", sep = "")
  saveRDS(.x, here(out, .y %p% ".RDS"))
})
fin()

anno_data <- lapply(sample_list, function(x){
  readRDS(paste0(out,x,".RDS"))
})
names(anno_data) <- sample_list

inicio("Save no_doublet objects")
iwalk(anno_data_updated, ~ {
  cat("Saving: ", .y, "\n", sep = "")
  saveRDS(.x, here(out, .y %p% "_no_doublets.RDS"))
})
fin()

anno_data_updated <- lapply(sample_list, function(x){
  readRDS(paste0(out,x,"_no_doublets.RDS"))
})
names(anno_data_updated) <- sample_list



# Extract embeddings ------------------------------------------------------

i <- seq(2, length(anno_data))

spca <- map(anno_data, ~ .[["ref.spca"]])
spca_upadted <- map(anno_data_updated, ~ .[["ref.spca"]])
umap <- map(anno_data, ~ .[["ref.umap"]])
umap_updated <- map(anno_data_updated, ~ .[["ref.umap"]])
meta_data <- map(anno_data, ~ .[[]] %>% rownames_to_column("barcode"))
meta_data_updated <- map(anno_data_updated, ~ .[[]] %>% rownames_to_column("barcode"))



inicio("Combine spcas")
spca <- merge(spca[[1]], spca[i])
fin()


inicio("Combine umaps")
umap <- merge(umap[[1]], umap[i])
fin()


inicio("Making all souporcell_LogProbSinglet columns double")
meta_data <- lapply(meta_data, function(x){
  x$souporcell_LogProbSinglet <- as.numeric(as.character(x$souporcell_LogProbSinglet))
  return(x)
})

meta_data_updated <- lapply(meta_data_updated, function(x){
  x$souporcell_LogProbSinglet <- as.numeric(as.character(x$souporcell_LogProbSinglet))
  return(x)
})

fin()

inicio("Combine metadata")

meta_data <- bind_rows(meta_data)
meta_data$pool <- gsub("_[ATCG]+-1", "", rownames(meta_data))
meta_data_updated <- bind_rows(meta_data_updated)
meta_data_updated$pool <- gsub("_[ATCG]+-1", "", rownames(meta_data_updated))

fin()



# Plot umap ---------------------------------------------------------------

em <- Embeddings(umap)
em_meta <- cbind(meta_data, em)
em_meta_update <- left_join(em_meta, meta_data_updated[,c("Barcode", "pool", "predicted.celltype.l2", grep("dropulation", colnames(meta_data_updated), value = TRUE), grep("demuxalot", colnames(meta_data_updated), value = TRUE))], by = c("Barcode","pool"))
colnames(em_meta_update) <- gsub("predicted.celltype.l2.y", "predicted.celltype.l2_no_doublets", colnames(em_meta_update))
colnames(em_meta_update) <- gsub("predicted.celltype.l2.x", "predicted.celltype.l2", colnames(em_meta_update))



write_delim(em_meta_update, paste0(out, "metadata_cell_classification.tsv"), delim = "\t")
em_meta_update <- read_delim(paste0(out, "metadata_cell_classification.tsv"), delim = "\t")

colnames(em_meta_update)

### Subset just singlets for all softwares ###
em_meta_sing <- em_meta_update[which(em_meta_update$demuxalot_DropletType == "singlet" &
                em_meta_update$demuxalot_refined_DropletType == "singlet" &
                em_meta_update$demuxlet_DropletType == "singlet" &
                em_meta_update$dropulation_DropletType == "singlet" &
								em_meta_update$freemuxlet_DropletType == "singlet" &
								em_meta_update$scSplit_DropletType == "singlet" &
								em_meta_update$souporcell_DropletType == "singlet" &
								em_meta_update$vireo_DropletType == "singlet" &
								em_meta_update$scrublet_DropletType == "singlet" &
								em_meta_update$scds_DropletType == "singlet" &
								em_meta_update$DoubletDetection_DropletType == "singlet" &
								em_meta_update$DoubletFinder_DropletType == "singlet" &
								em_meta_update$solo_DropletType == "singlet" &
								em_meta_update$scDblFinder_DropletType == "singlet" &
								em_meta_update$DoubletDecon_DropletType == "singlet"),]



n <- em_meta_sing$predicted.celltype.l2 %>% unique() %>% length()
pal <- colorRampPalette(brewer.pal(9, "Set1"))(n)


p <- ggplot(em_meta_sing, aes(refUMAP_1, refUMAP_2, color = predicted.celltype.l2)) +
  geom_point(shape = 16, size = 0.01) +
  xlab("ref UMAP 1") +
  ylab("ref UMAP 2") +
  scale_color_manual(values = pal) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(here(out, "ref_umap.png"), p, width = 9, height = 4.4)



# Save data ---------------------------------------------------------------

saveRDS(umap, here(out, "umap.RDS"))
saveRDS(spca, here(out, "spca.RDS"))
saveRDS(em_meta, here(out, "metadata.RDS"))


##### Make figure that contains doublets in black for supplementary comparison #####
em_meta <- readRDS(here(out, "metadata.RDS"))
umap <- readRDS(here(out, "umap.RDS"))
em_meta <- cbind(em_meta, Embeddings(umap))



em_meta_sing$predicted.celltype.l2 <- factor(em_meta_sing$predicted.celltype.l2, levels = c(sort(unique(em_meta_sing$predicted.celltype.l2)[(!(unique(em_meta_sing$predicted.celltype.l2) %in% "Doublet"))]), "Doublet"))

n <- em_meta$predicted.celltype.l2 %>% unique() %>% length()
pal <- c(colorRampPalette(brewer.pal(9, "Set1"))(n-1), "black")


pDoub <- ggplot(em_meta_sing, aes(refUMAP_1, refUMAP_2, color = predicted.celltype.l2)) +
  geom_point(shape = 16, size = 0.01) +
  xlab("ref UMAP 1") +
  ylab("ref UMAP 2") +
  scale_color_manual(values = pal) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 5)))



ggsave(here(out, "/ref_umap_doub_black.png"), pDoub, width = 9, height = 4.4)



# Session info ------------------------------------------------------------

print_session(here(out))

