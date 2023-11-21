library(data.table)
library(tidyverse)
library("RColorBrewer")    
library(ggpubr)



##### Set up variables #####
dir <- "/path/to/output/PBMC/simulation/SimulatedPools/SimulatedOverlap/"
outdir <- paste0(dir,"chord/")

dir.create(outdir)

pools <- gsub("/chord_metrics.tsv", "", (list.files(dir, pattern = "chord_metrics.tsv", recursive = TRUE)))



##### Read in data #####
### Doublet Detecting Restuls ###
dd_list <- lapply(pools, function(pool){
    fread(paste0(dir,pool, "/doublet_detecting_metrics.tsv"))
})

dd_dt <- do.call(rbind, dd_list)


### Chord results ###
chord_list <- lapply(pools, function(pool){
    fread(paste0(dir,pool, "/chord_metrics.tsv"))
})

chord_dt <- do.call(rbind, chord_list)


colnames(chord_dt) <- paste0(colnames(chord_dt), "_chord")



##### Combine together in one dataframe #####
combined_dt <- chord_dt[dd_dt, on = c("pool_chord" = "pool")]



##### Calculate different for all metrics between combinations  and chord #####
## positive => Demuxafy is better
## negative => chord is better
combined_dt$true_singlet_detla <- combined_dt$true_singlet - combined_dt$true_singlet_chord
combined_dt$false_singlet_detla <- combined_dt$false_singlet - combined_dt$false_singlet_chord
combined_dt$true_doublet_detla <- combined_dt$true_doublet - combined_dt$true_doublet_chord
combined_dt$false_doublet_detla <- combined_dt$false_doublet - combined_dt$false_doublet_chord
combined_dt$accuracy_detla <- combined_dt$accuracy - combined_dt$accuracy_chord
combined_dt$balanced_accuracy_detla <- combined_dt$balanced_accuracy - combined_dt$balanced_accuracy_chord 
combined_dt$sensitivity_detla <- combined_dt$sensitivity - combined_dt$sensitivity_chord 
combined_dt$specificity_detla <- combined_dt$specificity - combined_dt$specificity_chord       
combined_dt$npv_detla <- combined_dt$npv - combined_dt$npv_chord       
combined_dt$ppv_detla <- combined_dt$ppv - combined_dt$ppv_chord 
combined_dt$precision_detla <- combined_dt$precision - combined_dt$precision_chord
combined_dt$mcc_detla <- combined_dt$mcc - combined_dt$mcc_chord



##### check if chord is better than the best combination #####
head(combined_dt[!is.na(mcc)][rev(order(mcc))], n = 100)
head(combined_dt[!is.na(balanced_accuracy)][rev(order(balanced_accuracy))], n = 100)


##### Check for the combination used by chord
combined_dt[grepl("^DoubletFinder-scds$", software)]


##### Make Figures #####
### make dataframe for plotting ###
chord_dt_same_cols <- chord_dt
colnames(chord_dt_same_cols) <- gsub("_chord","",colnames(chord_dt))
chord_dt_same_cols$method <- "chord"
chord_dt_same_cols$paper <- "Chord"
dd_dt$paper <- "Demuxafy"
col_order <- colnames(dd_dt)
chord_dt_same_cols <- chord_dt_same_cols[,..col_order]

combined_dt_long <- rbind(dd_dt, chord_dt_same_cols)

### Compare same method inclusion

pDFscds_mcc <- ggplot(combined_dt_long[(grepl("^DoubletFinder-scds$", software) & method == "majority_singlet") | paper == "Chord"], aes(x = factor(paper, levels = c("Demuxafy", "Chord")), y = mcc)) +
    geom_line(aes(group=pool), color = "black") +
    geom_point()+ 
    theme_classic() +
    facet_wrap(vars(size), nrow = 1) +
    xlab(NULL) +
    ylab("MCC") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(pDFscds_mcc, filename =paste0(outdir,"dobuletfinder_scds_comparison_mcc.png"), height = 2)
ggsave(pDFscds_mcc, filename =paste0(outdir,"dobuletfinder_scds_comparison_mcc.pdf"), height = 2)


paired_ttest_results <- data.table(size = unique(combined_dt$size),p = as.numeric(NA))

for (size_pool in unique(combined_dt$size)){
    if (nrow(combined_dt[size == size_pool &(grepl("^DoubletFinder-scds$", software) & method == "majority_singlet")]) >= 3) {
        paired_ttest_results[size == size_pool]$p <- t.test(combined_dt[size == size_pool &(grepl("^DoubletFinder-scds$", software) & method == "majority_singlet")]$mcc, combined_dt[size == size_pool &(grepl("^DoubletFinder-scds$", software) & method == "majority_singlet")]$mcc_chord, paired = TRUE, alternative = "two.sided")$p.value
    }
}



combined_dt_long_DFscds_long <- melt(combined_dt_long[(grepl("^DoubletFinder-scds$", software) & method == "majority_singlet") | paper == "Chord"], id.vars =c("pool", "size", "paper"),
                measure.vars =c("true_singlet", "false_singlet", "true_doublet", "false_doublet"),
								variable.name = "Metric", value.name = "Count")



pDFscds_prop_celltypes <- ggplot(combined_dt_long_DFscds_long, aes(fill = Metric, x = factor(paper, levels = c("Demuxafy", "Chord")), y = Count)) +
    geom_bar(position="fill", stat="identity") +
    theme_classic() +
    facet_wrap(vars(size), nrow = 1) +
    xlab(NULL) +
    ylab("Proportion") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_brewer(palette = "Paired")

ggsave(pDFscds_prop_celltypes, filename =paste0(outdir,"dobuletfinder_scds_comparison_bar.png"), height = 2)
ggsave(pDFscds_prop_celltypes, filename =paste0(outdir,"dobuletfinder_scds_comparison_bar.pdf"), height = 2)


pDFscds_prop_celltypes_nolegend <- ggplot(combined_dt_long_DFscds_long, aes(fill = Metric, x = factor(paper, levels = c("Demuxafy", "Chord")), y = Count)) +
    geom_bar(position="fill", stat="identity") +
    theme_classic() +
    facet_wrap(vars(size), nrow = 1) +
    xlab(NULL) +
    ylab("Proportion") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_brewer(palette = "Paired") +
    theme(legend.position = "none")

ggsave(pDFscds_prop_celltypes_nolegend, filename =paste0(outdir,"dobuletfinder_scds_comparison_bar_nolegend.png"), height = 2)
ggsave(pDFscds_prop_celltypes_nolegend, filename =paste0(outdir,"dobuletfinder_scds_comparison_bar_nolegend.pdf"), height = 2)


combined_dt_long_DFscds_summary <- combined_dt_long[,c("pool","paper", "size")]


combined_dt_long_DFscds_summary$true_singlet_proportion <- combined_dt_long$true_singlet/(rowSums(combined_dt_long[,c("true_singlet", "false_singlet", "true_doublet", "false_doublet")]))
combined_dt_long_DFscds_summary$false_singlet_proportion <- combined_dt_long$false_singlet/(rowSums(combined_dt_long[,c("true_singlet", "false_singlet", "true_doublet", "false_doublet")]))
combined_dt_long_DFscds_summary$true_doublet_proportion <- combined_dt_long$true_doublet/(rowSums(combined_dt_long[,c("true_singlet", "false_singlet", "true_doublet", "false_doublet")]))
combined_dt_long_DFscds_summary$false_doublet_proportion <- combined_dt_long$false_doublet/(rowSums(combined_dt_long[,c("true_singlet", "false_singlet", "true_doublet", "false_doublet")]))

combined_dt_long_DFscds_summary$paper <- factor(combined_dt_long_DFscds_summary$paper, levels = c("Demuxafy", "Chord"))


combined_dt_long_DFscds_summary_long <- melt(combined_dt_long_DFscds_summary, id.vars =c("pool", "size", "paper"),
                measure.vars =c("true_singlet_proportion", "false_singlet_proportion", "true_doublet_proportion", "false_doublet_proportion"),
								variable.name = "Metric", value.name = "Count")


pDFscds_prop_celltypes_error <- ggbarplot(combined_dt_long_DFscds_summary_long, 
								x = "paper", y = "Count", fill =  "Metric", 
								facet.by = c("size"), nrow = 1, 
								position = position_stack(),
								add = "mean_se", 
								size = 0.25) +
							theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
							ylab("Proportion") +
							theme_classic() +
                            scale_fill_brewer(palette = "Paired") +
                            theme(legend.position = "none")


ggsave(pDFscds_prop_celltypes_error, filename =paste0(outdir,"dobuletfinder_scds_comparison_bar_error.png"), height = 2)
ggsave(pDFscds_prop_celltypes_error, filename =paste0(outdir,"dobuletfinder_scds_comparison_bar_error.pdf"), height = 2)


