##### Run locally with qrsh with conda activate baseR402 #####

##### Load libraries #####
library(tidyverse)
library(ggplot2)
library(elpatron)
library(data.table)

##### Set up directories #####
dir <- "/path/to/output/PBMC/benchmarks/"
data_dir <- paste0(dir,"/qstat_usage/")
outdir <- paste0(dir,"/qstat_usage_figures/"

dir.create(outdir)

files <- list.files(data_dir)

unique_rules <- gsub("usage_OneK1K_scRNA_Sample\\d+_", "", files) %>% gsub("_qstat.txt","",.) %>% unique()
unique_rules2 <- gsub("usage_OneK1K_scRNA_Sample\\d+.", "", files2) %>% gsub("_qstat.txt","",.) %>% gsub("demuxalot_script.txt", "demuxalot_script", .) %>% unique()


rules_files_list <- lapply(unique_rules, function(x){
    files[grep(paste0("\\d_",x,"_qstat.txt"), files)]
})
names(rules_files_list) <- unique_rules

rules_files_list2 <- lapply(unique_rules2, function(x){
    if (x == "demuxalot_script"){
        files2[grep(paste0(x,".txt"), files2)]
    } else {
        files2[grep(paste0(x,"_qstat.txt"), files2)]
    }
})
names(rules_files_list2) <- unique_rules2


### Read in file ###
results <- lapply(rules_files_list, function(x){
    lapply(x, function(y){
        read_delim(paste0(data_dir,y), delim = ",", col_names = c("CPU_time", "MEM", "io", "vmem", "Total_mem"))
    })
})


### Fix names ###
results <- lapply(names(results), function(x){
    names(results[[x]]) <- gsub(paste0("_",x), "", rules_files_list[[x]]) %>% gsub("usage_", "", .) %>% gsub("_qstat.txt","",.)
    return(results[[x]])
})
names(results) <- unique_rules


### Combine list to dataframe ###
results <- lapply(results, function(x){
    lapply(names(x), function(y){
        x[[y]]$Pool <- y
        return(x[[y]])
    })
})

results_df <- lapply(results, function(x){
    do.call(rbind, x)
})



### Get rid of extra text in entries ###
results_df <- lapply(results_df, function(x){
    x$CPU_time <- gsub("cpu=","",x$CPU_time)
    x$MEM <- gsub("mem=","",x$MEM)
    x$io <- gsub("io=","",x$io)
    x$vmem <- gsub("vmem=","",x$vmem)
    x$Total_mem <- gsub("maxvmem=","",x$Total_mem)
    return(x)
})



### Convert Total_mem to same scale and time to seconds
results_df <- lapply(results_df, function(x){
    x$Total_mem_G <- ifelse(grepl("M", x$Total_mem), as.numeric(gsub("M","",x$Total_mem))/1024, 
                        ifelse(grepl("G", x$Total_mem), as.numeric(gsub("G","",x$Total_mem)),
                            ifelse(is.na(x$Total_mem), 0, 
                                ifelse(grepl("N/A", x$Total_mem),0,"missing"))))
    x$CPU_time_h <- ifelse(grepl("N/A", x$CPU_time), 0, as.numeric(hms_to_sec(x$CPU_time)/60/60))
    return(x)
})


## check for nas
lapply(results_df, function(x) any(is.na(x$Total_mem_G)))
lapply(results_df, function(x) any(is.na(x$CPU_time_h)))

lapply(results2_df, function(x) any(is.na(x$Total_mem_G)))
lapply(results2_df, function(x) any(is.na(x$CPU_time_h)))

### Combine jobs for same software ###
demuxlet <- c("popscle_pileup", "popscle_demuxlet")
freemuxlet <- c("popscle_pileup", "freemuxlet_pool_vcf", "popscle_freemuxlet")
scSplit <- c("scSplit_allele_matrices", "scSplit_bgzip", "scSplit_demultiplex", "scSplit_freebayes", "scSplit_genotypes", "scSplit_pool_vcf", "scSplit_regions", "scSplit_rmdupe", "scSplit_sam_body", "scSplit_sam_combine", "scSplit_sam_header", "scSplit_sort", "scSplit_subset_vcf", "scSplit_vcf_qual_filt")
souporcell <- c("souporcell", "souporcell_pool_vcf")
vireo <- c("cellSNP", "subset_vcf", "vireo")
DoubletDecon <- c("rhop0.6_DoubletDecon", "rhop0.7_DoubletDecon", "rhop0.8_DoubletDecon", "rhop0.9_DoubletDecon", "rhop1_DoubletDecon", "rhop1.1_DoubletDecon")
DoubletDetection <- c("DoubletDetection")
DoubletFinder <- c("DoubletFinder")
scdblfinder <- c("scdblfinder")
scds <- c("scds")
scrublet <- c("80_scrublet", "85_scrublet", "90_scrublet", "95_scrublet")
solo <- c("solo")
demuxalot <- c("demuxalot_script")
dropulation <- c("dropulation_assign", "dropulation_doublet", "dropulation_tag")

software_list <- list("Demuxlet" = demuxlet, "Freemuxlet" = freemuxlet, "ScSplit" = scSplit, "Souporcell" = souporcell, "Vireo" = vireo, "DoubletDecon" = DoubletDecon, "DoubletDetection" = DoubletDetection, "DoubletFinder" = DoubletFinder, "ScDblFinder" = scdblfinder, "Scds" = scds, "Scrublet" = scrublet, "Solo" = solo)
software_list2 <- list(Demuxalot = demuxalot, Dropulation = dropulation)


results_df_grouped_list <- list()
results_df_grouped_list <- lapply(names(software_list), function(x){
    results_df_grouped_list[[x]] <- list()
    lapply(software_list[[x]], function(y){
        results_df_grouped_list[[x]][[y]] <- results_df[[y]]
        results_df_grouped_list[[x]][[y]]$rule <- y
        results_df_grouped_list[[x]][[y]]$software_group <- x
        return(results_df_grouped_list[[x]][[y]])
    })
})


results2_df_grouped_list <- list()
results2_df_grouped_list <- lapply(names(software_list2), function(x){
    results2_df_grouped_list[[x]] <- list()
    lapply(software_list2[[x]], function(y){
        results2_df_grouped_list[[x]][[y]] <- results2_df[[y]]
        results2_df_grouped_list[[x]][[y]]$rule <- y
        results2_df_grouped_list[[x]][[y]]$software_group <- x
        return(results2_df_grouped_list[[x]][[y]])
    })
})


results_df_grouped_list_merged <- lapply(results_df_grouped_list, function(x){
    do.call(rbind, x)
})
names(results_df_grouped_list_merged) <- names(software_list)


results_df_grouped <- do.call(rbind, results_df_grouped_list_merged)


##### Make figures #####
pAllRules_CPUtime <- ggplot(results_df_grouped, aes(rule, CPU_time_h)) +
    geom_boxplot() +
    theme_classic() +
    facet_wrap(~ software_group,nrow = 1, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(pAllRules_CPUtime, filename = paste0(outdir, "AllRules_CPUtime.eps"), width = 15)


pAllRules_Memtime <- ggplot(results_df_grouped, aes(rule, as.numeric(Total_mem_G))) +
    geom_boxplot() +
    theme_classic() +
    facet_wrap(~ software_group,nrow = 1, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(pAllRules_Memtime, filename = paste0(outdir, "AllRules_mem.eps"), width = 15)



##### Combine Measures from same pool, same software supergroup
results_df_grouped_wide_mem_list <- lapply(names(results_df_grouped_list_merged), function(x){
    if (x == "DoubletDecon" | x == "scrublet"){
        temp <- results_df_grouped_list_merged[[x]][,c("Pool", "software_group", "Total_mem_G", "CPU_time_h")]
        temp$Total_mem_G <- as.numeric(as.character(temp$Total_mem_G))
        temp$CPU_time_h <- as.numeric(as.character(temp$CPU_time_h))
    } else {
        temp1 <- pivot_wider(results_df_grouped_list_merged[[x]][,c("Pool", "Total_mem_G", "rule","software_group")], names_from = rule, values_from = Total_mem_G )
        temp1$Total_mem_G <- rowSums(mutate_all(temp1[,3:ncol(temp1)], function(x) as.numeric(as.character(x))))
        temp1 <- temp1[,c("Pool", "software_group", "Total_mem_G")]
        temp2 <- pivot_wider(results_df_grouped_list_merged[[x]][,c("Pool", "CPU_time_h", "rule","software_group")], names_from = rule, values_from = CPU_time_h )
        temp2$CPU_time_h <- rowSums(mutate_all(temp2[,3:ncol(temp2)], function(x) as.numeric(as.character(x))))
        temp2 <- temp2[,c("Pool", "software_group", "CPU_time_h")] 
        temp <- inner_join(temp1, temp2)
    }
    return(temp)
})
names(results_df_grouped_wide_mem_list) <- names(results_df_grouped_list_merged)



##### Combine Measures from same pool, same software supergroup
results2_df_grouped_wide_mem_list <- lapply(names(results2_df_grouped_list_merged), function(x){
    temp1 <- pivot_wider(results2_df_grouped_list_merged[[x]][,c("Pool", "Total_mem_G", "rule","software_group")], names_from = rule, values_from = Total_mem_G )
    temp1$Total_mem_G <- rowSums(mutate_all(temp1[,3:ncol(temp1)], function(x) as.numeric(as.character(x))))
    temp1 <- temp1[,c("Pool", "software_group", "Total_mem_G")]
    temp2 <- pivot_wider(results2_df_grouped_list_merged[[x]][,c("Pool", "CPU_time_h", "rule","software_group")], names_from = rule, values_from = CPU_time_h )
    temp2$CPU_time_h <- rowSums(mutate_all(temp2[,3:ncol(temp2)], function(x) as.numeric(as.character(x))))
    temp2 <- temp2[,c("Pool", "software_group", "CPU_time_h")] 
    temp <- inner_join(temp1, temp2)
    return(temp)
})
names(results2_df_grouped_wide_mem_list) <- names(results2_df_grouped_list_merged)



results_df_grouped_wide_mem <- do.call(rbind, c(results_df_grouped_wide_mem_list, results2_df_grouped_wide_mem_list))

NumberSteps <- lapply(names(software_list), function(x){
    data.frame(software_group = x, NumberSteps = length(software_list[[x]]))
})


NumberSteps2 <- lapply(names(software_list2), function(x){
    data.frame(software_group = x, NumberSteps = length(software_list2[[x]]))
})

NumberSteps_df <- do.call(rbind, c(NumberSteps, NumberSteps2))



results_df_grouped_wide_mem <- inner_join(results_df_grouped_wide_mem, NumberSteps_df)

software_list_complete <- c("Demuxalot", "Demuxlet", "Dropulation", "Freemuxlet", "ScSplit", "Souporcell", "Vireo", "DoubletDecon", "DoubletDetection", "DoubletFinder", "ScDblFinder", "Scds", "Scrublet", "Solo")

results_df_grouped_mem <- pivot_longer(results_df_grouped_wide_mem, cols = c("Total_mem_G", "CPU_time_h","NumberSteps"), names_to = "Metric", values_to = "Value")
results_df_grouped_mem$software_group <- factor(results_df_grouped_mem$software_group, levels = software_list_complete)
results_df_grouped_mem$Metric <- factor(results_df_grouped_mem$Metric, levels = c("CPU_time_h", "Total_mem_G", "NumberSteps"))

results_df_grouped_mem$software_type <- ifelse(results_df_grouped_mem$software_group %in% c("Demuxalot", "Demuxlet", "Dropulation", "Freemuxlet", "ScSplit", "Souporcell", "Vireo"), "Demultiplexing", "Doublet Detecting")

results_df_grouped_mem <- unique(results_df_grouped_mem)


stat_box_sum <- function(y, upper_limit = max(results_df_grouped_mem$Metric)) {
  DF <- data.frame(
    y = max(y),
    label = round(median(y), 2)
  )
  DF
}

pCPUmem <- ggplot(results_df_grouped_mem, aes(software_group, Value)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    facet_wrap(Metric ~ software_type, scales = "free", switch = "y", nrow = 3,
            labeller = labeller("Metric" = c(CPU_time_h = "Time (h)", Total_mem_G_log = "Memory (GB)", NumberSteps = "Number of Computational Step"))) +
        ylab(NULL) +
        theme(strip.background = element_blank(),
            strip.placement = "outside") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_summary(
    fun.data = stat_box_sum, 
               geom = "text", 
    #hjust = 1,
    vjust = -0.5)

ggsave(pCPUmem, filename = paste0(outdir, "CPU_mem_usage.eps"))
ggsave(pCPUmem, filename = paste0(outdir, "CPU_mem_usage.pdf"), width = 6, height = 6)




results_df_grouped_wide_mem$CPU_time_h_log <- log(results_df_grouped_wide_mem$CPU_time_h)


results_df_grouped_mem_log <- pivot_longer(results_df_grouped_wide_mem, cols = c("Total_mem_G", "CPU_time_h_log","NumberSteps"), names_to = "Metric", values_to = "Value")
results_df_grouped_mem_log$software_group <- factor(results_df_grouped_mem_log$software_group, levels = software_list_complete)
results_df_grouped_mem_log$Metric <- factor(results_df_grouped_mem_log$Metric, levels = c("CPU_time_h_log", "Total_mem_G", "NumberSteps"))

results_df_grouped_mem_log$software_type <- ifelse(results_df_grouped_mem_log$software_group %in% c("Demuxalot", "Demuxlet", "Dropulation", "Freemuxlet", "ScSplit", "Souporcell", "Vireo"), "Demultiplexing", "Doublet Detecting")

results_df_grouped_mem_log <- unique(results_df_grouped_mem_log)


pCPUmem_log <- ggplot(results_df_grouped_mem_log, aes(software_group, Value)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    facet_wrap(Metric ~ software_type, scales = "free", switch = "y", nrow = 3,
            labeller = labeller("Metric" = c(CPU_time_h = "log10(Time (h))", Total_mem_G_log = "Memory (GB)", NumberSteps = "Number of Computational Step"))) +
        ylab(NULL) +
        theme(strip.background = element_blank(),
            strip.placement = "outside") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pCPUmem_log, filename = paste0(outdir, "CPU_mem_usage_log.eps"))
ggsave(pCPUmem_log, filename = paste0(outdir, "CPU_mem_usage_log.pdf"), width = 6, height = 6)



