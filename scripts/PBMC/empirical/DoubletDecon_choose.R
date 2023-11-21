##### Load in libraries #####
library(tidyverse)


##### Set up directories #####
dir <- "/path/to/output/fibroblasts/"
meta <- read_delim("/path/to/Demuxafy_manuscript/files/PBMC/PBMC_sample_meta.tsv", delim = "\t")
pools <- meta$Pool



##### Read in data #####
DoubletDecon_singlets <- lapply(pools, function(x){
    lapply(seq(0.4,1.1,0.1), function(y){
        read_delim(paste0(dir,x,"/DoubletDecon_rhop",y,"/Final_nondoublets_groups_DoubletDecon_results.txt"), delim ="\t")
    })
})

DoubletDecon_singlets <- lapply(DoubletDecon_singlets, function(x){
    names(x) <- paste0("DoubletDecon_rhop",seq(0.4,1.1,0.1))
    return(x)
})
names(DoubletDecon_singlets) <- pools


barcodes <- lapply(pools, function(x){
    read_delim(paste0(dir,x,"/matrix_out/barcodes.tsv"), delim = "\t",  col_names = FALSE)
})
names(barcodes) <- pools




##### Calculate expected singlets #####
expected_singlets_number <- lapply(barcodes, function(x){
    nrow(x) - (((nrow(x)^2)*0.008)/1000)
})



##### Get actual predicted doublet numbers for each condition #####
DoubletDecon_singlets_number <- lapply(DoubletDecon_singlets, function(x){
    lapply(x, function(y){
        nrow(y)
    })
})

DoubletDecon_singlets_number_dfs <- lapply(DoubletDecon_singlets_number, function(x){
    do.call(rbind,x)
})

DoubletDecon_singlets_number_dfs <- lapply(names(DoubletDecon_singlets_number_dfs), function(x){
    colnames(DoubletDecon_singlets_number_dfs[[x]]) <- c("Singlets")
    DoubletDecon_singlets_number_dfs[[x]] <- as.data.frame(DoubletDecon_singlets_number_dfs[[x]])
    DoubletDecon_singlets_number_dfs[[x]]$Difference <- abs(DoubletDecon_singlets_number_dfs[[x]]$Singlets - expected_singlets_number[[x]])
    return(DoubletDecon_singlets_number_dfs[[x]])
})


##### Choose best condition based on closest to expected number of doublets #####
best <- lapply(DoubletDecon_singlets_number_dfs, function(x){
    rownames(x)[which.min(x$Difference)]
})
names(best) <- pools

best_df <- as.data.frame(do.call(rbind, best))
colnames(best_df) <- "rhop"
best_df$rhop <- gsub("DoubletDecon_rhop", "", best_df$rhop)
best_df$Pool <- rownames(best_df)

write_delim(best_df, paste0(dir,"DoubletDecon_rhops.tsv"), delim = "\t")
