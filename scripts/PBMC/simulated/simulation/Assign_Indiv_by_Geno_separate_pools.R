.libPaths("/usr/local/lib/R/site-library")
library(tidyr)
library(tidyverse)
library(dplyr)
library(vcfR)
library(lsa)
library(ComplexHeatmap)

### Read in arguments
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
meta_file <- as.character(arguments[1,])
dir <- as.character(arguments[2,])
out <- as.character(arguments[3,])
pool <- as.character(arguments[4,])

dir <- paste0(dir,"/")
out <- paste0(out,"/")

print(out)


########## Set up functions ##########
calculate_DS <- function(GP_df){
    columns <- c()
    for (i in 1:ncol(GP_df)){
        columns <- c(columns, paste0(colnames(GP_df)[i],"-0"), paste0(colnames(GP_df)[i],"-1"), paste0(colnames(GP_df)[i],"-2"))
    }
    df <- GP_df
    colnames(df) <- paste0("c", colnames(df))
    colnames_orig <- colnames(df)
    for (i in 1:length(colnames_orig)){
        df <- separate(df, sep = ",", col = colnames_orig[i], into = columns[(1+(3*(i-1))):(3+(3*(i-1)))])
    }
    df <- mutate_all(df, function(x) as.numeric(as.character(x)))
    for (i in 1: ncol(GP_df)){
        GP_df[,i] <- df[,(2+((i-1)*3))] + 2* df[,(3+((i-1)*3))]
    }
    return(GP_df)
}

pearson_correlation <- function(df, ref_df, clust_df){
    for (col in colnames(df)){
        for (row in rownames(df)){
            df[row,col] <- cor(as.numeric(pull(ref_df, col)), as.numeric(pull(clust_df, row)), method = "pearson", use = "complete.obs")
        }
    }
    return(df)
}


if (!file.exists(paste0(out,"/reference_cluster_genotype_pearson_correlations.rds"))){

    print(paste0("Correlations for ", pool, " DO NOT yet exist. Starting generation."))
    ref_geno_list <- list()
    cluster_geno_list <- list()
    ref_geno_tidy_list <- list()
    cluster_geno_tidy_list <- list()

    for (software in c("freemuxlet", "scSplit", "souporcell")){
        if (software == "freemuxlet"){
            # Get file information
            file_info <- file.info(paste0(dir,pool,"/popscle/freemuxlet/Individual_genotypes_subset.vcf.gz"))

            # Check if the file is empty
            if (file_info$size == 0 | is.na(file_info$size)){
                print("The file is empty and skipping.")
            } else {
                ref_geno_list[[software]] <- read.vcfR(paste0(dir,pool,"/popscle/freemuxlet/Individual_genotypes_subset.vcf.gz"))
                cluster_geno_list[[software]] <- read.vcfR(paste0(dir,pool,"/popscle/freemuxlet/freemuxletOUT.clust1.vcf.gz"))
                ref_geno_tidy_list[[software]] <- as_tibble(extract.gt(element = "DS",ref_geno_list[[software]], IDtoRowNames =F))
                ref_geno_tidy_list[[software]]$ID <- paste0(ref_geno_list[[software]]@fix[,'CHROM'],":", ref_geno_list[[software]]@fix[,'POS'],"_", ref_geno_list[[software]]@fix[,'REF'], "_",ref_geno_list[[software]]@fix[,'ALT'])
                cluster_geno_tidy_list[[software]] <- as_tibble(extract.gt(element = "GP",cluster_geno_list[[software]], IDtoRowNames =F))
                cluster_geno_tidy_list[[software]] <- calculate_DS(cluster_geno_tidy_list[[software]])
                cluster_geno_tidy_list[[software]]$ID <- paste0(cluster_geno_list[[software]]@fix[,'CHROM'],":", cluster_geno_list[[software]]@fix[,'POS'],"_", cluster_geno_list[[software]]@fix[,'REF'], "_",cluster_geno_list[[software]]@fix[,'ALT'])
                cluster_geno_tidy_list[[software]] <- cluster_geno_tidy_list[[software]][colSums(!is.na(cluster_geno_tidy_list[[software]])) > 0]
                # cluster_geno_tidy_list[[software]] <- cluster_geno_tidy_list[[software]][complete.cases(cluster_geno_tidy_list[[software]]),]
            }
        } else if (software == "scSplit"){
            # Get file information
            file_info <- file.info(paste0(dir,pool,"/scSplit/scSplit.vcf"))

            # Check if the file is empty
            if (file_info$size == 0 | is.na(file_info$size)){
                print("The file is empty and skipping.")
            } else {
                ref_geno_list[[software]] <- read.vcfR(paste0(dir,pool,"/scSplit/Individual_genotypes_subset.vcf.gz"))
                cluster_geno_list[[software]] <- read.vcfR(paste0(dir,pool,"/scSplit/scSplit.vcf"))
                ref_geno_tidy_list[[software]] <- as_tibble(extract.gt(element = "GP", ref_geno_list[[software]], IDtoRowNames =F))
                ref_geno_tidy_list[[software]] <- calculate_DS(ref_geno_tidy_list[[software]])
                ref_geno_tidy_list[[software]]$ID <- paste0(ref_geno_list[[software]]@fix[,'CHROM'],":", ref_geno_list[[software]]@fix[,'POS'])
                ref_geno_tidy_list[[software]] <- ref_geno_tidy_list[[software]][!(ref_geno_tidy_list[[software]]$ID %in% ref_geno_tidy_list[[software]]$ID[duplicated(ref_geno_tidy_list[[software]]$ID)]),]
                cluster_geno_tidy_list[[software]] <- as_tibble(extract.gt(element = "GP",cluster_geno_list[[software]], IDtoRowNames =F))
                cluster_geno_tidy_list[[software]] <- calculate_DS(cluster_geno_tidy_list[[software]])
                cluster_geno_tidy_list[[software]]$ID <- paste0(cluster_geno_list[[software]]@fix[,'CHROM'],":", cluster_geno_list[[software]]@fix[,'POS'])
                cluster_geno_tidy_list[[software]] <- cluster_geno_tidy_list[[software]][colSums(!is.na(cluster_geno_tidy_list[[software]])) > 0]
                # cluster_geno_tidy_list[[software]] <- cluster_geno_tidy_list[[software]][complete.cases(cluster_geno_tidy_list[[software]]),]
            }
        } else if (software == "souporcell"){
            # Get file information
            file_info <- file.info(paste0(dir,pool,"/souporcell/cluster_genotypes.vcf"))

            # Check if the file is empty
            if (file_info$size == 0 | is.na(file_info$size)){
                print("The file is empty and skipping.")
            } else {
                ref_geno_list[[software]] <- read.vcfR(paste0(dir,pool,"/souporcell/Individual_genotypes_subset.vcf.gz"))
                ref_geno_tidy_list[[software]] <- as_tibble(extract.gt(element = "DS",ref_geno_list[[software]], IDtoRowNames =F))
                ref_geno_tidy_list[[software]]$ID <- paste0(ref_geno_list[[software]]@fix[,'CHROM'],":", ref_geno_list[[software]]@fix[,'POS'],"_", ref_geno_list[[software]]@fix[,'REF'], "_",ref_geno_list[[software]]@fix[,'ALT'])
                cluster_geno_list[[software]] <- read.vcfR(paste0(dir,pool,"/souporcell/cluster_genotypes.vcf"))
                cluster_geno_tidy_list[[software]] <- as_tibble(extract.gt(element = "GT",cluster_geno_list[[software]], IDtoRowNames =F))
                cluster_geno_tidy_list[[software]] <- as_tibble(lapply(cluster_geno_tidy_list[[software]], function(x) {gsub("0/0",0, x)}) %>%
                                            lapply(., function(x) {gsub("0/1",1, x)}) %>%
                                            lapply(., function(x) {gsub("1/0",1, x)}) %>%
                                            lapply(., function(x) {gsub("1/1",2, x)}))
                cluster_geno_tidy_list[[software]]$ID <- paste0(cluster_geno_list[[software]]@fix[,'CHROM'],":", cluster_geno_list[[software]]@fix[,'POS'],"_", cluster_geno_list[[software]]@fix[,'REF'], "_",cluster_geno_list[[software]]@fix[,'ALT'])
                cluster_geno_tidy_list[[software]] <- cluster_geno_tidy_list[[software]][colSums(!is.na(cluster_geno_tidy_list[[software]])) > 0]
                # cluster_geno_tidy_list[[software]] <- cluster_geno_tidy_list[[software]][complete.cases(cluster_geno_tidy_list[[software]]),]

            }
        } else {
            print("something went wrong and the software is not freemuxlet, scSplit or souporcell)?")
        }
    }

    ########## Get a unique list of SNPs that is in both the reference and cluster genotypes ##########
    methods <- names(ref_geno_tidy_list)
    print(methods == names(cluster_geno_tidy_list))

    locations_list <- lapply(methods, function(method){
        temp <- inner_join(ref_geno_tidy_list[[method]][,"ID"],cluster_geno_tidy_list[[method]][,"ID"])
        return(temp[!(temp$ID %in% temp[duplicated(temp),"ID"]),])
    })
    names(locations_list) <- methods


    ########## Keep just the SNPs that overlap ##########
    ref_geno_tidy_list <- lapply(names(ref_geno_tidy_list), function(x){
        left_join(locations_list[[x]], ref_geno_tidy_list[[x]])
    })
    names(ref_geno_tidy_list) <- methods

    cluster_geno_tidy_list <- lapply(names(cluster_geno_tidy_list), function(x){
        left_join(locations_list[[x]], cluster_geno_tidy_list[[x]])
    })
    names(cluster_geno_tidy_list) <- methods


    ########## Correlate all the cluster genotypes with the individuals genotyped ##########
    ##### Make a dataframe that has the clusters as the row names and the individuals as the column names #####
    pearson_correlations <- lapply(names(ref_geno_tidy_list), function(x){
        temp <- as.data.frame(matrix(nrow = (ncol(cluster_geno_tidy_list[[x]]) -1), ncol = (ncol(ref_geno_tidy_list[[x]]) -1)))
        colnames(temp) <- colnames(ref_geno_tidy_list[[x]])[2:(ncol(ref_geno_tidy_list[[x]]))]
        rownames(temp) <- colnames(cluster_geno_tidy_list[[x]])[2:(ncol(cluster_geno_tidy_list[[x]]))]
        temp <- pearson_correlation(temp, ref_geno_tidy_list[[x]], cluster_geno_tidy_list[[x]])

        return(temp)
    })
    names(pearson_correlations) <- methods


    ########## Save the correlation dataframes ##########
    saveRDS(pearson_correlations, file= paste0(out,"/reference_cluster_genotype_pearson_correlations.rds"))
    lapply(names(pearson_correlations), function(x){
        write_delim(pearson_correlations[[x]],path = paste0(out,x,"_pearson_correlations.tsv"), delim = "\t" )
    })
} else {
    print(paste0("Correlations for ", pool, " EXIST. Reading file."))
    print(paste0(out,"/reference_cluster_genotype_pearson_correlations.rds"))
    
    pearson_correlations <- readRDS(paste0(out,"/reference_cluster_genotype_pearson_correlations.rds"))

    methods <- names(pearson_correlations)
}


########## Create correlation figures ##########
print("Creating correlation figures")
col_fun = colorRampPalette(c("white", "red"))
pPearsonCorrelations_list <- lapply(names(pearson_correlations), FUN = function(x){
    Heatmap(as.matrix(pearson_correlations[[x]]), cluster_rows = T, col = col_fun(101), column_title = paste0(x))
}) 
names(pPearsonCorrelations_list) <- methods


########## Save the correlation figures ##########
print("Saving correlation figures")
lapply(names(pPearsonCorrelations_list), FUN = function(x){
    png(filename = paste0(out,"/",x,"_pearson_correlation.png"), width = 500)
        print(pPearsonCorrelations_list[[x]])
    dev.off()
})


########## Assign individual to cluster based on highest correlating individual ##########
key_list_no_filtering <- lapply(pearson_correlations, function(x){
    df <- as.data.frame(matrix(nrow = ncol(x), ncol = 3))
    colnames(df) <- c("Genotype_ID","Cluster_ID","Correlation")
    df$Genotype_ID <- colnames(x)
    for (id in df$Genotype_ID){
        df$Cluster_ID[which(df$Genotype_ID == id)] <- rownames(x)[which.max(x[,id])]
        df$Correlation[which(df$Genotype_ID == id)] <- max(x[,id])
    }
    return(df)
})

key_list <- lapply(pearson_correlations, function(x){
    df <- as.data.frame(matrix(nrow = ncol(x), ncol = 3))
    colnames(df) <- c("Genotype_ID","Cluster_ID","Correlation")
    df$Genotype_ID <- colnames(x)
    for (id in df$Genotype_ID){
        if (max(x[,id]) == max(x[rownames(x)[which.max(x[,id])],])){
            df$Cluster_ID[which(df$Genotype_ID == id)] <- rownames(x)[which.max(x[,id])]
            df$Correlation[which(df$Genotype_ID == id)] <- max(x[,id])
        } else {
            df <- df[-c(which(df$Genotype_ID == id)),]
        }
    }
    return(df)
})

key_list_no_filtering <- lapply(names(key_list_no_filtering), function(x){
    key_list_no_filtering[[x]]$Software <- x
    return(key_list_no_filtering[[x]])
})
names(key_list_no_filtering) <- methods

key_no_filtering <- do.call(rbind, key_list_no_filtering)

key_list <- lapply(names(key_list), function(x){
    key_list[[x]]$Software <- x
    return(key_list[[x]])
})
names(key_list) <- methods

key <- do.call(rbind, key_list)


saveRDS(key_no_filtering, paste0(out,"PoolKeys_preFiltering.rds"))
saveRDS(key, paste0(out,"PoolKeys.rds"))



########## Check that all assignments are unique ##########
print(key_list)

unique_assignments_list <- lapply(names(key_list), function(x){
    unique_assignments <- as.data.frame(matrix(ncol=3, nrow = 1))
    colnames(unique_assignments) <- c("Pool","Software","AllClusterIDsUnique")
    unique_assignments$Pool <- pool
    unique_assignments$Software <- x
    if(all(!duplicated(key_list[[x]]$Cluster_ID))){
        unique_assignments$AllClusterIDsUnique <- TRUE
    } else{
        unique_assignments$AllClusterIDsUnique <- FALSE
    }
    return(unique_assignments)
})


unique_assignments <- do.call(rbind, unique_assignments_list)

print(unique_assignments_list)
print(unique_assignments)

write_delim(unique_assignments, file = paste0(out,"UniqueClusterAssignmentCheck.tsv"), delim = "\t" )



########## Make a list of those that would have failed pre-filtering
unique_assignments_prefilter_list <- lapply(names(key_list_no_filtering), function(x){
    unique_assignments <- as.data.frame(matrix(ncol=3, nrow = 1))
    colnames(unique_assignments) <- c("Pool","Software","AllClusterIDsUnique")
    unique_assignments$Pool <- pool
    unique_assignments$Software <- x
    if(all(!duplicated(key_list_no_filtering[[x]]$Cluster_ID))){
        unique_assignments$AllClusterIDsUnique <- TRUE
    } else{
        unique_assignments$AllClusterIDsUnique <- FALSE
    }
    return(unique_assignments)
})


unique_assignments_prefilter <- do.call(rbind, unique_assignments_prefilter_list)
print(unique_assignments_prefilter)


write_delim(unique_assignments_prefilter, file = paste0(out,"UniqueClusterAssignmentCheck_prefiltering.tsv"), delim = "\t" )
