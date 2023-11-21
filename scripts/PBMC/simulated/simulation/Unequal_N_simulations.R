### Author: Drew Neavin
### Date: 21 October, 2020
### Reason: This script is to simulate different numbers of ambient RNA for pools -> 5% ,10% or 20% of droplets with high MT % (avg 30%, sd = 3%)
### Requirements: Must have samtools in environment when running


##### Load library paths
library(knitr)
library(tidyverse)
library(stringi)
library(data.table)



##### Load arguments from command line
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)

pool <- as.character(arguments[1,])
pctl <- as.numeric(as.character(arguments[2,]))
outdir <- as.character(arguments[3,])
threads <- as.numeric(as.character(arguments[4,]))
sample_file <- as.character(arguments[5,])
barcode_file <- as.character(arguments[6,])
out <- as.character(arguments[7,])


##### Make a dataframe with doublets for each set of singlets #####
Ndoublets_dt <- data.table("Ntotal" = seq(1:150000), "Ndoublets" = as.numeric(NA), "Nsinglets" = as.numeric(NA))
Ndoublets_dt$Ndoublets <- apply(Ndoublets_dt, 1, function(x) {
    return(round(min((0.2*x['Ntotal']), ((x['Ntotal']^2)*0.008)/1000)))
})
Ndoublets_dt$Nsinglets <- Ndoublets_dt$Ntotal - Ndoublets_dt$Ndoublets



##### 1. Read in the sample file and get the bams in a list
sample_meta <- fread(sample_file, sep = "\t")
bam_list <- strsplit(sample_meta[Pool == pool]$Bams, ",")[[1]]


##### 2. Get the number of UMIs per cell #####
barcodes <- fread(barcode_file, sep = "\t", header = FALSE, col.names = c("Barcodes"))


##### 3. Identify the number of droplets in the pool #####
number_droplets <- nrow(barcodes)



##### 4. Calculate the number of barcodes of the selected individual and rest of individuals
N_selected <- round(number_droplets*pctl)
N_rest <- number_droplets - N_selected



##### 5. Select individual with highest number of barcodes to be higher percentage individual (less confounding with the replicate droplets hopefully than one with less cells)
individual <- names(sort(table(gsub("[ACTG]+-", "", barcodes$Barcode)), decreasing = TRUE))[1]



##### 6. Pull cells from that individual
barcodes_assignment <- barcodes
barcodes_assignment$Individual <- gsub("[ACTG]+-", "", barcodes$Barcode)

barcodes_assignment_selected <- barcodes_assignment[Individual == individual]
barcodes_selected <- barcodes_assignment_selected[data.table(Barcodes = sample(barcodes_assignment_selected$Barcodes, N_selected, replace = TRUE)), on = "Barcodes"]



if (pctl == 0.95){
    print(paste0("pctl is ", pctl))
    N_selected <- nrow(barcodes_assignment[Individual == individual])


    if (!file.exists(paste0(outdir, "/barcodes4simulation.tsv"))){
        spikein_individual <- sample(names(sort(table(gsub("[ACTG]+-", "", barcodes$Barcode)), decreasing = TRUE))[2:length(names(sort(table(gsub("[ACTG]+-", "", barcodes$Barcode)), decreasing = TRUE)))], 1)


        ##### 7. Randomly pull cells from spikein individual
        barcodes_assignment_rest <- barcodes_assignment[Barcodes %in% sample(barcodes_assignment[Individual == spikein_individual]$Barcodes, round(0.05*(nrow(barcodes_assignment_selected)/0.95)), replace = FALSE)]
        
        
        ##### 8. Combine barcodes 
        barcodes_combined <- rbind(barcodes_selected, barcodes_assignment_rest)


        ##### 9. Write dataframes
        write_delim(data.frame("Barcodes" = barcodes_assignment_rest$Barcodes, "Individual" = barcodes_combined[Individual == spikein_individual]$Individual), paste0(outdir, "/barcode_list_",spikein_individual, ".tsv"), delim = "\t", col_names = FALSE)
        write_delim(data.frame("Barcodes" = barcodes_combined$Barcodes), paste0(outdir, "/barcodes4simulation.tsv"), delim = "\t", col_names = FALSE)

        dir.create(paste0(outdir, "/individual_separated"))

        saveRDS(spikein_individual, paste0(outdir, "/individual_separated/spikein_individual.rds"))


        N_dt <- data.table(Pool = pool, pctl = pctl, N = Ndoublets_dt[Nsinglets == nrow(barcodes_combined)]$Ndoublets[1])


    } else {
        barcodes_combined <- fread(paste0(outdir, "/barcodes4simulation.tsv"), sep = "\t", header = FALSE, col_names = c("Barcodes", "Individual"))
        spikein_individual <- readRDS(paste0(outdir, "/individual_separated/spikein_individual.rds"))
    }


    ##### 11. Use sinto to pull the spikein individual barcodes
    system(paste0("sinto filterbarcodes -b ", bam_list[grep(paste0(spikein_individual, "_updated.bam"), bam_list)], " --cells ", outdir, "/barcode_list_", spikein_individual, ".tsv --barcodetag CB --outdir ", outdir, "/individual_separated/ --nproc ", threads))



    ##### 12. Combine the two bams together, sort and index
    system(paste0("samtools merge ", outdir, "/unequalN.bam ", bam_list[grep(paste0(individual, "_updated.bam"), bam_list)], " ", outdir, "/individual_separated/", spikein_individual, ".bam"))

    
    ##### Remove used files
    system(paste0("rm ", paste0(outdir, "/individual_separated/", spikein_individual, ".bam"), collapse = " "))



} else {
    print(paste0("pctl is ", pctl))
    if (!file.exists(paste0(outdir, "/barcodes_selected_duplicated.tsv")) | !file.exists(paste0(outdir, "/barcodes4simulation.tsv"))){
        ##### 7. Identify barcodes duplicated, then replace with new barcodes
        barcodes_selected_duplicated <- barcodes_selected

        barcodes_selected_duplicated$Duplicated <- duplicated(barcodes_selected_duplicated$Barcodes)
        barcodes_selected_duplicated <- barcodes_selected_duplicated[order(Duplicated)]
        barcodes_selected_duplicated$Updated_Barcodes <-c(barcodes_selected_duplicated[Duplicated == FALSE]$Barcodes, paste0(stri_rand_strings(nrow(barcodes_selected_duplicated[Duplicated == TRUE]), 16, pattern = "[J-S]"), "-", unique(barcodes_selected_duplicated$Individual)))


        
        ##### 8. Randomly pull cells from other barcodes
        barcodes_assignment_rest <- barcodes_assignment[data.table(Barcodes = sample(barcodes_assignment[Individual != individual]$Barcodes, min(N_rest, nrow(barcodes_assignment[Individual != individual])), replace = FALSE)), on = "Barcodes"]



        ##### 9. Write barcode list files for sinto                                 
        ## Selected individual barcodes
        write_delim(data.frame("Barcodes" = barcodes_selected_duplicated[Duplicated == FALSE]$Updated_Barcodes, "Barcodes2" = barcodes_selected_duplicated[Duplicated == FALSE]$Updated_Barcodes), paste0(outdir, "/barcode_list_",individual, ".tsv"), delim = "\t", col_names = FALSE)
        
        ## Other individuals barcodes
        lapply(unique(barcodes_assignment_rest$Individual), function(individual_rest){
            write_delim(data.frame("Barcodes" = barcodes_assignment_rest[Individual == individual_rest]$Barcodes, "Individual" = barcodes_assignment_rest[Individual == individual_rest]$Individual), paste0(outdir, "/barcode_list_",individual_rest, ".tsv"), delim = "\t", col_names = FALSE)
        })


        N_dt <- data.table(Pool = pool, pctl = pctl, N = Ndoublets_dt[Nsinglets == length(c(barcodes_selected_duplicated$Updated_Barcodes, barcodes_assignment_rest$Barcodes))]$Ndoublets)

        ## Combined 
        write_delim(data.frame("Barcodes" = c(barcodes_selected_duplicated$Updated_Barcodes, barcodes_assignment_rest$Barcodes)), paste0(outdir, "/barcodes4simulation.tsv"), delim = "\t", col_names = FALSE)


        ## Dataframes needed if rerun from middle
        fwrite(barcodes_assignment_rest, paste0(outdir,"barcodes_assignment_rest.tsv"), sep = "\t")
        fwrite(barcodes_selected_duplicated, paste0(outdir,"barcodes_selected_duplicated.tsv"), sep = "\t")


    } else {
        ## Read in needed dataframes
        barcodes_assignment_rest <- fread(paste0(outdir,"barcodes_assignment_rest.tsv"), sep = "\t")
        barcodes_selected_duplicated <- fread(paste0(outdir,"barcodes_selected_duplicated.tsv"), sep = "\t")
    }


    ##### 10. Apply sinto to the selected individual to separate each barcode into its own bam
    dir.create(paste0(outdir, "/individual_separated"))

    if (!file.exists(paste0(barcodes_selected$Barcode[nrow(barcodes_selected)], ".bam"))){
        system(paste0("sinto filterbarcodes -b ", bam_list[grep(paste0(individual, "_updated.bam"), bam_list)], " --cells ", outdir, "/barcode_list_", individual, ".tsv --barcodetag CB --outdir ", outdir, "/individual_separated/ --nproc ", threads))
    }



    ##### 11. Apply sinto to each of the rest of the individuals in the pool separately to split the barcodes into a bam per individual
    ### Make a dataframe that has the selected individual and other barcodes for easy processing during sinto
    cells4separation <- barcodes_assignment_rest[, c("Barcodes", "Individual")]

    for (i in unique(cells4separation$Individual)){
        if (!file.exists(paste0(outdir, "/individual_separated/", i, ".bam"))){
            system(paste0("sinto filterbarcodes -b ", bam_list[grepl(paste0(i, "_updated.bam"), bam_list)], " --cells ", outdir, "/barcode_list_", i, ".tsv --barcodetag CB --outdir ", outdir, "/individual_separated/ --nproc ", threads))
        }
    }


    ##### 12. Merge the rest individuals bams together into one sam
    if (!file.exists(paste0(outdir, "/unequalN.bam"))){
        if (!file.exists(paste0(outdir, "/unequalN.sam"))){
            system(paste0("samtools merge -O SAM ", outdir, "/unequalN.sam ", paste(paste0(outdir, "/individual_separated/*.bam"), collapse = " ")))



            ##### 13. Copy the duplicated barcodes to the new barcode name and update name in file
            for (updated_barcode in barcodes_selected_duplicated[Duplicated == TRUE]$Updated_Barcodes){
                original_barcode <- barcodes_selected_duplicated[Updated_Barcodes == updated_barcode]$Barcodes
                system(paste0("samtools view ", outdir, "/individual_separated/", barcodes_selected_duplicated[Updated_Barcodes == updated_barcode]$Barcodes, ".bam | awk 'BEGIN{FS=\"\tCB:Z:[ACTG]+-[0-9_:]+\t\"}{print $1}' > ", outdir, "/individual_separated/", updated_barcode, "_part1.tsv"))
                system(paste0("samtools view ", outdir, "/individual_separated/", barcodes_selected_duplicated[Updated_Barcodes == updated_barcode]$Barcodes, ".bam | awk 'BEGIN{FS=\"\tCB:Z:[ACTG]+-[0-9_:]+\t\"}{print $2}' > ", outdir, "/individual_separated/", updated_barcode, "_part3.tsv"))
                f <- file(paste0(outdir, "/individual_separated/", updated_barcode, "_part1.tsv"), open="rb")
                n_count <- 0L
                while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
                    n_count <- n_count + sum(chunk == as.raw(10L))
                }
                close(f)
                fwrite(data.table(newBarcodes = paste0("CB:Z:", rep(updated_barcode,n_count))), paste0(outdir,"/individual_separated/", updated_barcode, "_part2.tsv"), col.names = FALSE)
                
                system(paste0("paste -d \"\t\" ", outdir, "/individual_separated/", updated_barcode, "_part1.tsv ", outdir,"/individual_separated/", updated_barcode, "_part2.tsv ", outdir, "/individual_separated/", updated_barcode, "_part3.tsv > ", outdir, "/individual_separated/", updated_barcode, ".tsv"))

                system(paste0("cat ", outdir, "/individual_separated/", updated_barcode, ".tsv >> ", outdir, "/unequalN.sam"))

                system(paste0("rm ", outdir, "/individual_separated/", updated_barcode, "_part1.tsv ",
                                        outdir, "/individual_separated/", updated_barcode, "_part2.tsv ",
                                        outdir, "/individual_separated/", updated_barcode, "_part3.tsv ",
                                        outdir, "/individual_separated/", updated_barcode, ".tsv"))
            }

        }
        ##### 14. convert to bam file and remove extra header lines in the process
        system(paste0("samtools view -H ", outdir, "/unequalN.sam | grep \"@HD\" > ", outdir, "/combined_singlets.sam"))
        system(paste0("samtools view -H ", outdir, "/unequalN.sam | grep \"@SQ\" >> ", outdir, "/combined_singlets.sam"))
        system(paste0("samtools view -H ", outdir, "/unequalN.sam | grep \"@RG\" | sort -u >> ", outdir, "/combined_singlets.sam"))
        system(paste0("samtools view -H ", outdir, "/unequalN.sam | grep \"@PG\" | sort -u >> ", outdir, "/combined_singlets.sam"))
        system(paste0("samtools view -H ", outdir, "/unequalN.sam | grep \"@CO\" | sort -u >> ", outdir, "/combined_singlets.sam"))
        system(paste0("samtools view ", outdir, "/unequalN.sam >> ", outdir, "/combined_singlets.sam"))


        system(paste0("samtools view -bS ", outdir, "/combined_singlets.sam > ", outdir, "/unequalN.bam"))

       
        ##### Remove some unnecessary files 
        system(paste0("rm ", outdir, "/unequalN.sam ", paste(paste0(outdir, "/individual_separated/", unique(barcodes_assignment_rest$Individual), ".bam"), collapse = " ")))
        system(paste0("rm ", outdir, "/individual_separated/*.bam"))
        system(paste0("rm ", outdir, "/combined_singlets.sam"))
    }

}


##### 15. sort bam #####
system(paste0("samtools sort -o ", outdir, "/combined_singlets.bam ", outdir, "/unequalN.bam"))
system(paste0("samtools index ", outdir, "/combined_singlets.bam"))


fwrite(N_dt, paste0(out,"/unequalN_doublets.tsv"), append = TRUE, sep = "\t")


print(Sys.time())

##### 16. Remove unnecessary files #####
system(paste0("rm ", outdir, "/unequalN.bam"))

