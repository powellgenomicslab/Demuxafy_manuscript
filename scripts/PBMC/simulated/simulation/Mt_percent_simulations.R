### Author: Drew Neavin
### Date: 21 October, 2020
### Reason: This script is to simulate different numbers of ambient RNA for pools -> 5% ,10% or 20% of droplets with high MT % (avg 30%, sd = 3%)
### Requirements: Must have samtools in environment when running


##### Load library paths
library(knitr)
library(tidyverse)
library(stringi)
library(Seurat)
library(data.table)
library(parallel)


##### Load arguments from command line
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)

pool <- arguments[1,]
pctl <- arguments[2,]
pctl <- as.numeric(as.character(pctl))
dir <- arguments[3,]
outdir <- arguments[4,]
print(outdir)
matrix_dir <- arguments[5,]
matrix_dir <- as.character(matrix_dir)
print(matrix_dir)
bamfile <- arguments[6,]
threads <- as.numeric(as.character(arguments[7,]))



##### 1. Get the number of UMIs per cell #####
counts <- Read10X(matrix_dir, gene.column = 1)
UMIs <- data.frame("Barcodes" = colnames(counts), "UMIs" = Matrix::colSums(counts, sparseResult = FALSE))



##### 2. Identify the number of droplets in the pool #####
number_droplets <- nrow(UMIs)



##### 3. Calculate the number of barcodes to increase Mt% in
number_barcodes <- round(number_droplets*pctl)


if (!file.exists(paste0(outdir, "/barcode_list.tsv"))){
    ##### 4. Randomly sample barcodes to have increased Mt% and generate normal distribution to assign mt% to each barcode
    selected_barcodes <- data.frame("Barcodes" = sample(UMIs$Barcodes, number_barcodes), "Mt_percent" = rnorm(number_barcodes, 30, 3))



    ##### 5 write the barcodes as a file to be used as input for grep to pull just these barcodes from MT file
    write_delim(data.frame("Barcodes" = selected_barcodes$Barcodes, "Barcodes2" = selected_barcodes$Barcodes), paste0(outdir, "/barcode_list.tsv"), delim = "\t", col_names = FALSE)




    ##### Split barcodes into smaller numbers if more than 3000 barcodes #####
    if (length(selected_barcodes$Barcodes) > 3000){
        for (i in 1:ceiling(length(selected_barcodes$Barcodes)/3000)){
            write_delim(data.frame("Barcodes" = selected_barcodes$Barcodes[(1+(3000*(i-1))):min((3000*i), length(selected_barcodes$Barcodes))], "Barcodes2" = selected_barcodes$Barcodes[(1+(3000*(i-1))):min((3000*i), length(selected_barcodes$Barcodes))]), paste0(outdir, "/barcode_list_",i,".tsv"), delim = "\t", col_names = FALSE)
        }
    }





    ##### 6. Use bash samtools to pull the MT reads #####
    system(paste0("samtools view -b ",bamfile, " MT > ", outdir, "/MT_reads.bam"))
    system(paste0("samtools index ",outdir, "/MT_reads.bam"))
} else {
    selected_barcodes <- fread(paste0(outdir, "/barcode_list.tsv"), header = FALSE, col.names = c("Barcodes", "Barcodes2"))
    selected_barcodes$Barcodes2 <- NULL
}




##### 7. Pull just barcodes of interest from Mt sam file
dir.create(paste0(outdir, "/high_mt_barcodes"))
if (length(selected_barcodes$Barcodes) > 3000){
    for (i in 1:ceiling(length(selected_barcodes$Barcodes)/3000)){
        print(i)
        if (!file.exists(paste0(outdir, "/high_mt_barcodes/", selected_barcodes$Barcodes[min(length(selected_barcodes$Barcodes), i*3000)], ".bam"))){
            # if (file.exists(paste0(outdir, "/high_mt_barcodes/", selected_barcodes$Barcodes[(1+(3000*(i-1)))], ".+"))){
            if (length(list.files(paste0(outdir, "/high_mt_barcodes/"), pattern = selected_barcodes$Barcodes[min((3000*i), length(selected_barcodes$Barcodes))])) > 0) {
                print(paste0("Removing files since last ", i, " run wasn't completed."))
                for (c in selected_barcodes$Barcodes[(1+(3000*(i-1))):min((3000*i), length(selected_barcodes$Barcodes))]){
                    system(paste0("rm ", paste(paste0(outdir, "/high_mt_barcodes/", c, "*"), collapse = " ")))
                }
            }
            system(paste0("sinto filterbarcodes -b ", outdir, "/MT_reads.bam --cells ", outdir, "/barcode_list_", i, ".tsv --barcodetag CB --outdir ", outdir, "/high_mt_barcodes/ --nproc ", threads))
        }
    }
} else {
    system(paste0("sinto filterbarcodes -b ", outdir, "/MT_reads.bam --cells ", outdir, "/barcode_list.tsv --barcodetag CB --outdir ", outdir, "/high_mt_barcodes/ --nproc ", threads))
}


##### 8. copy the normal bam to the high mt bam
system(paste0("samtools view -H ", bamfile, " > ", outdir, "/high_mt.sam"))
system(paste0("samtools view ", bamfile, " >> ", outdir, "/high_mt.sam"))



##### 9. Identify the number of lines in the bam file to identify the total number of reads to possibly pull from using bash
for (barcode in selected_barcodes$Barcodes) {
    system(paste0("samtools view ",outdir, "/high_mt_barcodes/", barcode, ".bam | wc -l | awk '{print $1}' >", outdir, "/high_mt_barcodes/number_mt_lines_", barcode, ".tsv"))
}

mt_sam_length <- lapply(selected_barcodes$Barcodes, function(x){
    read_delim(paste0(outdir,"/high_mt_barcodes/number_mt_lines_", x, ".tsv"), delim ="\t", col_names = c("N_lines"))
})
names(mt_sam_length) <- selected_barcodes$Barcodes





##### 11. Pull randomly chosen lines to new file and replace UMIs with new codes
lapply(names(mt_sam_length), function(x){
    print(x)
    fwrite(data.table(newUMIs = paste0("UB:Z:", stri_rand_strings(mt_sam_length[[x]]$N_lines, 10, pattern = "[J-S]"))), paste0(outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part2.tsv"), col.names = FALSE)
    system(paste0("samtools view ", outdir,"/high_mt_barcodes/", x, ".bam | shuf -r -n ", mt_sam_length[[x]]$N_lines, " > ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add.sam"))
    system(paste0("awk 'BEGIN{FS=\"\tUB:Z:\"}{print $1}' ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add.sam > ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part1.tsv"))
    system(paste0("awk 'BEGIN{FS=\"\tUB:Z:[ACTG]+-[0-9_:]+\t\"}{print $2}' ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add.sam > ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part3.tsv"))
    system(paste0("paste -d \"\t\" ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part1.tsv ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part2.tsv ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part3.tsv > ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_updated.tsv"))
    ##### Add Mt edited lines to high mt file
    system(paste0("cat ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_updated.tsv >> ", outdir, "/high_mt.sam"))
})

lapply(names(mt_sam_length), function(x){
    ##### Delete files not needed anymore
    system(paste0("rm ", outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part1.tsv ", 
                        outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part2.tsv ", 
                        outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add_part3.tsv ", 
                        outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_add.sam ",  
                        outdir,"/high_mt_barcodes/number_mt_lines_", x, ".tsv ",
                        outdir, "/high_mt_barcodes/barcode_", x, "_MT_reads_2_updated.tsv"))
})

##### 12. convert sam to bam and sort #####
system(paste0("samtools view -b ", outdir, "/high_mt.sam > ", outdir, "/high_mt.bam"))
system(paste0("samtools sort -o ", outdir, "/pooled.sorted.bam ", outdir, "/high_mt.bam"))
system(paste0("samtools index ", outdir, "/pooled.sorted.bam"))

print(Sys.time())

##### 13. Remove unnecessary files #####
system(paste0("rm ", outdir, "/high_mt.bam ", 
                        outdir, "/high_mt.sam ", 
                        outdir, "/MT_reads.bam ", 
                        outdir, "/MT_reads.bam.bai "))







