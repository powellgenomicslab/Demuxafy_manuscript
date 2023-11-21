### Author: Drew Neavin
### Date: 5 November, 2020
### Reason: This script is to simulate different numbers of mt % RNA for pools
### Requirements: Must have samtools in environment when running

library(knitr)
library(tidyverse)
library(stringi)
library(data.table)


args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)

pool <- arguments[1,]
pctl <- arguments[2,]
pctl <- as.numeric(as.character(pctl))
dir <- arguments[3,]
outdir <- arguments[4,]
matrix_dir <- arguments[5,]
matrix_dir <- as.character(matrix_dir)
print(matrix_dir)
bamfile <- arguments[6,]



##### 1. Get the number of UMIs per cell #####
counts <- Read10X(matrix_dir, gene.column = 1)
UMIs <- data.frame("Barcodes" = colnames(counts), "UMIs" = Matrix::colSums(counts, sparseResult = FALSE))



##### 2. Identify the number of droplets in the pool #####
number_droplets <- nrow(UMIs)



##### 3. Calculate the number of reads to reassign to simulate various ambient RNA percentages
number_ambient <- round(sum(UMIs$UMIs)*pctl)



##### 4. Randomly choose lines to pull reads from (number = total to pull, range = 1:length of bam)
lines <- sample(1:sum(UMIs$UMIs), number_ambient, replace = TRUE)



##### 5. use bash to copy bam file to sam
system(paste0("samtools view -H ", bamfile, " > ", outdir,"/ambient_rna.sam"))
system(paste0("samtools view ", bamfile, " >> ", outdir,"/ambient_rna.sam"))



##### 6. Get list of random barcodes to replace original barcodes
new_barcodes <- paste0("CB:Z:",UMIs$Barcodes[sample(1:nrow(UMIs), number_ambient, replace = TRUE)])
new_umis <- paste0("UB:Z:", stri_rand_strings(number_ambient, 10, pattern = "[J-S]"))


##### 7. Combine the UMIs and barcodes
new_dt <- data.table(UMIs = new_umis, Barcodes = new_barcodes)

fwrite(new_dt, paste0(outdir,"/ambient_2_add_part2.tsv"), sep = "\t", col.names = FALSE)



##### 8. Get the lines that will be given random ambient assignments, too long to be done in one go, loop it
system(paste0("samtools view ", bamfile, " | shuf -r -n ", number_ambient, " > ", outdir, "/ambient_2_add.tsv"))


##### 9. Get lines before UMIs and Barcodes and after
## Later realized the doublets are ordered differently so need to change this
system(paste0("awk 'BEGIN{FS=\"\tCB:Z:[A-Z0-9:_-]+\t\"}{print $1}' ", outdir, "/ambient_2_add.tsv > ", outdir, "/ambient_2_add_part1.tsv"))
system(paste0("awk 'BEGIN{FS=\"\tCB:Z:[A-Z0-9:_-]+\t\"}{print $2}' ", outdir, "/ambient_2_add.tsv > ", outdir, "/ambient_2_add_part3.tsv"))

##### 10. Paste together parts
system(paste0("paste -d \"\t\" ", outdir, "/ambient_2_add_part1.tsv ", outdir, "/ambient_2_add_part2.tsv ", outdir, "/ambient_2_add_part3.tsv > ", outdir, "/ambient_2_add_updated.tsv"))



##### 11. Add the new ambient reads to the file
system(paste0("cat ", outdir, "/ambient_2_add_updated.tsv >> ", outdir, "/ambient_rna.sam"))



##### 12. convert sam to bam and sort #####
system(paste0("samtools view -b ", outdir, "/ambient_rna.sam > ", outdir, "/ambient_rna.bam"))
system(paste0("samtools sort -o ", outdir, "/pooled.sorted.bam ", outdir, "/ambient_rna.bam"))
system(paste0("samtools index ", outdir, "/pooled.sorted.bam"))



##### 12. Remove unnecessary files #####
system(paste0("rm ", outdir, "/ambient_rna.bam ",
                        outdir, "/ambient_rna.sam ",
                        outdir, "/ambient_2_add_updated.tsv ", 
                        outdir, "/ambient_2_add_part3.tsv ", 
                        outdir, "/ambient_2_add_part1.tsv ", 
                        outdir, "/ambient_2_add.tsv ", 
                        outdir, "/ambient_2_add_part2.tsv"))
