library("scrunchy")


args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
dir <- arguments[1,]

print(dir)

umitools_to_mtx(paste0(dir,"/counts.tsv.gz"), output_path = paste0(dir,"/matrix_out/"))

