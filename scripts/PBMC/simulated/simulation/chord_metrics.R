library(tidyverse)
library(plyr)
library(data.table)
library(future.apply)
suppressMessages(suppressWarnings(library(argparse)))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-p", "--pool", required = FALSE, type = "character", help = "the pool ID")

args <- parser$parse_args()


###### Set up Directories #####
##### Read in Results #####
doublets_dt <- fread(paste0(args$out, "/chord/chord_results_doublet.csv"), header = TRUE, col.names = c("row", "Barcode"))

doublets <- doublets_dt$Barcode


##### Read in the cell_info barcode files #####
simulated_barcodes <- fread(paste0(args$out,"/matrix_out/barcodes.tsv.gz"), sep = "\t", header = FALSE, col.names = c("Barcode"))

simulated_barcodes$DropletType <- ifelse(grepl(":", simulated_barcodes$Barcode), "doublet","singlet")
simulated_barcodes$DoubletType_DoubletType <- ifelse(gsub("[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", gsub(":.+", "", simulated_barcodes$Barcode)) == gsub(".+:[A,C,G,T,J,K,L,M,N,O,P,Q,R,S]+-", "", simulated_barcodes$Barcode), "homogenic_doublet", ifelse(grepl(":", simulated_barcodes$Barcode), "heterogenic_doublet","singlet"))


singlets <- simulated_barcodes$Barcode[!(simulated_barcodes$Barcode %in% doublets_dt$Barcode)]


print("metrics")

##### Calculate necessary metrics for individual softwares ######
## True singlet (TP)
## False singlet (FP)
## True doublet (TN)
## False doublet (FN)
## Accuracy
## Balanced Accuracy
## Sensitivity (TPR)
## Specificity (TNR)
## Precision
## MCC

# soft <- "DoubletDetection-scds-solo"
# pool <- "size128_SimulatedPool3"
# method <- "majority_doublet"

##### Set up dataframes #####
Metrics <- data.frame(matrix(nrow = 1, ncol = 13))
colnames(Metrics) <- c("software", "true_singlet", "false_singlet", "true_doublet", "false_doublet", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "npv", "ppv", "precision", "mcc")
Metrics$software <- "chord"
rownames(Metrics) <- "chord"
Metrics$pool <- args$pool

## True singlet (TP)
Metrics$true_singlet <- as.numeric(length(which(simulated_barcodes$DropletType == "singlet" & simulated_barcodes$Barcode %in% singlets)))

## False singlet (FP)
Metrics$false_singlet <- as.numeric(length(singlets) - Metrics$true_singlet)

### Update NA for pools where could not run (ie DoubletDecon)
if (Metrics$true_singlet == 0 & Metrics$false_singlet == 0){
    Metrics$true_doublet <- NA
    Metrics$false_doublet <- NA
    Metrics$accuracy <- NA
    Metrics$sensitivity <- NA
    Metrics$specificity <- NA
    Metrics$npv <- NA
    Metrics$ppv <- NA
    Metrics$balanced_accuracy <- NA
    Metrics$precision <- NA
    Metrics$mcc <- NA
} else{

    ## True doublet (TN)
    Metrics$true_doublet <-as.numeric(length(which(simulated_barcodes$DropletType == "doublet" & simulated_barcodes$Barcode %in% doublets)))


    ## False doublet (FN)
    Metrics$false_doublet <- as.numeric(length(doublets) - Metrics$true_doublet)


    ## Accuracy
    Metrics$accuracy <- as.numeric((Metrics$true_singlet + Metrics$true_doublet)/nrow(simulated_barcodes))


    ## Sensitivity (TPR)
    Metrics$sensitivity <- as.numeric(Metrics$true_singlet/length(which(simulated_barcodes$DropletType == "singlet")))


    ## Specificity (TNR)
    Metrics$specificity <- as.numeric(Metrics$true_doublet/length(which(simulated_barcodes$DropletType == "doublet")))


    ## Negative Predictive Value
    Metrics$npv <- as.numeric(Metrics$true_doublet/length(doublets))


    ## Positive Predictive Value
    Metrics$ppv <- as.numeric(Metrics$true_singlet/length(singlets))


    ## Balanced Accuracy
    Metrics$balanced_accuracy <- as.numeric((Metrics$sensitivity + Metrics$specificity)/2)


    ## Precision
    Metrics$precision <- as.numeric(Metrics$true_singlet/(Metrics$true_singlet + Metrics$false_singlet))

    
    ## MCC
    Metrics$mcc <- as.numeric(((Metrics$true_singlet * Metrics$true_doublet) - (Metrics$false_singlet * Metrics$false_doublet))/sqrt((Metrics$true_singlet + Metrics$false_singlet) * (Metrics$true_singlet + Metrics$false_doublet) * (Metrics$true_doublet + Metrics$false_singlet) * (Metrics$true_doublet + Metrics$false_doublet)))
}





rownames(Metrics) <- NULL
Metrics <- as.data.table(Metrics)

Metrics$size <- factor(as.numeric(gsub("size","", Metrics$pool) %>% gsub("_.+", "", .)), levels = c(2,4,8,16,32,64,128))
Metrics$mt_percent <- ifelse(grepl("mt", Metrics$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_mt_","", Metrics$pool) %>% gsub("pctl", "", .)), 0)
Metrics$ambient_percent <- ifelse(grepl("ambient", Metrics$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_ambient_","", Metrics$pool) %>% gsub("pctl", "", .)), 0)
Metrics$subsampled <- ifelse(grepl("subsampled", Metrics$pool), "subsampled", "normal")
Metrics$uneven <- ifelse(grepl("unevenN", Metrics$pool), as.numeric(gsub("size\\d+_SimulatedPool\\d_unevenN_","", Metrics$pool) %>% gsub("pctl", "", .)), 0)

fwrite(Metrics, paste0(args$out, "/chord_metrics.tsv"), sep = "\t")
