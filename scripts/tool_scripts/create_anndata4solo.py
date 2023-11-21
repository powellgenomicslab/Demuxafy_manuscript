#!/usr/bin/env python3

import anndata
import os
import argparse
import sys
import pandas as pd

# Load read10x function from mods directory
mods_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(mods_path)
import read10x

parser = argparse.ArgumentParser(
    description="Script to read in simulated 10x data and make anndata for solo.")
parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix (matrix.mtx)")
parser.add_argument("-b", "--barcodes_file", required = True, help = "cell ranger barcodes file")
parser.add_argument("-g", "--genes_file", required = True, help = "cell ranger genes file")
parser.add_argument("-d", "--matrix_dir", required = True, help = "matrix directory to save the anndata to")
args = parser.parse_args()

### Read in data ###
raw_counts = read10x.import_cellranger_mtx(args.counts_matrix)
barcodes_df = read10x.read_barcodes(args.barcodes_file)
genes_df = read10x.read_genes(args.genes_file)


### Create AnnData ###
data = anndata.AnnData(X = raw_counts, obs = barcodes_df, var = genes_df)


### Save AnnData ###
data.write(args.matrix_dir + "/anndata.h5ad")