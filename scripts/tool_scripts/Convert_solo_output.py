#!/usr/bin/env python

##### Set up the packages #####
import argparse
import pandas as pd
import numpy
import sys
import os

mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods" ## Do not change - this is the path to the mods folder in the singularity image with custom script for loading 10x data in python
sys.path.append(mods_path)
import read10x

##### Parse the variables passed to python #####
parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
parser.add_argument("-b", "--barcodes", required = True, help = "barcodes.tsv or barcodes.tsv.gz from cellranger")
parser.add_argument("-s", "--solo_output", required = True, help = "solo output directory")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory")
args = parser.parse_args()

##### Read in the barcodes and the solo results #####
barcodes_df = read10x.read_barcodes(args.barcodes)
doublet_list = numpy.load(args.solo_output + "/is_doublet.npy")
scores_list = numpy.load(args.solo_output + "/logit_scores.npy")

##### make a final dataframe of results + barcodes #####
dataframe = barcodes_df
dataframe["solo_DropletType"] = doublet_list
dataframe["solo_DropletScore"] = scores_list

##### Replace True and False with singlet and doublet #####
dataframe.solo_DropletType = dataframe.solo_DropletType.replace(True, "doublet")
dataframe.solo_DropletType = dataframe.solo_DropletType.replace(False, "singlet")

##### Write results
dataframe.to_csv(os.path.join(args.outdir,'solo_results.txt'), sep = "\t", index = False)
