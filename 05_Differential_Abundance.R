# --------------------------------------------------------------#
# Porites 16S Analyses - Depends on "01_Process_Raw_Reads.R"
# and "01_Cleanup_Processed_Data.R"
#
# Diffrential abundance and dispersion of taxa based on CoralAge
#
# Author: Geoffrey Zahn
# --------------------------------------------------------------#


# Packages, functions, data ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(corncob)

# custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44","#b3bf2d",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000")

# Read cleaned data
ps_ra <- readRDS("./output/phyloseq_cleaned_relabund.RDS")
