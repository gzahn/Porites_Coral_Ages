# --------------------------------------------------------------#
# Porites 16S Analyses - Depends on "01_Process_Raw_Reads.R"
# Author: Geoffrey Zahn
# --------------------------------------------------------------#

# Load packages and functions ####
library(phyloseq)
library(tidyverse)
library(vegan)

# functions
source("./R/plot_bar2.R")
source("./R/summarize_taxa_Joey711.R")
source("./R/heatmap_left.R")

# custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44","#b3bf2d",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000")
colorblindr::palette_plot(pal)


# Load cleaned data ####
ps <- readRDS("./output/phyloseq_object_16S.RDS")
meta = as.data.frame(sample_data(ps))

# Housekeeping ####
# Rename metadata columns
names(ps@sam_data) <- c("Multiplex ID","Library ID","SampleID","Location","Country","Species","CoralAge","CoralAgeBinned","Average_LE_mm","GPS","Control")
# Add Lat/Lon columns
LAT <- unlist(map(str_split(ps@sam_data$GPS,pattern = " "),1))
LON <- unlist(map(str_split(ps@sam_data$GPS,pattern = " "),2))
LAT <- str_remove(LAT,"N")
LON <- str_remove(LON,"E")
ps@sam_data$LAT <- LAT
ps@sam_data$LON <- LON

# Subset taxa to bacteria only (remove eukaryotes)
table(tax_table(ps)[,1])
ps <- subset_taxa(ps, Kingdom == "Bacteria")
colnames(tax_table(ps)) # No species-level assignments (see previous script)


# Normalize with relative abundance
ps_ra <- transform_sample_counts(ps, function(x) x / sum(x) )

# coral age vs community structure
# include location / size in model?
# mantel test to rule out? location as factor instead of age

dim(otu_table(ps_ra))
dim(sample_data(ps_ra))
adonis((otu_table(ps_ra)) ~ ps_ra@sam_data$CoralAgeBinned)
