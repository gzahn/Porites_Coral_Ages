# --------------------------------------------------------------#
# Porites 16S Analyses - Depends on "01_Process_Raw_Reads.R"
# and "01_Cleanup_Processed_Data.R"
#
# Combine P luteus data with previously published P. acuta 
# microbiome data from same sites - Define and investigate
# core bacterial microbiome of both host corals
#
# Author: Geoffrey Zahn
# --------------------------------------------------------------#



# Packages
library(phyloseq)
library(tidyverse)
library(vegan)
# library(devtools) # Load the devtools package
# install_github("microbiome/microbiome") # Install the package
library(microbiome)
library(RColorBrewer)

# functions
source("./R/plot_bar2.R")


# Load data ####
pluteus <- readRDS("./output/phyloseq_object_16S_cleaned.RDS")
pacuta <- readRDS("./output/Pacuta_phyloseq_object.Rds")


# fix column names
pacuta@sam_data$Location <- as.character(pacuta@sam_data$Island.Collected.From)
pacuta@sam_data$Library.ID <-  as.character(pacuta@sam_data$SampleName)

# Parse components for merging
pacuta.otu <- otu_table(pacuta)
pacuta.tax <- tax_table(pacuta)
pacuta.sam <- as(sample_data(pacuta), "data.frame")
pluteus.otu <- otu_table(pluteus) 
pluteus.tax <- tax_table(pluteus)
pluteus.sam <- as(sample_data(pluteus),"data.frame")

# Merge phyloseq objects ####
ps <- merge_phyloseq(pluteus.otu,sample_data(pluteus.sam),pluteus.tax,pacuta.otu,sample_data(pacuta.sam),pacuta.tax)

# Fix metadata
ps@sam_data$Species[ps@sam_data$Species=="1"] <- "P acuta"

ggplot(mapping=aes(x=1:nrow(otu_table(ps)),y=rowSums(otu_table(ps)))) + geom_point(aes(color=ps@sam_data$Species)) + 
  labs(x="Sample number",y="Read count",color="Species") + theme_bw()

# Add read count to sampledata
sample_data(ps)$ReadDepth <- rowSums(otu_table(ps))



# Collapse based on taxonomy at genus level ####
ps2 = tax_glom(ps,taxrank = rank_names(ps)[6])
?tax_glom

# Remove low-abundance OTUs ####

# Look at taxa sums 
summary(taxa_sums(ps2))

# Keep OTUs with at least 10 occurrances
ps.min <- subset_taxa(ps2, taxa_sums(ps2) >= 10)



# Ordination ####

ord1=ordinate(ps.min,method = "NMDS",color="Species")
plot_ordination(ps.min,ord1,color="Species") + stat_ellipse()
ggsave("./output/Pacuta_vs_Plutea_NMDSPlot_draft.png")


mds <- metaMDS(otu_table(ps.min))
stressplot(mds)
plot(nmds$points[,1])
mds1 <- nmds$points[,1]
mds2 <- nmds$points[,2]

ord.df <- data.frame(MDS1=mds1,MDS2=mds2,Location=sample_data(ps.min)$Location,HostSpecies=sample_data(ps.min)$Species)

ggplot(ord.df,aes(x=MDS1,y=MDS2,color=HostSpecies)) + geom_point() + stat_ellipse()


ord1 = ordinate(ps.min, method = "NMDS")
plot_ordination(ps.min,ord1,color="Species")

# Calculate distance matrix
comm.dist <- vegdist(otu_table(ps.min),method = "bray")

plot(comm.dist)

# Convert to presence-absence ####
ps_pa <- transform_sample_counts(ps.min, function(abund) 1*(abund>0))


# prepare data for heatmap and plot ####
pa = as(t(otu_table(ps_pa)),"matrix")

cols = plyr::mapvalues(ps_pa@sam_data$Species,from=unique(ps_pa@sam_data$Species),to=c("Blue","Red"))

heatmap(pa,Rowv = NA,ColSideColors = cols,Colv = NA)

# Blue = P lutea, Red = P acuta
heatmap(t(as.matrix(otu_table(ps_pa))),Rowv = NA,ColSideColors = cols,Colv = NA, labRow = NA, col = gray.colors(2))


# Find genera tha overlap between species of corals ####

Plutea.cols = grep("ABB",x = colnames(pa))
Pacuta.cols = grep("AOO",x = colnames(pa))
Plutea.matrix = pa[,Plutea.cols]
Pacuta.matrix = pa[,Pacuta.cols]

shared.genera = rowSums(Plutea.matrix) > 0 & rowSums(Pacuta.matrix) > 0
shared.genera = which(shared.genera == TRUE)
sink("./output/shared_Genera_Pactua-Plutea.txt")
ps_pa@tax_table[shared.genera,]
sink(NULL)


# Prep shared taxa data frame
shared.genera.names = c(ps_pa@tax_table[shared.genera,"Genus"])
shared.family.names = c(ps_pa@tax_table[shared.genera,"Family"])
shared.order.names = c(ps_pa@tax_table[shared.genera,"Order"])
shared.phylum.names = c(ps_pa@tax_table[shared.genera,"Phylum"])

df.shared.taxa <- data.frame(Genus = shared.genera.names,Family = shared.family.names,Order=shared.order.names,Phylum=shared.phylum.names)

# Bar plot, colored by phylum
ggplot(df.shared.taxa) + geom_bar(aes(x=reorder(Order,Order,function(x)-length(x)),fill=Phylum)) +
  theme(axis.text.x = element_text(angle=90),
        axis.title = element_text(face="bold")) + 
  labs(x="Order",y="Shared Count")
ggsave("./output/figs/Shared_Genera_Pacuta-Plutea_within_each-order.png",dpi=300)



# PermANOVA ####
mod1 <- adonis(otu_table(ps) ~ sample_data(ps)$Species + sample_data(ps)$ReadDepth)
mod1



# Core Microbiome(s) ####
pluteus # P. luteus ... this study
pacuta # P. acuta ... previous study
ps # Combined Pluteus and Pacuta
ps2 # combined and taxa merged to genus level


# Let's start with a heatmap of all genera relative abundance in full merged dataset
ps2_ra <- transform_sample_counts(ps2,function(x) x/sum(x))
cols = plyr::mapvalues(ps2@sam_data$Species,from=unique(ps2@sam_data$Species),to=c("Blue","Red"))

    # Blue = P luteus, Red = P acuta
    heatmap((as.matrix(otu_table(ps2_ra))),Rowv = NA,RowSideColors = cols,Colv = NA, labRow = NA, col = gray.colors(20),labCol = NA)

# ... Not very helpful...try to find the core of each host, then merge

# Compositional versions of the data sets, individually, and combined...
ps2.rel <- microbiome::transform(ps2, "clr")
head(prevalence(ps2.rel, detection = 0.01, sort = TRUE))

ps.rel <- microbiome::transform(ps, "compositional")
pluteus.rel <- microbiome::transform(pluteus, "compositional")
pacuta.rel <- microbiome::transform(pacuta, "compositional")

# tax-glommed to genus level
pluteus2 <- tax_glom(pluteus,taxrank = rank_names(ps)[6])
pluteus2.rel <- comp_transform(pluteus2, "compositional")

pacuta2 <- tax_glom(pacuta,taxrank = rank_names(ps)[6])
pacuta2.rel <- microbiome::transform(pacuta2, "compositional")

    # Examine core thresholds:
prevalences <- seq(.05, .5, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Investigate visually
plot_core(ps2.rel, prevalences = prevalences, detections = detections, plot.type = "lineplot") + xlab("Relative Abundance (%)")
ggsave("./output/figs/core_lineplot_thresholds_full.png")

plot_core(ps2.rel, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, horizontal = TRUE) +
  theme(axis.text.x = element_blank())
ggsave("./output/figs/core_heatmap_thresholds_full.png", dpi=300)

plot_core(pacuta2.rel, prevalences = prevalences, detections = detections, plot.type = "lineplot") + xlab("Relative Abundance (%)")



# Define core as prevalence = 0.5, detection threshold = 0.1% relabund
ps2.core <- core(ps2.rel, detection = 0.001,prevalence = 0.5, include.lowest = TRUE)

# Same for separate host coral species...and for non-glommed full taxa
pluteus.core <- core(pluteus2.rel, detection = 0.001,prevalence = 0.5, include.lowest = TRUE)
pacuta.core <- core(pacuta2.rel, detection = 0.001,prevalence = 0.5, include.lowest = TRUE)
    # remove phy_tree for convenience
pacuta.core <- phyloseq(otu_table(pacuta.core),tax_table(pacuta.core),sample_data(pacuta.core))

full_core <- core(ps.rel, detection = 0.001,prevalence = 0.5, include.lowest = TRUE) # error

# No shared core taxa (ESVs) between both species of corals unless merged to genus level!

tax <- merge_phyloseq_pair(pacuta.core@tax_table,pluteus.core@tax_table)
otu <- merge_phyloseq_pair(pacuta.core@otu_table,pluteus.core@otu_table)



# Make new sample data frame
Species = c(as.character(pacuta.core@sam_data$Species), pluteus.core@sam_data$Species)
Location = c(pacuta.core@sam_data$Location,pluteus.core@sam_data$Location)
SampleNames = c(as.character(pacuta.core@sam_data$SampleID), pluteus@sam_data$SampleID)
IDs = c(as.character(pacuta.core@sam_data$SampleName),as.character(pluteus.core@sam_data$`Library ID`))

sam = data.frame(SampleID,Species,Location,row.names = IDs)

# Merge phyloseq objects manually
ps.core_ra <- phyloseq(otu_table(otu),tax_table(tax),sample_data(sam))

# Merge samples based on location AND species
newvar = paste0(ps.core_ra@sam_data$Species,"_",ps.core_ra@sam_data$Location)
ps.core_ra@sam_data$NewVar <- newvar
ps.core.location <- merge_samples(ps.core_ra,group = "NewVar",fun = "mean")
# repair merged sample data
spp = unlist(map(str_split(sample_names(ps.core.location),"_"),1))
loc = unlist(map(str_split(sample_names(ps.core.location),"_"),2))

ps.core.location@sam_data$Location <- loc
ps.core.location@sam_data$Species <- spp

ps.core.location.ra <- transform_sample_counts(ps.core.location,function(x) x/sum(x))
plot_bar(ps.core.location.ra,fill="Family") + facet_wrap(~Location)

# remove locations where both species weren't collected
ps.core.compare <- subset_samples(ps.core.location.ra, !Location %in% c("St John","Tanah Merah","TPT"))

# Barplot
plot_bar(ps.core.compare,x="sample_Species",fill="Class") + facet_grid(~Location) +
  theme_bw() + theme(axis.text.x = element_text(face = "italic",angle=90)) +
  labs(x="Host species",y="Relative abundance")
ggsave("./output/figs/Core_Class_relabund_by_coral-species_and_Location.png",dpi=300,width = 10,height = 6)

# Heatmap of core genera of both host coral species
sp.cols <- c(rep("Red",97),rep("Blue",160)) # Blue = P luteus, Red = P acuta
heatmap(t(as(otu,"matrix")),col=gray.colors(10),
        Rowv = NA,ColSideColors = sp.cols,Colv = NA,labCol = NA,
        labRow = tax[,6],margins = c(5,5))
# No real overlap in core microbiome between coral hosts, as defined here!



# Should I build core microbiome for each sampling location???