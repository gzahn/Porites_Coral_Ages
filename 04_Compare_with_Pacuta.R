# Combine data with 
# Packages
library(phyloseq)
library(tidyverse)
library(vegan)

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


# Remove low-abundance OTUs ####

# Look at taxa sums 
summary(taxa_sums(ps2))

# Keep OTUs with at least 10 occurrances
ps.min <- subset_taxa(ps2, taxa_sums(ps2) >= 10)



# Ordination ####

ord1=ordinate(ps.min,method = "NMDS",color="Species")
plot_ordination(ps.min,ord1,color="Species") + stat_ellipse()



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

heatmap(t(as.matrix(otu_table(ps_pa))))


# Calculate distance matrix
comm.dist <- vegdist(otu_table(ps_pa),method = "jaccard")

mds <- metaMDS(otu_table(ps_pa),distance = "jaccard")
plot(comm.dist)
stressplot(mds)

# PermANOVA
mod1 <- adonis(otu_table(ps) ~ sample_data(ps)$Species + sample_data(ps)$ReadDepth)
mod1
