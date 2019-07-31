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
