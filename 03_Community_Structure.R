# --------------------------------------------------------------#
# Porites 16S Analyses - Depends on "01_Process_Raw_Reads.R"
# Author: Geoffrey Zahn
# --------------------------------------------------------------#

# Load packages, functions, and data ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(corrplot)
library(ecodist)
library(ade4)
library(splines)
library(modelr)
library(lme4)
library(patchwork)
library(igraph)

# functions
source("./R/plot_bar2.R")
source("./R/summarize_taxa_Joey711.R")
source("./R/heatmap_left.R")

# custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44","#b3bf2d",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000")

# Read cleaned data
ps_ra <- readRDS("./output/phyloseq_cleaned_relabund.RDS")

# convert OTU table and metadata to data.frame for easier downstream access
otu = as.data.frame(as(otu_table(ps_ra),"matrix"))
meta = as.data.frame(sample_data(ps_ra))

# Calculate community distance measures ####
bray = vegdist(otu,method = "bray")
jaccard = vegdist(otu,method = "jaccard", binary = TRUE)

# Plot community dist vs age dist ####
age <- dist(meta$CoralAge)
ggplot(mapping=aes(x=age,y=bray)) + geom_point(alpha=.2) + 
  geom_smooth(method = "lm") + theme_bw() + labs(x="Coral Age Distance",y="Bray-Curtis Distance")


# Mantel Test ####
spatial.dist = vegdist(cbind(meta$LON, meta$LAT))
mantel.bray = mantel.rtest(spatial.dist, bray, nrepet = 999)
mantel.jaccard = mantel.rtest(spatial.dist, jaccard, nrepet = 999)

plot(mantel.bray)
plot(mantel.jaccard)


# Multiple Regression on distance matrices ####
dist_MRM <- MRM(bray ~ spatial.dist,  nperm = 9999)
age_MRM <- MRM(bray ~ age,  nperm = 999)

sink("./output/MRM_Table.txt")
print("Bray-Curtis distance regressed against spatial distance:")
print(dist_MRM)
sink(NULL)


# Ordinations ####
ord <- ordinate(ps_ra,method = "PCoA")
plot_ordination(ps_ra,ord,color="CoralAgeBinned")
ggsave("./output/figs/PCoA_CoralAgeGroups.png")

plot_ordination(ps_ra,ord,color="Location")
ggsave("./output/figs/PCoA_Location.png")

growthrate <- cut(meta$Average_LE_mm, breaks = 3)
growthrate <- plyr::mapvalues(growthrate,from = levels(growthrate), to=c("Low","Med","High"))
ps_ra@sam_data$GrowthRateCat <- growthrate

plot_ordination(ps_ra,ord,color="GrowthRateCat") + labs(color="Growth Rate")
ggsave("./output/figs/PCoA_GrowthRate.png")


# Non-metric multidimensional scaling ####
bray.nmds <- monoMDS(bray)
stressplot(bray.nmds)
jaccard.nmds <- monoMDS(jaccard)
stressplot(jaccard.nmds)

# Build data frame
bray.x <- bray.nmds$points[,1]
bray.y <- bray.nmds$points[,2]
jaccard.x <- jaccard.nmds$points[,1]
jaccard.y <- jaccard.nmds$points[,2]

nmds <- data.frame(Bray.X = bray.x,Bray.Y=bray.y,Jaccard.X=jaccard.x,Jaccard.Y=jaccard.y)
nmds.df <- (cbind(meta,nmds))

# plot NMDS results

ggplot(nmds.df, aes(x=Bray.X,y=Bray.Y,color=Location)) +
  geom_point()
ggsave("./output/figs/NMDS_Location.png")

# ADONIS ####
perm.mod <- adonis(otu ~ meta$Location)

sink("./output/PermANOVA_Table.txt")
perm.mod
print("Location is a significant factor in bacterial community structure.")
sink(NULL)

# Network plot ####
ig=make_network(ps_ra, max.dist = .9)
plot_network(ig, physeq = ps_ra, color = "Location",label = NULL,point_size = 2)
ggsave("./Output/Network_Jaccard.png", dpi=300, height = 12, width = 18)



# coral age vs community structure
# plot the unbinned coral age against ESV richness or something similar

# In addition to all the stuff you did before can you do an 
# ordination coloured by binned age, and a separate ordination 
# coloured by the average LE (mm) rate? 
# Essentially what I’m trying to see is whether or not age or growth rate 
# structure bacterial community in addition to the the standard biogeog stuff we have been doing.
# x=age,y=richness,facet=location

# include location / size in model?
# mantel test to rule out? location as factor instead of age

dim(otu_table(ps_ra))
dim(sample_data(ps_ra))
adonis((otu_table(ps_ra)) ~ ps_ra@sam_data$CoralAgeBinned)
