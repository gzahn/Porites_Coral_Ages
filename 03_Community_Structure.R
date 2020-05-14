# -----------------------------------------------------------------------------------#
# Porites 16S Analyses - Beta diversity and community structure
# 
# Depends on "00_Process_Raw_Reads.R" and "01_Cleanup_Processed_Data.R"
#
# Investigate Beta Diversity along geographic, genotype, and age variation of hosts
# 
# Author: Geoffrey Zahn
# -----------------------------------------------------------------------------------#

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
library(rsq)
library(microbiome)
library(RColorBrewer)

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
  geom_smooth(method = "lm") + theme_bw() + labs(x="Coral Age Distance (Years)",y="Bray-Curtis Distance")
ggsave("See Coral_Age_Dist_vs_Community_Dist_Plot.png",dpi=300)


# Quick glm model
mod1 = glm(bray~age)

# write results to file
sink("./output/Age_Distance_vs_Community_Distance_GLM-Table.txt")
print("Bacterial community dissimilarity as a function of coral age distance...")
summary(mod1)
print("Corals of more similar ages have more similar bacterial communities.")
paste0("This is a significant (P = ", signif(coef(summary(mod1))[2,4],4),"), but weak (R-sq = ",signif(rsq(mod1),4),") relationship.")
print("See Coral_Age_Dist_vs_Community_Dist_Plot.png")
sink(NULL)


# Mantel Test ####
spatial.dist = vegdist(cbind(meta$LON, meta$LAT))
mantel.bray = mantel.rtest(spatial.dist, bray, nrepet = 999)
mantel.jaccard = mantel.rtest(spatial.dist, jaccard, nrepet = 999)

sink("./output/Mantel_Test.txt")
mantel.bray
sink(NULL)

plot(mantel.bray)
plot(mantel.jaccard)


# Multiple Regression on distance matrices ####
dist_MRM <- MRM(bray ~ spatial.dist,  nperm = 9999)
age_MRM <- MRM(bray ~ age,  nperm = 999)

sink("./output/MRM_Table.txt")
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices):")
print(dist_MRM)
print("Geographic distance is positively associated with bacterial community distance.")
sink(NULL)


# Ordinations ####

# Subset ps_ra to only samples with defined coral ages...need to re-load fresh RDS later
ps_ra <- subset_samples(ps_ra,!is.na(ps_ra@sam_data$CoralAgeBinned))

ord <- ordinate(ps_ra,method = "NMDS")
stress <- ord$stress
stress_R2 <- 1 - ord$stress^2

png("./output/figs/Stressplot_NMDS_CoralAgeBinned.png")
stressplot(ord)
dev.off()


plot_ordination(ps_ra,ord,color="CoralAgeBinned") +
  scale_color_discrete(limits=c("10 to 20","21 to 30","31 to 40",
                            "41 to 50","51 to 60","61 to 70",
                            "81 to 90","82 to 90","91 to 100",
                            "101 to 110")) + 
  theme_bw() + labs(color="Coral age group",caption = paste0("Stress = ",signif(stress,3),"    Non-metric fit, R^2 = ",signif(stress_R2,3)))
ggsave("./output/figs/NMDS_CoralAgeGroups.png",dpi=300,height = 8,width = 8)

plot_ordination(ps_ra,ord,color="Location") +
  theme_bw() +
  labs(color="Location",caption = paste0("Stress = ",signif(stress,3),"    Non-metric fit, R^2 = ",signif(stress_R2,3)))
ggsave("./output/figs/NMDS_Location.png")


# re-load ps_ra 
ps_ra <- readRDS("./output/phyloseq_cleaned_relabund.RDS")

# subset to remove samples with missing coral age data
aged_samples <- meta$CoralAge != "NA"
aged_samples <- which(aged_samples == TRUE)
otu_aged <- otu[aged_samples,]
meta_aged <- meta[aged_samples,]
 
   
growthrate <- cut(meta$Average_LE_mm, breaks = 3)
growthrate <- plyr::mapvalues(growthrate,from = levels(growthrate), to=c("Low","Med","High"))
ps_ra@sam_data$GrowthRateCat <- growthrate

# subset to remove NA values .. re-load later
ps_ra <- subset_samples(ps_ra,!is.na(ps_ra@sam_data$GrowthRateCat))
ord <- ordinate(ps_ra,method = "NMDS")
stress <- ord$stress
stress_R2 <- 1 - ord$stress^2

png("./output/figs/Stressplot_NMDS_GrowthRate.png")
stressplot(ord)
dev.off()


plot_ordination(ps_ra,ord,color="GrowthRateCat") + labs(color="Growth Rate") + theme_bw() +
  labs(caption = paste0("Stress = ",signif(stress,3),"    Non-metric fit, R^2 = ",signif(stress_R2,3)))
ggsave("./output/figs/NMDS_GrowthRate.png")

# re-load
ps_ra <- readRDS("./output/phyloseq_cleaned_relabund.RDS")


# Non-metric multidimensional scaling ####
bray.nmds <- monoMDS(bray)
stressplot(bray.nmds)
stress.bray <- bray.nmds$stress
stress_R2.bray <- 1 - bray.nmds$stress^2


jaccard.nmds <- monoMDS(jaccard)
stressplot(jaccard.nmds)
stress.jaccard <- jaccard.nmds$stress
stress_R2.jaccard <- 1 - jaccard.nmds$stress^2

# Build data frame
bray.x <- bray.nmds$points[,1]
bray.y <- bray.nmds$points[,2]
jaccard.x <- jaccard.nmds$points[,1]
jaccard.y <- jaccard.nmds$points[,2]

nmds <- data.frame(Bray.X = bray.x,Bray.Y=bray.y,Jaccard.X=jaccard.x,Jaccard.Y=jaccard.y)
nmds.df <- (cbind(meta,nmds))

# plot NMDS results

ggplot(nmds.df, aes(x=Bray.X,y=Bray.Y,color=Location)) +
  geom_point() + theme_bw() +
  labs(title = "NMDS using Bray-Curtis dissimilarity",
       caption = paste0("Stress = ",signif(stress.bray,3),"    Non-metric fit, R^2 = ",signif(stress_R2.bray,3)),
       x="NMDS1",y="NMDS2")
ggsave("./output/figs/NMDS_Location_Bray.png",dpi=300,height = 8,width = 8)


ggplot(nmds.df, aes(x=Jaccard.X,y=Jaccard.Y,color=Location)) +
  geom_point() + theme_bw() +
  labs(title = "NMDS using Jaccard dissimilarity",
       caption = paste0("Stress = ",signif(stress.jaccard,3),"    Non-metric fit, R^2 = ",signif(stress_R2.jaccard,3)),
       x="NMDS1",y="NMDS2")
ggsave("./output/figs/NMDS_Location_Jaccard.png",dpi=300,height = 8,width = 8)

# ADONIS ####

perm.mod1 <- adonis(otu_aged ~ meta_aged$Location * meta_aged$CoralAge)
perm.mod2 <- adonis(otu ~ meta$Location)
sink("./output/PermANOVA_Table.txt")
perm.mod1
print("Location is a significant factor in bacterial community structure. Coral age is not.")
perm.mod2
print("When ALL islands are included in the model, location is still a significant factor.")
sink(NULL)

# Network plot ####
ig=make_network(ps_ra, max.dist = .8)
set.seed(13)
plot_network(ig, physeq = ps_ra, color = "Location",label = NULL,point_size = 2)
ggsave("./output/figs/Network_Jaccard_Location.png", dpi=300, height = 12, width = 18)

set.seed(13)
plot_network(ig, physeq = ps_ra, color = "CoralAgeBinned",label = NULL,point_size = 2)
ggsave("./output/figs/Network_Jaccard_CoralAge.png", dpi=300, height = 12, width = 18)

# 
# # Core microbiome ####
# 
# # find core members of oral microbiome ####
# detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
# 
# # Also define gray color palette
# gray <- gray(seq(0,1,length=10))
# plot_core(ps_ra, plot.type = "heatmap", colours = gray,
#                prevalences = prevalences, detections = detections) +
#   xlab("Detection Threshold (Relative Abundance (%))")
# print(p)    
# 
# 
# # Core heatmap
# prevalences <- seq(.05, .5, .05)
# detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
# 
# p <- plot_core(ps_ra, plot.type = "heatmap", 
#                prevalences = prevalences,
#                detections = detections,
#                colours = rev(brewer.pal(5, "Spectral")),
#                min.prevalence = .1, horizontal = TRUE)
# print(p)
# ggsave("./16S/output/Core_Heatmap.png", dpi=300)
# 
# # prevalence
# taxa_prevalence = (prevalence(ps_ra, detection = 0.001, sort = TRUE))
# 
# # core members at >= 0.2 sample prevalence 0.01% detection threshold
# core.taxa.standard <- core_members(ps_ra, detection = 0.001, prevalence = .1)
# 
# # Total core abundance in each sample (sum of abundances of the core members):
# ps_core <- core(ps_ra, detection = 0.001, prevalence = .1)
# # subset to only samples containing core microbiome
# ps_core_samples = subset_samples(ps_core,rowSums(otu_table(ps_core)) > 0)
# 
# nmds_core = ordinate(ps_core_samples, "NMDS")
# plot_ordination(ps_core_samples, nmds_core, color = "Diet") +
#   stat_ellipse()
# 
# # coral age vs community structure
# # plot the unbinned coral age against ESV richness or something similar
# 
# # In addition to all the stuff you did before can you do an 
# # ordination coloured by binned age, and a separate ordination 
# # coloured by the average LE (mm) rate? 
# # Essentially what I’m trying to see is whether or not age or growth rate 
# # structure bacterial community in addition to the the standard biogeog stuff we have been doing.
# # x=age,y=richness,facet=location
# 
# 
# dim(otu_table(ps_ra))
# dim(sample_data(ps_ra))
# adonis((otu_table(ps_ra)) ~ ps_ra@sam_data$CoralAgeBinned)


# Core microbiome ####

# find core members of oral microbiome ####
# Core with compositionals:
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(ps_ra, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance")
ggsave("./output/figs/CoreSize_vs_RelativeAbundance.png", dpi = 300)

detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=10))

# Core heatmap
prevalences <- seq(.05, .5, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# rename ESVs
colnames(ps_ra@otu_table) <- paste0("ESV_",1:length(colnames(ps_ra@otu_table)))

plot_core(ps_ra, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, horizontal = TRUE)
ggsave("./output/figs/Core_Heatmap.png", dpi=300, height = 8,width = 12)

# prevalence
taxa_prevalence = (prevalence(ps_ra, detection = 0.001, sort = TRUE))

# core members at >= 0.2 sample prevalence 1% detection threshold
core.taxa.standard <- core_members(ps, detection = 0.01, prevalence = .2)
core.taxa.numbers <- as.numeric(unlist(purrr::map(str_split(core.taxa.standard,"_"),2)))

# subset to only samples containing core microbiome
ps_core <- subset_taxa(ps, taxa_names(ps) %in% core.taxa.standard)


# Save core bacterial community taxonomic info
write.csv(ps_core@tax_table, "./output/Core_Microbiome_Taxonomy.csv",quote = FALSE)


# Merge core phyloseq by Location
psm_core <- merge_samples(ps_core,"Location")
psm_core@sam_data$Location <- row.names(psm_core@sam_data)


# convert to relative abundance for merged core microbiome
psm_core_ra <- transform_sample_counts(psm_core, function(x) x / sum(x))
sample_data(psm_core_ra)

# Barplot
plot_bar2(psm_core_ra,fill="Phylum") + 
  theme(axis.text.x = element_text(angle=90,hjust = -.01),
        axis.title = element_text(face="bold",size=12),
        legend.title = element_text(size=12,face="bold"),
        axis.title.x = element_text(vjust = -1.5))
ggsave("./output/figs/Core_Microbiome_Phylum-Level_by_Location.png", dpi=300, width = 10,height = 8)

# *** Try core microbiome with coral age dist vs community dist !!!!!



# coral age vs community structure
# plot the unbinned coral age against ESV richness or something similar

# In addition to all the stuff you did before can you do an 
# ordination coloured by binned age, and a separate ordination 
# coloured by the average LE (mm) rate? 
# Essentially what I’m trying to see is whether or not age or growth rate 
# structure bacterial community in addition to the the standard biogeog stuff we have been doing.
# x=age,y=richness,facet=location


dim(otu_table(ps_ra))
dim(sample_data(ps_ra))
adonis((otu_table(ps_ra)) ~ ps_ra@sam_data$CoralAgeBinned)

