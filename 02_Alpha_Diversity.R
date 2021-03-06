# -------------------------------------------------------------------#
# Porites 16S Analyses - Depends on "01_Process_Raw_Reads.R"
# and "01_Cleanup_Processed_Data.R"
#
# Investigate Alpha Diversity of bacteria associated with P lutea
#
# Author: Geoffrey Zahn
# -------------------------------------------------------------------#


# Load packages and functions ####
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
library(ggmap)
library(maps)
library(rsq)
library(ggpmisc)
library(viridis)

# functions
source("./R/plot_bar2.R")
source("./R/summarize_taxa_Joey711.R")
source("./R/heatmap_left.R")

# custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000")

# Read cleaned data
ps_ra <- readRDS("./output/phyloseq_cleaned_relabund.RDS")

# remove Chloroplast reads
ps_ra <- subset_taxa(ps_ra,Order != "Chloroplast")

# convert OTU table and metadata to data.frame for easier downstream access
otu = as.data.frame(as(otu_table(ps_ra),"matrix"))
meta = as.data.frame(sample_data(ps_ra))

# convert OTU table to data.frame for easier downstream access
otu = as.data.frame(as(otu_table(ps_ra),"matrix"))

# Look for covariance and correlation between coral age and location ####

# just numeric columns
meta.numeric <- meta[,names(meta) %in% c("CoralAge","LAT","LON")]
meta.numeric <- meta.numeric[complete.cases(meta.numeric),]

cov(meta.numeric)
corrplot(cor(meta.numeric), method = "color",title = "Correlation Between Coral Age and Location")


# Calculate alpha diversity measures ####
shannon = diversity((otu),index = "shannon")
richness = specnumber(otu)

# Does genotype affect bacterial alpha diversity? ####
mod.gen <- aov(shannon$shannon ~ meta$Genotype)
summary(mod.gen)

# Write t-test results to file

sink("./output/alpha_div_vs_genotype_t-test_tables.txt")
print("Shannon Diversity")
t.test(shannon$shannon ~ meta$Genotype)
print("Species Richness")
t.test(richness ~ meta$Genotype)
sink(NULL)


# Does bacterial richness change with coral age? ####
aged = !is.na(meta$CoralAge)

mean.ages <- meta %>% group_by(Location) %>% drop_na(CoralAge) %>%
  summarise(MeanCoralAge=mean(CoralAge))


mod = aov(richness[aged] ~ CoralAge * Location, data = meta[aged,])
summary(mod) # probable interaction between coral age and location

sink("./output/richness_vs_CoralAge_and_Location_ANOVA.txt")
summary(mod)
sink(NULL)
                    # Richness
p1 = ggplot(meta[aged,],aes(x=CoralAge,y=richness[aged])) + 
  geom_smooth(method = "lm",se=FALSE,color="Red",alpha=.5,linetype=2) +
  geom_point() + 
  stat_fit_glance(method = "lm", 
                  label.x = "right",
                  label.y = "top",
                  aes(label = sprintf('R^2~"="~%.3f~"\n"~italic(P)~"="~%.2f',
                                      stat(r.squared), stat(p.value))),
                  parse = TRUE) +
  facet_grid(~Location) + theme_bw() +
  labs(x="Coral Age (years)",y="ESV Richness") + theme(axis.title = element_text(face="bold",size = 14),
                                                       strip.text = element_text(face="bold",size=10))

  ggsave(p1,filename = "./output/figs/RichnessLinear_Fit.png",dpi=300,device = "png")

  

# explore models
meta.div <- cbind(meta,shannon,richness)
mod1 = glm(data = meta.div, richness ~ CoralAge*Location)
mod2 = glm(data = meta.div, richness ~ CoralAge+Location)
mod3 = glm(data = meta.div, richness ~ ns(CoralAge,2)*Location)
# mod4 = lmer(data = meta.div, richness ~ CoralAge * (1|Location))

grid <- meta.div[aged,] %>% 
  data_grid(CoralAge,Location) %>% 
  gather_predictions(mod1, mod2, mod3)

ggplot(meta.div[aged,], aes(x=CoralAge, y=richness,group=Location,color=Location)) + 
  geom_point() +
  geom_line(data = grid,aes(y=pred)) +
  facet_wrap(~ model) + theme_bw() + scale_color_manual(values=pal)
ggsave("./output/figs/RichnessModel_Comparison.png",dpi=300)

# Plot residuals
gather_residuals(meta.div[aged,],mod1,mod2,mod3) %>%
  ggplot(aes(x=CoralAge,y=resid,color=Location)) + geom_point() + geom_ref_line(h=0,size = .5,colour = "Black") +
  facet_grid(model~Location) + theme_bw() + scale_color_manual(values=pal)
ggsave("./output/figs/RichnessModel_Residuals.png",dpi=300)

                    # Shannon
p2 <- ggplot(meta[aged,],aes(x=CoralAge,y=shannon[aged])) + 
  geom_smooth(method = "lm",se=FALSE,color="DodgerBlue",alpha=.5,linetype=2) +
  geom_point() + 
  stat_fit_glance(method = "lm", 
                  label.x = "right",
                  label.y = "top",
                  aes(label = sprintf('R^2~"="~%.3f~"\n"~italic(P)~"="~%.2f',
                                      stat(r.squared), stat(p.value))),
                  parse = TRUE) +
  facet_grid(~Location) + theme_bw() +
  labs(x="Coral Age (Years)",y="Shannon Diversity") + theme(axis.title = element_text(face="bold",size = 14),
                                                            strip.text = element_text(face="bold",size=10))
ggsave(p2,filename = "./output/figs/ShannonLinear_Fit.png",dpi=300,device = "png")

# explore models
mod1 = glm(data = meta.div, shannon ~ CoralAge*Location)
mod2 = glm(data = meta.div, shannon ~ CoralAge+Location)
mod3 = glm(data = meta.div, shannon ~ ns(CoralAge,2)*Location)

# mod4 = lmer(data = meta.div, shannon ~ CoralAge * (1|Location))

grid <- meta.div[aged,] %>% 
  data_grid(CoralAge,Location) %>% 
  gather_predictions(mod1, mod2, mod3)

ggplot(meta.div[aged,], aes(x=CoralAge, y=shannon,group=Location,color=Location)) + 
  geom_point() +
  geom_line(data = grid,aes(y=pred)) +
  facet_wrap(~ model) + theme_bw() + scale_color_manual(values=pal) + labs(y="Shannon diversity",x="Coral age")
ggsave("./output/figs/ShannonModel_Comparison.png",dpi=300)

# Plot residuals
gather_residuals(meta.div[aged,],mod1,mod2,mod3) %>%
  ggplot(aes(x=CoralAge,y=resid,color=Location)) + geom_point() + geom_ref_line(h=0,size = .5,colour = "Black") +
  facet_grid(model~Location) + theme_bw() + scale_color_manual(values=pal)
ggsave("./output/figs/ShannonModel_Residuals.png",dpi=300)


                    # Simpson
p3 <- ggplot(meta[aged,],aes(x=CoralAge,y=simpson[aged])) + 
  geom_smooth(method = "lm",se=FALSE,color="Purple",alpha=.5,linetype=2) +
  geom_point() + 
  facet_grid(~Location) + theme_bw() +
  labs(x="Coral Age (years)",y="Simpson Diversity") + theme(axis.title = element_text(face="bold",size = 14))
ggsave(p3,filename = "./output/figs/SimpsonLinear_Fit.png",dpi=300,device = "png")

# explore models
mod1 = glm(data = meta.div, simpson ~ CoralAge*Location)
mod2 = glm(data = meta.div, simpson ~ CoralAge+Location)
mod3 = glm(data = meta.div, simpson ~ ns(CoralAge,2)*Location)
# mod4 = lmer(data = meta.div, simpson ~ CoralAge * (1|Location))

grid <- meta.div[aged,] %>% 
  data_grid(CoralAge,Location) %>% 
  gather_predictions(mod1, mod2, mod3)

ggplot(meta.div[aged,], aes(x=CoralAge, y=simpson,group=Location,color=Location)) + 
  geom_point() +
  geom_line(data = grid,aes(y=pred)) +
  facet_wrap(~ model) + theme_bw() + scale_color_manual(values=pal)
ggsave("./output/figs/SimpsonModel_Comparison.png",dpi=300)

# Plot residuals
gather_residuals(meta.div[aged,],mod1,mod2,mod3) %>%
  ggplot(aes(x=CoralAge,y=resid,color=Location)) + geom_point() + geom_ref_line(h=0,size = .5,colour = "Black") +
  facet_grid(model~Location) + theme_bw() + scale_color_manual(values=pal)
ggsave("./output/figs/SimpsonModel_Residuals.png",dpi=300)


# Best Model ####

# Merge all linear model fit plots for diversity measures
p3 = p1+labs(x=NULL) + ggtitle("ESV Alpha Diversity") + 
  theme(title = element_text(size=16,face="bold",hjust = .5)) +
  p2+
  plot_layout(ncol=1) 

ggsave(p3,filename = "./output/figs/AlphaDiversity_over_CoralAge.png",dpi=600,height = 8,width = 14,device = "png")


# Best models:
mod <- glm(richness ~ CoralAge*Location,data = meta.div)
summary(mod)
mod2 <- lmer(richness ~ CoralAge/Location + (1|Location), data = meta.div)
summary(mod2)
pred <- add_predictions(meta.div[aged,],mod2)


# lmer model residual plot
ggplot(meta.div[aged,],aes(x=CoralAge,y=residuals(mod2))) + geom_point() + geom_hline(yintercept = 0)

ggplot(meta.div[aged,], aes(x=CoralAge,y=richness,color=Location,group=Location)) +
  geom_point(alpha=.25) + geom_line(aes(y=pred$pred),size=1) + theme_bw() +
  scale_color_manual(values=pal) + ggtitle("Mixed-Effect Model Predictions") + 
  labs(x="Coral Age (Years)",y="ESV Richness",subtitle = "Points are observed values")
ggsave("./output/figs/lmer_model_predictions_richness_v_age.png",dpi=300, width = 8,height = 10)

# Barplots of diversity ####
# Load raw data
ps <- readRDS("./output/phyloseq_object_16S_cleaned.RDS")
ps <- subset_taxa(ps, Order != "Chloroplast")
# fix metadata typo
ps@sam_data$CoralAgeBinned[ps@sam_data$CoralAgeBinned == "82 to 90"] <- "81 to 90"

summary(colSums(otu_table(ps)))

# Drop empty and low-abundance taxa
keepers <- colSums(otu_table(ps)) >= 1000
otu_table(ps) <- otu_table(ps)[,keepers]

# Merge samples by location
ps.Location <- merge_samples(ps,group = "Location")
# Repair values
ps.Location@sam_data$Location <- unique(meta$Location)

# Merge by Age Class
# remove missing coral ages
ps.ages <- subset_samples(ps,!is.na(CoralAgeBinned))

ps.ages@sam_data$CoralAgeBinned <- factor(ps.ages@sam_data$CoralAgeBinned,
                                          levels = c("10 to 20","21 to 30","31 to 40",
                                                     "41 to 50","51 to 60","61 to 70",
                                                     "81 to 90","91 to 100",
                                                     "101 to 110"))


ps.Age <- merge_samples(ps,group="CoralAgeBinned")
ps.Age@sam_data$CoralAgeBinned <- levels(ps.ages@sam_data$CoralAgeBinned)
# Convert to relabund
ps.Location.ra <- transform_sample_counts(ps.Location,fun = function(x) x/sum(x))
ps.Age.ra <- transform_sample_counts(ps.Age,fun = function(x) x/sum(x))
sample_data(ps.Age.ra)


# Change NA to "Unassigned"
tax_table(ps.Location.ra)[,2][is.na(tax_table(ps.Location.ra)[,2])] <- "Unassigned"
tax_table(ps.Age.ra)[,2][is.na(tax_table(ps.Age.ra)[,2])] <- "Unassigned"


plot_bar2(ps.Location.ra,fill = "Phylum") + theme_bw() +
  scale_fill_manual(values=pal) +
  labs(x="Relative Abundance",y="Site")
ggsave("./output/figs/Barplot_Phylum_Location.png",height = 10,width = 12,dpi=300)

plot_bar2(ps.Age.ra,fill = "Phylum") + theme_bw() +
  scale_fill_manual(values=pal) +
  scale_x_discrete(limits=c("10 to 20","21 to 30","31 to 40",
                            "41 to 50","51 to 60","61 to 70",
                            "81 to 90","91 to 100",
                            "101 to 110")) + 
  labs(x="Age Class",y="Relative Abundance")
ggsave("./output/figs/Barplot_Phylum_AgeClass.png",height = 10,width = 12,dpi=300)


plot_bar2(ps.Age.ra,fill = "Class") + theme_bw() +
  scale_fill_viridis_d() + 
  scale_x_discrete(limits=c("10 to 20","21 to 30","31 to 40",
                            "41 to 50","51 to 60","61 to 70",
                            "81 to 90","91 to 100",
                            "101 to 110")) + 
  labs(x="Age Class",y="Relative Abundance")
ggsave("./output/figs/Barplot_Class_AgeClass.png",height = 10,width = 12,dpi=300)

plot_bar2(ps.Age.ra,fill = "Order") + theme_bw() +
  scale_fill_viridis_d() + 
  scale_x_discrete(limits=c("10 to 20","21 to 30","31 to 40",
                            "41 to 50","51 to 60","61 to 70",
                            "81 to 90","91 to 100",
                            "101 to 110")) + 
  labs(x="Age Class",y="Relative Abundance")

ggsave("./output/figs/Barplot_Order_AgeClass.png",height = 10,width = 12,dpi=300)

plot_bar2(ps.Age.ra,fill = "Family") + theme_bw() +
  scale_fill_viridis_d() + 
  scale_x_discrete(limits=c("10 to 20","21 to 30","31 to 40",
                            "41 to 50","51 to 60","61 to 70",
                            "81 to 90","91 to 100",
                            "101 to 110")) + 
  labs(x="Age Class",y="Relative Abundance")

ggsave("./output/figs/Barplot_Family_AgeClass.png",height = 10,width = 12,dpi=300)


# Boxplots of alpha diversity
meta.boxplot <- meta.div[meta.div$Location %in% c("Jong","Kusu","Raffles Lighthouse","Semakau","Sisters"),]

ggplot(meta.boxplot, aes(x=CoralAgeBinned,y=shannon,fill=CoralAgeBinned)) + 
  geom_boxplot() + theme_bw() + labs(fill="Coral age group",x="Coral age group",y="Shannon diversity") +
  facet_wrap(~Location) + theme(axis.text.x = element_text(angle=90),
                                legend.title = element_text(size=12,face="bold"),
                                axis.title = element_text(size=12,face="bold"),strip.text = element_text(face="bold")) +
  scale_x_discrete(limits=c("10 to 20","21 to 30","31 to 40",
                            "41 to 50","51 to 60","61 to 70",
                            "81 to 90","91 to 100",
                            "101 to 110"))
ggsave("./output/figs/Shannon_Diversity_Boxplot_by_Age_and_Location.png", dpi=300, height = 8,width = 10)


ggplot(meta.div, aes(x=Location,y=shannon, fill=Location)) +
  geom_boxplot() + theme_bw() + labs(fill="Location",x="Location",y="Shannon diversity") +
  theme(axis.title = element_text(face="bold",size=12), legend.title = element_text(face="bold",size=12),
        axis.text = element_text(size = 8,face="bold"))
ggsave("./output/figs/Shannon_Diversity_by_Island.png",dpi=300,height = 6,width = 8)

ggplot(meta.boxplot, aes(x=Location,y=richness, fill=Location)) +
  geom_boxplot() + theme_bw() + labs(fill="Location",x="Location",y="ESV Richness") +
  theme(axis.title = element_text(face="bold",size=12), legend.title = element_text(face="bold",size=12),
        axis.text = element_text(size = 8,face="bold"))
ggsave("./output/figs/ESV_Richness_by_Island.png",dpi=300,height = 6,width = 8)


ggplot(meta.boxplot, aes(x=CoralAgeBinned,y=richness,color=CoralAgeBinned)) + 
  geom_boxplot() + theme_bw() + labs(color="Coral age group",x="Coral age group",y="Shannon diversity") +
  facet_wrap(~Location) + theme(axis.text.x = element_text(angle=90),
                                legend.title = element_text(size=12,face="bold"),
                                axis.title = element_text(size=12,face="bold"),strip.text = element_text(face="bold")) +
  scale_x_discrete(limits=c("10 to 20","21 to 30","31 to 40",
                            "41 to 50","51 to 60","61 to 70",
                            "81 to 90","82 to 90","91 to 100",
                            "101 to 110"))
ggsave("./output/figs/Richness_Boxplot_by_Age_and_Location.png", dpi=300, height = 8,width = 10)





# Maps of site locations ####
register_google(key = "AIzaSyDjTK7ZjYKGr1nE5PmaSncAg4g9L913C_o")
maplocation <- c(lon = mean(ps@sam_data$LON), lat = mean(ps@sam_data$LAT))
mygooglemap <- get_googlemap(center = maplocation,
              zoom = 11, scale = 2,
              maptype ='terrain')

mystamenmap <- get_map(location = maplocation,source="stamen",
                       maptype = "watercolor",zoom = 11,color = "bw", crop = FALSE)

mytonermap2 <- get_map(location = maplocation,source="stamen",maptype = "toner",zoom = 11)

ggmap(mygooglemap) + 
  geom_point(aes(x = LON, y = LAT, color=Location), data = meta, size = 4) +
  scale_color_manual(values=pal) + labs(x="Longitude",y="Latitude")
ggsave("./output/figs/googlemap.png", dpi=300, height = 8,width = 10)

ggmap(mystamenmap) + 
  geom_point(aes(x = LON, y = LAT, color=Location), data = meta, size = 4) +
  scale_color_manual(values=pal)+ labs(x="Longitude",y="Latitude")
ggsave("./output/figs/stamenmap.png",dpi=300, height = 8,width = 10)

ggmap(mytonermap2) + 
  geom_point(aes(x = LON, y = LAT,color=Location), data = meta, size = 4) +
  scale_color_manual(values=pal) + labs(x="Longitude",y="Latitude")
ggsave("./output/figs/tonermap.png",dpi=300,height = 8, width = 10)


# Number of ESVs per colony


meta$ReadCount <- rowSums(otu_table(ps))
meta$Richness <- richness
names(meta)
colkeepers <- c("Library ID","Location","CoralAge","GPS","Genotype","LAT","LON","Richness","ReadCount")
meta_concise <- meta %>% select(colkeepers)
write.csv(meta_concise,"./output/SampleData_with_ESVRichness.csv",quote = FALSE,row.names = FALSE)

meta_concise_by_site <- meta_concise %>% group_by(Location) %>%
  summarize(N=n(),Mean_ESV_Richness = mean(Richness))
write.csv(meta_concise_by_site,"./output/Mean_Site_ESV_Richness.csv",quote = FALSE,row.names = FALSE)

