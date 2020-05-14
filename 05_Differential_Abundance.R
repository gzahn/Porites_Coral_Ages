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
library(broom)

source("./R/plot_bbdml2.R")

summary.bbdml <- function(object, ...) {
  # For now, Wald test
  coef.table <- waldt(object)
  # keep <- match(c("call", "df.model", "df.residual", "logL", "link", "phi.link", "formula", "phi.formula", "np.mu", "np.phi", "sep_da", "sep_dv"),
  #               names(object), 0L)
  keep <- match(c("call", "df.model", "df.residual", "logL", "link", "phi.link", "formula", "phi.formula", "np.mu", "np.phi", "sep_da", "sep_dv", "mu.resp", "phi.resp"),
                names(object), 0L)
  ans <- c(object[keep],
           list(coefficients = coef.table))
  
  class(ans) <- "summary.bbdml"
  return(ans)
}

# custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44","#b3bf2d",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000")

# Read cleaned data
ps <- readRDS("./output/phyloseq_object_16S_cleaned.RDS")

# remove Chloroplasts
ps <- subset_taxa(ps, Order != "Chloroplast")

# metadata
meta = as.data.frame(sample_data(ps))
glimpse(meta)

# Differential Abundance Based on Coral Age ####


# Find differentially-abundant taxa
set.seed(123)
da_analysis <- differentialTest(formula = ~ CoralAgeBinned, #abundance
                                phi.formula = ~ CoralAgeBinned, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps,
                                fdr_cutoff = 0.05)

# possibly not enough replication within coral age groups
# some taxa were not present at all in certain coral age groupings

# Try it with merged taxa at class/order/family levels
ps <- clean_taxa_names(ps)

rank_names(ps)
ps_family <- tax_glom(ps,"Family")

set.seed(123)
da_analysis <- differentialTest(formula = ~ CoralAgeBinned, #abundance
                                phi.formula = ~ CoralAgeBinned, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_family,
                                fdr_cutoff = 0.05)

da_analysis$significant_taxa



stax_1 <- paste(tax_table(ps)[da_analysis$significant_taxa,2:5][1],sep="_")
stax_1 <- paste(stax_1[1],stax_1[2],stax_1[3],stax_1[4],sep="_")
stax_2 <- paste(tax_table(ps)[da_analysis$significant_taxa,2:5][2],sep="_")
stax_2 <- paste(stax_2[1],stax_2[2],stax_2[3],stax_2[4],sep="_")
stax_3 <- paste(tax_table(ps)[da_analysis$significant_taxa,2:5][3],sep="_")
stax_3 <- paste(stax_3[1],stax_3[2],stax_3[3],stax_3[4],sep="_")


names(da_analysis$significant_models) <- c(stax_1,stax_2,stax_3)

sink("./output/Differential_abundance_model_stats_tables.txt")
print("Family-level taxonomic comparisons...")
print("No sig. diffs in dispersion for these taxa, but differences in abundance are shown below:")
da_analysis$significant_models
sink(NULL)


daplot_1 <- plot(da_analysis) + labs(y="Differentially abundant families\n(relative abundance)") +
  theme(axis.text.y = element_text(face="bold.italic"),
        axis.title.y = element_text(face="bold",size=16),
        strip.text = element_text(face="bold",size=12))
ggsave(daplot_1,filename="./output/figs/Diff_Abund_Family_by_BinnedAge.png",width = 44,height = 6,device = "png",dpi=400)


#individual taxa deeper analyses ####

set.seed(123)
corncob_da1 <- bbdml(formula = OTU8 ~ CoralAgeBinned,
                     phi.formula = ~ CoralAgeBinned,
                     data = ps_family)


# pull out model results into df
corncob_da1_wald <- waldt(corncob_da1)
corncob_da1_wald <- corncob_da1_wald[grep("mu.",row.names(corncob_da1_wald)),]
corncob_da1_wald <- tidy(corncob_da1_wald)
corncob_da1_wald$OTU <- "Pseudomonadales - Moraxellaceae"

set.seed(123)
corncob_da2 <- bbdml(formula = OTU14 ~ CoralAgeBinned,
                     phi.formula = ~ CoralAgeBinned,
                     data = ps_family)


# pull out model results into df
corncob_da2_wald <- waldt(corncob_da2)
corncob_da2_wald <- corncob_da2_wald[grep("mu.",row.names(corncob_da2_wald)),]
corncob_da2_wald <- tidy(corncob_da2_wald)
corncob_da2_wald$OTU <- "Pseudomonadales - Pseudomonadaceae"

set.seed(123)
corncob_da3 <- bbdml(formula = OTU20 ~ CoralAgeBinned,
                     phi.formula = ~ CoralAgeBinned,
                     data = ps_family)


# pull out model results into df
corncob_da3_wald <- waldt(corncob_da3)
corncob_da3_wald <- corncob_da3_wald[grep("mu.",row.names(corncob_da3_wald)),]
corncob_da3_wald <- tidy(corncob_da3_wald)
corncob_da3_wald$OTU <- "Betaproteobacteriales - Burkholderiaceae"

# join all 3 together
full_corncob <- rbind(corncob_da1_wald,corncob_da2_wald,corncob_da3_wald)
 
# plot
tidy_corncob <- full_corncob %>% select(AgeGroup = .rownames, Estimate, StdErr = Std..Error, t.value, P.val = Pr...t..,OTU) %>%
  filter(AgeGroup != "mu.(Intercept)") %>%
  # arrange(AgeGroup) %>%
  mutate(ymin = Estimate - StdErr, ymax=Estimate + StdErr)
  
tidy_corncob$AgeGroup <-  str_remove(tidy_corncob$AgeGroup,pattern = "mu.CoralAgeBinned") 
  



ggplot(tidy_corncob, aes(x=AgeGroup,y=Estimate)) +
  geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = .2) + theme_bw() +
  geom_hline(yintercept =0,linetype=2,alpha=.5) +
  coord_flip() +
  facet_grid(~OTU) + 
  labs(x="Coral age group (years)", y= "Wald test estimate") + theme(strip.text = element_text(size=14,face="bold"),
                                    axis.title = element_text(size=14,face="bold"),
                                    axis.text = element_text(size=12,face = "bold"))
ggsave("./output/figs/Diff_Abund_Family_by_BinnedAge_Individual_Taxa.png", width = 20,height = 8,dpi=300)
