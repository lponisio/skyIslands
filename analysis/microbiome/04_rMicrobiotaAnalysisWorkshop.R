


##following 
# https://rpubs.com/dillmcfarlan/R_microbiotaSOP

library(phyloseq)

library(ggpubr)

#Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq
library(ape)

#This package will also help us more easily manipulate our data
library(dplyr)

#Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.
library(ggplot2)

#This package is used to calculate and plot Venn diagrams as well as heatmaps
library(gplots)

#Linear mixed-effects models like repeated measures analysis
library(lme4)

library(phangorn)
library(picante)

#The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.
library(phyloseq)

#A package to create interactive web graphics of use in 3D plots
library(plotly)

#This package will help us more easily manipulate our data, which are matrices
library(tidyr)

#The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.
library(vegan)

#Pretty Venn disgrams
library(VennDiagram)


library(tidyverse)
library(ggplot2)

setwd('~/Dropbox (University of Oregon)/PonisioLocal/skyIslands/')

spec16s <- read.csv('spec_RBCL_16s.csv') %>%
  tibble::column_to_rownames('UniqueID') %>%
  filter(Apidae == 1) %>%
  filter(!Genus == 'Agapostemon')
  

meta_cols <- c('Family', 'Genus', 'Species', 'Sex', 'GeographyFK', 'Site', 'Meadow')

microbes <- spec16s %>%
  select(Site, Genus, Species, starts_with('X16s')) %>%
  #group_by(Genus) %>%
  filter(!Genus == 'Agapostemon') %>%
  mutate(gensp = paste(Genus, Species)) %>%
  arrange(gensp) %>%
  na.omit()

species_features <- microbes %>%
  select(starts_with('X16s')) 

meta <- spec16s %>%
  mutate(gensp = paste(Genus, Species)) %>%
  # arrange(gensp) %>%
  select(all_of(meta_cols), gensp) 


species.16s <- as.data.frame(species.16s)

tax.clean <- separate(as.data.frame(species.16s), V1, into = c("Phylum", "Class", "Order", "Family", "Genus", "Species", 'Strain'), sep=";")


######
#calculating PD

otu_table_ps1 <- as.data.frame(physeq16sR0@otu_table) 

metadata_table_ps1  <- as.data.frame(physeq16sR0@sam_data)

df.pd <- picante::pd(t(otu_table_ps1), physeq16sR0@phy_tree,include.root=F) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
df.pd <- df.pd[row.names(df.pd) %in% row.names(meta),]
meta <- meta[row.names(meta) %in% row.names(df.pd),]

# now we need to plot PD
# check above how to get the metadata file from phyloseq object.
# We will add the results of PD to this file and then plot.
meta$PhyloDiversity <- df.pd$PD 

meta$SpeciesRichness <- df.pd$SR

library(rstatix)
library(ggarrange)

stat.test.pd <- meta %>% 
  pairwise_wilcox_test(PhyloDiversity ~ Genus, p.adjust.method="fdr") %>%
  filter(p.adj.signif != 'ns')


plot.pd.genus <- ggplot(meta, aes(Genus, PhyloDiversity)) + 
  geom_violin(aes(fill = Genus)) + 
  geom_jitter(size = 2, alpha=0.2) + 
  #coord_flip() + 
  theme_classic() +
  scale_fill_manual(values=plasma(4)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
               colour = "black", width = 1)  +
  stat_summary(fun = "median", geom = "crossbar",
               colour = "red", width = 1) +
  theme(legend.position='none', text=element_text(size=21)) +
  labs(y='Phylogenetic Diversity (Branch Length)') + 
  stat_pvalue_manual(stat.test.pd, label = "p.adj.signif", tip.length = 0.01, y.position = 10, step.increase = 0.05)

print(plot.pd.genus)

stat.test.sr <- meta %>% 
  pairwise_wilcox_test(SpeciesRichness ~ Genus, p.adjust.method="fdr") %>%
  filter(p.adj.signif != 'ns')


plot.sr.genus <- ggplot(meta, aes(Genus, SpeciesRichness)) + 
  geom_violin(aes(fill = Genus)) + 
  geom_jitter(size = 2, alpha=0.2) + 
  #coord_flip() + 
  theme_classic() +
  scale_fill_manual(values=plasma(4)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
               colour = "black", width = 1)  +
  stat_summary(fun = "median", geom = "crossbar",
               colour = "red", width = 1) +
  theme(legend.position='none', text=element_text(size=21)) +
  labs(y='OTU Richness') + 
  stat_pvalue_manual(stat.test.sr, label = "p.adj.signif", tip.length = 0.01, y.position = 100, step.increase = 0.05)
print(plot.sr.genus)

ggarrange(plot.sr.genus, plot.pd.genus)



#lets look at distrubutions -- neither are normal

hist(meta$PhyloDiversity, main="PhyloDiverity", xlab="", breaks=10)
shapiro.test(meta$PhyloDiverity)

#transform to make distribution normal
meta$TransPhyloDiversity <- log(meta$PhyloDiverity+1)
hist(meta$TransPhyloDiversity, main="TransPhyloDiverity", xlab="", breaks=10)


x=log(x+1)

hist(meta$SpeciesRichness, main="SpeciesRichness", xlab="", breaks=10)
meta$TransSpeciesRichness <- log(meta$SpeciesRichness)
hist(meta$TransSpeciesRichness, main="TransPhyloDiverity", xlab="", breaks=10)
library(ggpubr)
ggdensity(meta$SpeciesRichness)
ggdensity(meta$PhyloDiversity)
ggdensity(meta$TransPhyloDiversity)
shapiro.test(meta$TransSpeciesRichness)

#non parametric test for if genus microbes are different
kruskal.test(PhyloDiversity ~ Genus, data=meta)
pairwise.wilcox.test(meta$PhyloDiversity, meta$Genus, p.adjust.method="fdr")

kruskal.test(SpeciesRichness ~ Genus, data=meta)
pairwise.wilcox.test(meta$SpeciesRichness, meta$Genus, p.adjust.method="fdr")

#lets see interaction between genus and site
aov.shannon.all = aov(SpeciesRichness ~ Site*Genus, data=meta)
summary(aov.shannon.all)
glm.shannon.all = glm(SpeciesRichness ~ Site*Genus, data=meta)
summary(glm.shannon.all)

#Create 2x2 plot environment so that we can see all 4 metrics at once. 
par(mfrow = c(2, 2))

#Then plot each metric.
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(meta$simpson, main="Simpson diversity", xlab="", breaks=10)
hist(meta$chao, main="Chao richness", xlab="", breaks=15)
hist(meta$ace, main="ACE richness", xlab="", breaks=15)




BC.nmds = metaMDS(species_features, distance="bray", k=2, trymax=1000)
BC.dist=vegdist(species_features, distance="bray")

#anosim
anosim(BC.dist, meta$Genus, permutations = 1000)

glm.shannon.all = glm(SpeciesRichness ~ Genus, data=meta)
summary(glm.shannon.all)












