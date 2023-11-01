library(tidyverse)
library(ggplot2)
library(factoextra)

setwd('~/Dropbox (University of Oregon)/skyIslands/')

rm(list=ls())

load('data/spec_RBCL_16s.Rdata')

#meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Sex', 'GeographyFK', 'Site', 'Meadow')

microbes <- spec.net %>%
  select(UniqueID, Site, Apidae, Genus, Species, starts_with('16s')) %>%
  #group_by(Genus) %>%
  filter(Genus == 'Apis' | Genus == 'Bombus' | Genus == 'Megachile' | Genus == 'Anthophora') %>%
  mutate(gensp = paste(Genus, Species)) %>%
  arrange(gensp) %>%
  na.omit()

pollen <- spec.net %>%
  select(UniqueID, Site, Genus, Species, starts_with('RBCL')) %>%
  #group_by(Genus) %>%
  filter(Genus == 'Apis' | Genus == 'Bombus' | Genus == 'Megachile' | Genus == 'Anthophora') %>%
  mutate(gensp = paste(Genus, Species)) %>%
  mutate(gensp = paste(Genus, Species)) %>%
  arrange(gensp) %>%
  na.omit()

# meta <- spec.net %>%
#   mutate(gensp = paste(Genus, Species)) %>%
#   # arrange(gensp) %>%
#   select(all_of(meta_cols), gensp)



rel_abund_assay<- microbes %>% select(starts_with('16s'))

# Calculates Bray-Curtis distances between samples. 
bray_curtis_dist <- vegan::vegdist(rel_abund_assay, method = "bray")

# PCoA 
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)

# All components could be found here: # bray_curtis_pcoa$vectors 
# But we only
#need the first two to demonstrate what we can do: 
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], pcoa2 = bray_curtis_pcoa$vectors[,2])



# Create a plot
bray_curtis_plot1 <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2)) + 
  geom_point(aes(fill=microbes$Genus), alpha=0.7, size=5, pch=21, color='black') +
  scale_fill_manual(values=plasma(4)) +
  stat_ellipse(aes(linetype= microbes$Genus)) +
  scale_linetype_manual("Genus", values=c("longdash", "solid", "dotted", "twodash")) +
  labs(x = paste("PC1 (", round(bray_curtis_pcoa$values[1], 2), "%)"), 
       y = paste("PC2 (", round(bray_curtis_pcoa$values[2], 2), "%)"), fill='Genus') + # makes titles smaller
  theme_classic() +
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.position='bottom', text = element_text(size = 20))

bray_curtis_plot1


bray_curtis_pcoa_df_2 <- data.frame(pcoa3 = bray_curtis_pcoa$vectors[,3],
                                    pcoa4 = bray_curtis_pcoa$vectors[,4])



# Create a plot
bray_curtis_plot2 <- ggplot(data = bray_curtis_pcoa_df_2, aes(x=pcoa3, y=pcoa4)) + 
  geom_point(aes(fill=microbes$Genus), alpha=0.7, size=5, pch=21, color='black') +
  scale_fill_manual(values=plasma(4)) +
  stat_ellipse(aes(linetype= microbes$Genus)) +
  scale_linetype_manual("Genus", values=c("longdash", "solid", "dotted", "twodash")) +
  labs(x = paste("PC3 (", round(bray_curtis_pcoa$values[3], 2), "%)"), 
       y = paste("PC4 (", round(bray_curtis_pcoa$values[4], 2), "%)"), fill='Genus') + # makes titles smaller
  theme_classic() +
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.position='bottom', text = element_text(size = 20))

bray_curtis_plot2


########### do the same for rbcl #######################


rel_abund_assay_poll <- pollen %>% select(starts_with('RBCL'))

# Calculates Bray-Curtis distances between samples. 
bray_curtis_dist_poll <- vegan::vegdist(rel_abund_assay_poll, method = "bray")

# PCoA 
bray_curtis_pcoa_poll <- ecodist::pco(bray_curtis_dist_poll)

# All components could be found here: # bray_curtis_pcoa$vectors 
# But we only
#need the first two to demonstrate what we can do: 
bray_curtis_pcoa_df_poll <- data.frame(pcoa1 = bray_curtis_pcoa_poll$vectors[,1], pcoa2 = bray_curtis_pcoa_poll$vectors[,2])



# Create a plot
bray_curtis_plot1_poll <- ggplot(data = bray_curtis_pcoa_df_poll, aes(x=pcoa1, y=pcoa2)) + 
  geom_point(aes(fill=pollen$Genus), alpha=0.7, size=5, pch=21, color='black') +
  scale_fill_manual(values=plasma(4)) +
  stat_ellipse(aes(linetype= pollen$Genus)) +
  scale_linetype_manual("Genus", values=c("longdash", "solid", "dotted", "twodash")) +
  labs(x = paste("PC1 (", round(bray_curtis_pcoa_poll$values[1], 2), "%)"), 
       y = paste("PC2 (", round(bray_curtis_pcoa_poll$values[2], 2), "%)"), fill='Genus') + # makes titles smaller
  theme_classic() +
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.position='bottom', text = element_text(size = 20))

bray_curtis_plot1_poll


bray_curtis_pcoa_df_2_poll <- data.frame(pcoa3 = bray_curtis_pcoa_poll$vectors[,3],
                                    pcoa4 = bray_curtis_pcoa_poll$vectors[,4])



# Create a plot
bray_curtis_plot2_poll <- ggplot(data = bray_curtis_pcoa_df_2_poll, aes(x=pcoa3, y=pcoa4)) + 
  geom_point(aes(fill=pollen$Genus), alpha=0.7, size=5, pch=21, color='black') +
  scale_fill_manual(values=plasma(4)) +
  stat_ellipse(aes(linetype= pollen$Genus)) +
  scale_linetype_manual("Genus", values=c("longdash", "solid", "dotted", "twodash")) +
  labs(x = paste("PC3 (", round(bray_curtis_pcoa_poll$values[3], 2), "%)"), 
       y = paste("PC4 (", round(bray_curtis_pcoa_poll$values[4], 2), "%)"), fill='Genus') + # makes titles smaller
  theme_classic() +
  guides(fill=guide_legend(nrow=2)) +
  theme(legend.position='bottom', text = element_text(size = 20))

bray_curtis_plot2_poll

