rm(list=ls())
library(tidyverse)
library(patchwork)
source("lab_paths.R")
local.path
dir.bombus <- file.path(local.path, "skyIslands/analysis/parasites")
setwd(dir.bombus)
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/standardize_weights.R")
vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "Net_BeeDiversity",
                 "Lat", "SRDoy"  
)
vars_yearsrsp <- "rare.degree"
vars_sp <- "MeanITD"


variables.to.log <- c("rare.degree", "MeanITD")

variables.to.log.1<- c("Net_HBAbundance", "Net_BombusAbundance", 
                       "Net_NonBombusHBAbundance")
source("src/init.R")

spec.orig <- filter(spec.orig, Site != "VC" & Site != "UK" & Site != "SS")

## Plots by meadow
## Bee abundances by meadows
###############################################################################
##Bombus
bombus_abundance<- spec.orig %>% 
ggplot(aes(reorder(Site, Lat, decreasing = TRUE), Net_BombusAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Sites", y = "Bombus Abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))+
  theme_bw()


##Apis
HB_abundance<- spec.orig %>% 
  ggplot(aes(reorder(Site, Lat, decreasing = TRUE), Net_HBAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Sites", y = "Apis Abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))+
  theme_bw()


##Other bees
melissodes_abundance<- spec.orig %>% 
  filter(Genus == "Melissodes") %>% 
  ggplot(aes(reorder(Site, Lat, decreasing = TRUE), Net_NonBombusHBAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Sites", y = "Melissodes Abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))+
  theme_bw()

################################################################################
## Bee Diversity by meadows
bee_diversity<- spec.orig %>% 
  ggplot(aes(reorder(Site, Lat, decreasing = TRUE), Net_BeeDiversity))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Sites", y = "Bee Diversity", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))+
  theme_bw()


###############################################################################
## Floral abundance by meadows
floral_abundance<- spec.orig %>% 
  ggplot(aes(reorder(Site, Lat, decreasing = TRUE), MeanFloralAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Sites", y = "Mean Floral Abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))+
  theme_bw()


## Floral diversity by meadows

floral_diversity<- spec.orig %>% 
  ggplot(aes(reorder(Site, Lat, decreasing = TRUE), MeanFloralDiversity))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Year", y = "Mean Floral Diversity", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))+
  theme_bw()

bombus_abundance + HB_abundance + melissodes_abundance + bee_diversity + 
  floral_abundance + floral_diversity + plot_layout(ncol = 3)+ 
  plot_annotation(tag_levels = "A")


###############################################################################
## Parasite prevalence by meadows

sick.totals <- spec.net %>%
  
  
## Subset to Weights == 1
  group_by(Site) %>%
  summarise(TestedTotals = length(UniqueID[!is.na(ParasitePresence)]),
            ParasitismRate=round(mean(ParasitePresence, na.rm=TRUE),2),
            InfectedIndividuals=round(sum(ParasitePresence, na.rm=TRUE),2),
            InfectedApicystisSpp=round(mean(ApicystisSpp, na.rm=TRUE),2),
            InfectedAscosphaeraSpp=round(mean(AscosphaeraSpp, na.rm=TRUE),2),
            InfectedCrithidiaBombi=round(mean(CrithidiaBombi, na.rm=TRUE),2),
            InfectedCrithidiaExpoeki=round(mean(CrithidiaExpoeki,
                                                na.rm=TRUE), 2),
            InfectedCrithidia=round(mean(CrithidiaPresence,
                                         na.rm=TRUE), 2),
            InfectedNosemaBombi=round(mean(NosemaBombi,
                                           na.rm=TRUE), 2),
            InfectedNosemaCeranae=round(mean(NosemaCeranae,
                                             na.rm=TRUE), 2))

ggplot(sick.totals, aes(x=Site, y=ParasitismRate)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip()

################################################################################
## Relationship between bee abundance 
ggplot(spec.all, aes(x=Net_BombusAbundance,
                     y=Net_HBAbundance, color=Site))+
  geom_point()

## Relationship between bombus abundance and other bees
ggplot(spec.all, aes(x=Net_BombusAbundance,
                     y=Net_NonBombusHBAbundance, color=Site))+
  geom_point()

## Relationship between bombus abundance and latitude
ggplot(spec.all, aes(x= Lat,
                     y=Net_BombusAbundance, color=Site))+
  geom_point()
## Relationship between Hb abundance and latitude
ggplot(spec.net, aes(x= Site,
                     y=Net_HBAbundance))+
  geom_boxplot()

ggplot(spec.net, aes(x= Year,
                     y=Net_HBAbundance))+
  geom_boxplot()

## Relationship between other bees abundance and latitude
ggplot(spec.all, aes(x= Lat,
                     y=Net_NonBombusHBAbundance, color=Site))+
  geom_point()

## Relationship between pollinator abundance and year
ggplot(spec.all, aes(x= Year,
                     y=Net_PollAbundance, color=Site))+
  geom_point()


################################################################################
## Relationships with bee diversity
## Relationship between bee diversity and year
ggplot(spec.net, aes(x= Lat,
                     y=Net_BeeDiversity, color=Site))+
  geom_point()

ggplot(spec.all, aes(Year, Net_BeeDiversity)) + geom_boxplot()

## Relationship between bee diversity and lat
ggplot(spec.all, aes(x= Lat,
                     y=Net_BeeDiversity, color=Site))+
  geom_point()

## Relationship between bee diversity and flower diversity
ggplot(spec.all, aes(x= MeanFloralDiversity,
                     y=Net_BeeDiversity, color=Site))+
  geom_point()

## Relationship between bee diversity and flower abundance 
ggplot(spec.all, aes(x= MeanFloralAbundance,
                     y=Net_BeeDiversity, color=Site))+
  geom_point()

## Boxplot of bumble bee species and their parasite rate

Bombus_ParasiteRate<- spec.net %>% 
  filter(Genus == "Bombus") %>% 
  group_by(GenusSpecies, Site) %>% 
  summarise(ParasitismRate = mean(ParasitePresence, na.rm = TRUE)) %>% 
  ggplot(aes(GenusSpecies, ParasitismRate))+
  geom_boxplot()+
  labs(x = "Bombus Species", y = "Parasitism Rate")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))

ggsave(Bombus_ParasiteRate, file="figures/Bombus_ParasiteRate.pdf",
       height=4, width=5)


## Parasite prevalence in different groups of bees. 

parasites_bees<-pivot_longer(spec.net, cols = c(ApicystisSpp, AscosphaeraSpp, CrithidiaBombi, CrithidiaExpoeki, 
                        CrithidiaPresence, NosemaCeranae, NosemaBombi), names_to = "Parasites") 
parasite_prevalence <- parasites_bees %>%   
  filter(Genus == c("Melissodes", "Anthophora", "Apis", "Bombus")) %>% 
  filter(!is.na(value)) %>% 
  group_by(Site, Genus, Parasites) %>% 
  summarise(ParasitismRate = mean(value, na.rm = TRUE)) %>% 
  ggplot(aes(Parasites, ParasitismRate))+
  geom_boxplot()+
  labs(x = "Parasites", y = "Parasitism Rate")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))+
  facet_wrap(~ Genus)

ggsave(parasite_prevalence, file="figures/parasite_prevalence_by_genus.jpg",
       height=4, width=5)

## Bar plot of parasite prevalence per sites for each genus
parasite_prevalence_sites <- parasites_bees %>%   
  filter(Genus == c("Bombus", "Melissodes", "Apis", "Anthophora")) %>% 
  filter(Parasites != "NosemaBombi" & Parasites != "NosemaCeranae") %>% 
  filter(Site != "VC" & Site != "UK" & Site != "SS")%>% 
  filter(!is.na(value)) %>% 
  group_by(Site, Genus, Parasites) %>% 
  summarise(ParasitismRate = mean(value, na.rm = TRUE)) %>%  
  ggplot(aes(x=Site, y=ParasitismRate, fill = Genus)) +
  geom_bar(stat="identity", position = "dodge") + theme_minimal() + coord_flip()+
  scale_fill_brewer(palette = "Oranges")

ggsave(parasite_prevalence_sites, file="figures/prevalence_by_sites.jpg",
       height=4, width=5)
## Looking at parasitism of each of the parasites per each genus

parasite_prevalence_genus <- parasites_bees %>%   
  filter(Genus == c("Bombus", "Melissodes", "Apis", "Anthophora")) %>% 
  filter(Parasites != "NosemaBombi" & Parasites != "NosemaCeranae") %>% 
  filter(Site != "VC" & Site != "UK" & Site != "SS")%>% 
  filter(!is.na(value)) %>% 
  group_by(Site, Genus, Parasites) %>% 
  summarise(ParasitismRate = mean(value, na.rm = TRUE)) %>%  
  ggplot(aes(x=Parasites, y=ParasitismRate, fill = Genus)) +
  geom_bar(stat="identity", position = "dodge") + theme_minimal() + coord_flip()+
  scale_fill_brewer(palette = "Oranges")

ggsave(parasite_prevalence_genus, file="figures/prevalence_by_genus.jpg",
       height=4, width=5)

## Looking at total screened bees per site and num of positives
parasites_bees$value <- as.factor(parasites_bees$value)
tested_pos_neg<- parasites_bees %>%   
  filter(Genus == c("Bombus", "Melissodes", "Apis", "Anthophora")) %>% 
  filter(Parasites != "NosemaBombi" & Parasites != "NosemaCeranae") %>% 
  filter(Site != "VC" & Site != "UK" & Site != "SS")%>% 
  filter(!is.na(value)) %>% 
  group_by(Site, value) %>% 
  summarise(TestedTotals = length(value), 
            Positives = length(which(value == 1))) %>% 
  ggplot(aes(x=Site)) +
  geom_bar(aes(y = TestedTotals, fill = value), stat ="identity")+
  theme_minimal() + coord_flip()+
  theme(legend.title=element_blank())+
  scale_fill_manual(col = c("grey", "#feb24c"), labels = c("Screened", "Positives"))
                     
ggsave(tested_pos_neg, file="figures/tested_totals_pos_neg.jpg",
       height=4, width=5)

## Looking at total screened bees per site and num of positives by genus groups

tested_pos_neg_genus <- parasites_bees %>%   
  filter(Genus == c("Bombus", "Melissodes", "Apis", "Anthophora")) %>% 
  filter(Parasites != "NosemaBombi" & Parasites != "NosemaCeranae") %>% 
  filter(Site != "VC" & Site != "UK" & Site != "SS")%>% 
  filter(!is.na(value)) %>% 
  group_by(Site, Genus, value) %>% 
  summarise(TestedTotals = length(value)) %>% 
  ggplot(aes(x=Site)) +
  geom_bar(aes(y = TestedTotals, fill = value), stat="identity", position = "dodge")+
  theme_minimal() + coord_flip()+
  theme(legend.title=element_blank())+
  scale_fill_discrete(labels = c("Negatives", "Positives"))+
  facet_wrap(~Genus)

ggsave(tested_pos_neg_genus, file="figures/tested_totals_by_genus.jpg",
       height=4, width=5)

spec.net %>% 
  filter(Apidae == 1) %>% 
summarize(TestedTotals = length(Apidae))

par.counts <- table(spec.net$ParasiteRichness[spec.net$Apidae == 1])
par.counts <- par.counts/sum(par.counts)
par.counts <- as.data.frame(par.counts)

p <- ggplot(data=par.counts, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill= "#feb24c") +
  theme_minimal() + coord_flip() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  xlab("Number of parasites/individual") +
  ylab("Proportion of screened bees")
ggsave(p, file="figures/proportion_pos.jpg",
       height=4, width=5)

# Boxplot of sites and  rate of crithidia in bombus species

parasite_prevalence_sites <- spec.net %>%  
  filter(Site != "VC" & Site != "UK" & Site != "SS") %>% 
  filter(Genus == "Bombus") %>% 
  ggplot(aes(x= reorder(Site, Lat, decreasing = TRUE), y= SpCrithidiaBombusParasitismRate, 
              color = GenusSpecies )) +
  geom_box_plot() +
 # geom_bar(stat="identity", position = "dodge") + theme_minimal() +
  labs(x = "Sites", y = "Crithidia Parasitism Rate")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))

ggsave(parasite_prevalence_sites, file="figures/prevalence_by_sites.jpg",
       height=4, width=5)
