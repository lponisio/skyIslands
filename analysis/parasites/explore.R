rm(list=ls())
library(tidyverse)
library(patchwork)
setwd("~/")
source("lab_paths.R")
local.path
dir.bombus <- file.path(local.path, "skyIslands/analysis/parasites")
setwd(dir.bombus)
load(file="saved/spec_weights.Rdata")
veg <- read_csv("../../../skyIslands_saved/data/relational/original/veg.csv")
## These sites were only visited once. 
## spec.orig <- filter(spec.orig, Site != "VC")

## Chiricahua and Sacramento has two meadows so renamed them to be able to show 
## both meadows in all the community graphs. 
for(i in 1:nrow(spec.orig)){
    if(spec.orig$Site[i] == "CH" ){
        spec.orig$MtRange[i] <- "Chiricahua A"
    } else if(spec.orig$Site[i] == "RP"){
        spec.orig$MtRange[i] <- "Chiricahua B"
    } else  if(spec.orig$Site[i] == "UK" ){
        spec.orig$MtRange[i] <- "Sacramento A"
    } else if(spec.orig$Site[i] == "SS"){
        spec.orig$MtRange[i] <- "Sacramento B"
    }
}

spec.orig %>% 
    group_by(MtRange, Meadow, Site, Lat, SampleRound) %>% 
    summarize (n = n())

## Used this to get the summary of number of species/morphospecies total, 
## genus with > 10 species, and species that were found in >7 sites.
num_per_site <- spec.orig %>% 
  filter(Order == "Hymenoptera" & Family != "Vespidae" & Family != "Sphecidae") %>% 
  #group_by(GenusSpecies) %>% 
  #summarize(n = n_distinct(Site)) %>% 
  #filter(n > 7) %>% 
  summarize(n = n())

## Summary numbers for veg
veg <- filter(veg, Site != "VC")
veg %>%
  group_by(PlantGenusSpecies) %>% 
  summarize(n = n_distinct(Site)) %>% 
  filter(n > 5)
## Plots by meadow
## Bee abundances by meadows
###############################################################################
##Bombus
bombus_abundance <- spec.orig %>%  
ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), Net_BombusAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Meadows", y = "Bombus Abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))
bombus_abundance  


##Apis
HB_abundance <- spec.orig %>% 
  ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), Net_HBAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Meadows", y = "Apis Abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))
HB_abundance  


## melissodes
melissodes_abundance <- spec.orig %>% 
  filter(Genus == "Melissodes") %>% 
  ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), Net_NonBombusHBAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Meadows", y = "Melissodes Abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))
  melissodes_abundance

################################################################################
## Bee Diversity by meadows
bee_diversity <- spec.orig %>% 
    ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), Net_BeeDiversity))+ 
    geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
    labs(x = "Meadows", y = "Bee Diversity", color = "Year")+
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
          axis.title.y = element_text(size=10),
          text = element_text(size=10))
bee_diversity


###############################################################################
## Floral abundance by meadows
floral_abundance <- spec.orig %>% 
    ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), MeanFloralAbundance))+ 
    geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
    labs(x = "Meadows", y = "Mean Floral Abundance", color = "Year")+
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
          axis.title.y = element_text(size=10),
          text = element_text(size=10))
floral_abundance


## Floral diversity by meadows

floral_diversity<- spec.orig %>% 
  ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), MeanFloralDiversity))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  labs(x = "Meadows", y = "Mean Floral Diversity", color = "Year")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))
floral_diversity

bombus_abundance + HB_abundance + melissodes_abundance + bee_diversity + 
  floral_abundance + floral_diversity + plot_layout(ncol = 3)+ 
  plot_annotation(tag_levels = "A")+ plot_layout(guides = "collect", axis_titles = "collect")


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


parasites_bees <- pivot_longer(spec.net,
                              cols = c(ApicystisSpp, AscosphaeraSpp,
                                       CrithidiaPresence), names_to =
                              "Parasites")

parasites_bees <- parasites_bees[!parasites_bees$Genus %in%
                                 c("Dufourea", "Pseudopanurgus",
                                 "Colletes", ""),]

parasite_prevalence <- parasites_bees %>%   
    ## filter(Genus == c("Melissodes", "Anthophora", "Apis", "Bombus")) %>% 
    filter(!is.na(value)) %>% 
    group_by(Site, Genus, Parasites) %>% 
    summarise(ParasitismRate = mean(value, na.rm = TRUE)) %>% 
    ggplot(aes(Parasites, ParasitismRate))+
    geom_boxplot()+
    labs(x = "Parasites", y = "Parasitism Rate")+
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10), 
          axis.title.y = element_text(size=10),
          text = element_text(size=10)) + facet_wrap(~ Genus)

ggsave(parasite_prevalence, file="figures/parasite_prevalence_by_genus.jpg",
       height=4, width=5)

## Bar plot of parasite prevalence per sites for each genus
parasite_prevalence_sites <- parasites_bees %>%   
  ## filter(Genus == c("Bombus", "Melissodes", "Apis", "Anthophora")) %>% 
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


tested_pos_neg<- spec.orig %>%   
  filter(Genus == "Melissodes" | Genus == "Bombus"| Genus == "Apis") %>% 
  filter(Apidae == 1) %>% 
  group_by(MtRange) %>% 
  summarize(TestedTotals = sum(Apidae, na.rm = TRUE),
            Positives = sum(ParasitePresence,na.rm=TRUE)) %>% 
  pivot_longer(cols = c(TestedTotals, Positives), 
               names_to = "Category", 
               values_to = "Screenings") %>% 
  ggplot(aes(x= MtRange, fill = Category)) +
  geom_bar(aes(y = Screenings), stat ="identity", 
           position = "dodge")+
  labs(x = "Meadows", y = "Number of Screened Bees")+
  theme_minimal() + coord_flip()+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#feb24c", "grey"), 
                    labels = c("Positives", "Screened")) 
                     
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
