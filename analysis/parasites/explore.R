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
 spec.orig <- filter(spec.orig, Site != "VC")

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
  #filter(Apidae == 1) %>% 
  group_by(GenusSpecies) %>% 
  summarize(n = n_distinct(Site)) %>% 
  filter(n == 1) 

## Summary numbers for veg
veg <- filter(veg, Site != "VC")
veg %>%
  group_by(PlantGenusSpecies) %>% 
  summarize(n = n_distinct(Site)) %>% 
  filter(n == 2)
## Plots by meadow
## Bee abundances by meadows
###############################################################################
##Bombus
bombus_abundance <- spec.orig %>%  
ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), Net_BombusAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  scale_color_manual(values = c("#08519c","#6baed6","#fd8d3c","darkgoldenrod3", "goldenrod1")) +
  labs(x = "Meadows", y = "Bombus abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size=12),
        text = element_text(size=12))
bombus_abundance  


##Apis
HB_abundance <- spec.orig %>% 
  ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), Net_HBAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  scale_color_manual(values = c("#08519c","#6baed6","#fd8d3c","darkgoldenrod3", "goldenrod1")) +
  labs(x = "Meadows", y = "Apis abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size=12),
        text = element_text(size=12))
HB_abundance  


## bee abundance
bee_abundance <- spec.orig %>% 
  ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), Net_BeeAbundance))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  scale_color_manual(values = c("#08519c","#6baed6","#fd8d3c","darkgoldenrod3", "goldenrod1")) +
  labs(x = "Meadows", y = "Bee abundance", color = "Year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size=12),
        text = element_text(size=12))
  bee_abundance

################################################################################
## Bee Diversity by meadows
bee_diversity <- spec.orig %>% 
    ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), Net_BeeDiversity))+ 
    geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
    scale_color_manual(values = c("#08519c","#6baed6","#fd8d3c","darkgoldenrod3", "goldenrod1")) +
    labs(x = "Meadows", y = "Bee diversity", color = "Year")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=12),
          text = element_text(size=12))
bee_diversity


###############################################################################
## Floral abundance by meadows
floral_abundance <- spec.orig %>% 
    ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), MeanFloralAbundance))+ 
    geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  scale_color_manual(values = c("#08519c","#6baed6","#fd8d3c","darkgoldenrod3", "goldenrod1")) +
    labs(x = "Meadows", y = "Mean floral abundance", color = "Year")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=12),
          text = element_text(size=12))
floral_abundance


## Floral diversity by meadows

floral_diversity<- spec.orig %>% 
  ggplot(aes(reorder(MtRange, Lat, decreasing = TRUE), MeanFloralDiversity))+ 
  geom_boxplot()+ geom_point(aes(color = as.factor(Year)))+
  scale_color_manual(values = c("#08519c","#6baed6","#fd8d3c","darkgoldenrod3", "goldenrod1")) +
  labs(x = "Meadows", y = "Mean floral diversity", color = "Year")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        text = element_text(size=12))
floral_diversity

meadows_summary <- ggarrange(floral_diversity, floral_abundance, bee_diversity, 
                             bee_abundance, bombus_abundance, HB_abundance,
                             nrow = 2, ncol = 3, 
                             labels = c("A", "B", "C", "D", "E", "F"), 
                             common.legend = T,
                             legend = "right")

ggsave(meadows_summary, file = "figures/meadows_summary.pdf", height = 8, width = 12)

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
            InfectedCrithdiaMellificae = round(mean(CrithidiaMellificae,
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
ggplot(spec.orig, aes(y=Net_BombusAbundance,
                     x=Net_HBAbundance))+
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

## Relationship between bombus abundance and diet breadth
bombus_abund_degree <- ggplot(spec.net, aes(x= rare.degree,
                     y=Net_BombusAbundance))+
  geom_point()
ggsave(bombus_abund_degree, file="figures/bombus_abund_degree.pdf",
       height=4, width=5)

## Relationship between hb abundance and diet breadth
hb_abund_degree <- ggplot(spec.net, aes(x= rare.degree,
                     y=Net_HBAbundance))+
  geom_point()
ggsave(hb_abund_degree, file="figures/hb_abund_degree.pdf",
       height=4, width=5)

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


parasites_bees <- pivot_longer(spec.orig,
                              cols = c(ApicystisSpp, AscosphaeraSpp,
                                       CrithidiaPresence, 
                                       CrithidiaExpoeki,
                                       CrithidiaMellificae, CrithidiaBombi), 
                              names_to = "Parasites")

parasites_bees <- parasites_bees[!parasites_bees$Genus %in%
                                 c("Dufourea", "Pseudopanurgus",
                                 "Colletes", ""),]
parasites_bees %>% select(Parasites, value) %>% 
  filter(!is.na(value)) %>% 
  group_by(Parasites) %>% 
  summarise(Parasitism = sum(value, na.rm = TRUE)) %>% 
  ggplot()+
  geom_bar(aes(Parasites))

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

parasite_prevalence_sites <- spec.orig %>%  
  filter(Site != "VC" & Site != "UK" & Site != "SS") %>% 
  filter(Genus == "Bombus") %>% 
  ggplot(aes(x= reorder(Site, Lat, decreasing = TRUE), y= SpCrithidiaBombusParasitismRate, 
              color = GenusSpecies )) +
  geom_boxplot() +
 # geom_bar(stat="identity", position = "dodge") + theme_minimal() +
  labs(x = "Sites", y = "Crithidia Parasitism Rate")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), 
        axis.title.y = element_text(size=10),
        text = element_text(size=10))

ggsave(parasite_prevalence, file="figures/parasite_boxplots.pdf",
       height=5, width=8)

## Parasite prevalence by genus

parasite_prevalence <- parasites_bees %>%   
  filter(!is.na(value)) %>% 
  filter(Order == "Hymenoptera" & Family != "Vespidae" & Family != "Sphecidae") %>% 
  filter(Apidae == 1) %>% 
  group_by(Site, Genus, Parasites) %>% 
  summarise(total_screened = length(UniqueID[!is.na(ParasitePresence)]),              # Total number of individuals screened
            positives = sum(value),         # Total positives
            prop_positive = (positives / total_screened)  # Percentage of positives
            ) %>% 
  ggplot(aes(Genus, prop_positive)) +
  geom_boxplot(aes(fill = Parasites)) +
  labs(x = "Genus", y = "Proportion tested positive")+
  theme(axis.text.x = element_text(size=10), 
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        text = element_text(size=15))+
  scale_fill_manual(values = c("darkgoldenrod3", "grey", "#3182bd"),
                    labels = c("Apicystis spp.", "Ascosphaera spp.", "Crithidia spp."))
  
ggsave(parasite_prevalence, file="figures/parasite_boxplots.pdf",
       height=4, width=7)


spec.orig <- spec.orig[!spec.orig$Genus%in%
           c("Dufourea", "Pseudopanurgus",
             "Colletes", ""),]

par_prev <- spec.orig %>%
  filter(WeightsPar == 1) %>% 
  ## Subset to Weights == 1
  #group_by(Genus) %>%
  summarise(TestedTotals = length(UniqueID[Apidae == 1]),
            InfectedApicystisSpp= round(sum(ApicystisSpp, na.rm=TRUE)/TestedTotals * 100, 2),
            InfectedAscosphaeraSpp=round(sum(AscosphaeraSpp, na.rm=TRUE)/TestedTotals *100, 2),
            InfectedPresence=round(sum(CrithidiaPresence, na.rm=TRUE)/TestedTotals *100, 2),
            InfectedCrithidiaBombi=round(sum(CrithidiaBombi, na.rm=TRUE)/TestedTotals *100, 2),
            InfectedCrithidiaExpoeki=round(sum(CrithidiaExpoeki, na.rm=TRUE)/TestedTotals *100, 2),
            InfectedCrithdiaMellificae = round(sum(CrithidiaMellificae, na.rm=TRUE)/TestedTotals *100, 2),
)
sp<-spec.orig %>% 
  filter(WeightsPar == 1) %>%
  group_by(GenusSpecies) %>% 
  summarise(n = n())

