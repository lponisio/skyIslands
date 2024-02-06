library(tidyverse)

## Plots by meadow
## Bee abundances by meadows
###############################################################################
##Bombus

#ggplot(spec.all, aes(Site, Net_BombusAbundance))+ geom_boxplot()

ggplot(spec.all, aes(Site, Net_BombusAbundance)) + geom_bar(stat = "identity")

##Apis
#ggplot(spec.all, aes(Site, Net_HBAbundance))+ geom_boxplot()

ggplot(spec.all, aes(Site, Net_HBAbundance)) + geom_bar(stat = "identity")

##Other bees
#ggplot(spec.all, aes(Site, Net_NonBombusHBAbundance))+ geom_boxplot()

ggplot(spec.all, aes(Site, Net_NonBombusHBAbundance)) + geom_bar(stat = "identity")
################################################################################
## Bee Diversity by meadows
ggplot(spec.all, aes(Site, Net_BeeDiversity)) + geom_boxplot()
## Bombus Diversity
ggplot(spec.all, aes(Site, Net_BombusDiversity)) + geom_boxplot()

###############################################################################
## Floral abundance by meadows
ggplot(spec.all, aes(Site, MeanFloralAbundance)) + geom_boxplot()

## Floral diversity by meadows

ggplot(spec.all, aes(Site, MeanFloralDiversity)) + geom_boxplot()

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
ggplot(spec.all, aes(x= Year,
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

ggsave(parasite_prevalence, file="figures/parasite_prevalence_by_genus.pdf",
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

ggsave(parasite_prevalence_sites, file="figures/prevalence_by_sites.pdf",
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

ggsave(parasite_prevalence_genus, file="figures/prevalence_by_genus.pdf",
       height=4, width=5)

## Looking at total screened bees per site and num of positives
parasites_bees$value <- as.factor(parasites_bees$value)
parasites_bees %>%   
  filter(Genus == c("Bombus", "Melissodes", "Apis", "Anthophora")) %>% 
  filter(Parasites != "NosemaBombi" & Parasites != "NosemaCeranae") %>% 
  filter(Site != "VC" & Site != "UK" & Site != "SS")%>% 
  filter(!is.na(value)) %>% 
  group_by(Site, Genus, value) %>% 
  summarise(TestedTotals = length(value)) %>% 
  ggplot(aes(x=Site)) +
  geom_bar(aes(y = TestedTotals, fill = value), stat="identity", position = "dodge")+
  theme_minimal() + coord_flip()+
  facet_wrap(~Genus)
  