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

sick.totals <- spec.all %>%
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
ggplot(spec.all, aes(x= Lat,
                     y=Net_HBAbundance, color=Site))+
  geom_point()

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

## 
