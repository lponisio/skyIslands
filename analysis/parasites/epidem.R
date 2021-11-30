## setwd('~/Dropbox (University of Oregon)/skyIslands/')
setwd('analysis/parasites')
library(car)
library(tidyverse)
load('../../data/spec_RBCL_16s.Rdata')
site.sum <- read.csv("../../data/sitestats.csv")
veg <- read.csv("../../data/veg.csv")


## setwd('analysis\\parasites')
## load('..\\..\\data\\spec_RBCL_16s.Rdata')
## site.sum <- read.csv("..\\..\\data\\sitestats.csv")

## get data in the correct format

parasites <- c("AspergillusSpp",
               "AscosphaeraSpp",
               "ApicystisSpp",
               "CrithidiaExpoeki",
               "CrithidiaBombi",
               "NosemaBombi",
               "NosemaCeranae")

## subset to the years we did the screenings
spec <- spec[spec$Year == "2018",]

## subset to the bees we screened and the screenings worked
spec <- spec[spec$Apidae == 1 &
             !is.na(spec$Apidae),]

spec <- merge(spec, site.sum)

spec <- merge(spec, veg, all.x=TRUE)

check1 <- spec[is.na(spec$FloralRichness)]

table(spec$GenusSpecies)
table(spec$Genus)

## split up by genus
bombus <- spec[spec$Genus == "Bombus",]
table(bombus$Species)
table(bombus$Site)

megachile <- spec[spec$Genus == "Megachile",]
table(megachile$Species)
table(megachile$Site)

anthophora <- spec[spec$Genus == "Anthophora",]
table(anthophora$Species)

apis <- spec[spec$Genus == "Apis",]
table(apis$Site)


## check for any NAs
apply(spec[, parasites], 2, function(x) sum(is.na(x)))

## take means for plotting
## take the mean parasite prevalence for each species
para.gensp <- spec %>%
    group_by(GenusSpecies, Site) %>%
    summarise(AspergillusSpp = mean(AspergillusSpp),
              AscosphaeraSpp = mean(AscosphaeraSpp),
              ApicystisSpp = mean(ApicystisSpp),
              CrithidiaExpoeki =mean(CrithidiaExpoeki),
              CrithidiaBombi=mean(CrithidiaBombi),
              NosemaBombi=mean(NosemaBombi),
              NosemaCeranae=mean(NosemaCeranae))


boxplot(para.gensp$AspergillusSpp~para.gensp$GenusSpecies)
boxplot(para.gensp$AscosphaeraSpp)
boxplot(para.gensp$ApicystisSpp)
boxplot(para.gensp$CrithidiaExpoeki)

## some preliminary models to check out data

library(lmer)
library(lmerTest)

all.mod <- lmer(FloralDiversity ~ scale(Elev) +
                    scale(Lat) + (1|Site),
                data=spec)

summary(all.mod)
vif(all.mod)

bomb.mod <- glm(ParasitePresence ~ scale(PollAbundance) +
                    scale(PollDiversity) + scale(FloralDiversity) +
                    scale(FloralAbundance),
    data=bombus, family="binomial")
summary(bomb.mod)

vif(bomb.mod)


apis.mod <- glm(ParasitePresence ~ scale(PollAbundance) +
                    scale(PollDiversity)  + scale(FloralDiversity) +
                    scale(FloralAbundance),
    data=apis, family="binomial")
summary(apis.mod)

vif(apis.mod)

anthophora.mod <- glm(ParasitePresence ~ PollAbundance +
                          PollDiversity  + scale(FloralDiversity) +
                    scale(FloralAbundance) ,
    data=anthophora, family="binomial")
summary(anthophora.mod)
vif(anthophora.mod)

megachile.mod <- glm(ParasitePresence ~ PollAbundance +
                         PollDiversity + scale(FloralDiversity) +
                    scale(FloralAbundance),
    data=megachile, family="binomial")
summary(megachile.mod)
vif(megachile.mod)



all.mod <- glm(ParasitePresence ~ PollAbundance + PollDiversity +
                   scale(FloralDiversity) +
                    scale(FloralAbundance),
    data=spec, family="binomial")
summary(all.mod)

vif(all.mod)
