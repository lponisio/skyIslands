## setwd('~/Dropbox (University of Oregon)/skyIslands/')
setwd('analysis/parasites')
library(car)
load('../../data/spec_RBCL_16s.Rdata')
site.sum <- read.csv("../../data/sitestats.csv")

parasites <- c("AspergillusSpp",
               "AscosphaeraSpp",
               "ApicystisSpp",
               "CrithidiaExpoeki",
               "CrithidiaBombi",
               "NosemaBombi",
               "NosemaCeranae")

spec <- spec[spec$Year == "2018",]
spec <- spec[spec$Apidae == 1 &
             !is.na(spec$Apidae),]

spec <- merge(spec, site.sum)

table(spec$GenusSpecies)
table(spec$Genus)

bombus <- spec[spec$Genus == "Bombus",]
table(bombus$Species)

megachile <- spec[spec$Genus == "Megachile",]
table(megachile$Species)

anthophora <- spec[spec$Genus == "Anthophora",]
table(anthophora$Species)

apis <- spec[spec$Genus == "Apis",]


## check for any NAs
apply(spec[, parasites], 2, function(x) sum(is.na(x)))

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


bomb.mod <- glm(ParasitePresence ~ TotalAbundance + Richness,
    data=bombus, family="binomial")
summary(bomb.mod)

vif(bomb.mod)


apis.mod <- glm(ParasitePresence ~ scale(TotalAbundance) +
                    scale(Richness),
    data=bombus, family="binomial")
summary(apis.mod)

vif(apis.mod)

anthophora.mod <- glm(ParasitePresence ~ TotalAbundance + Richness,
    data=anthophora, family="binomial")
summary(anthophora.mod)
vif(anthophora.mod)

megachile.mod <- glm(ParasitePresence ~ TotalAbundance + Richness,
    data=megachile, family="binomial")
summary(megachile.mod)
vif(megachile.mod)



all.mod <- glm(ParasitePresence ~ TotalAbundance + Richness,
    data=spec, family="binomial")
summary(all.mod)

vif(all.mod)
