rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
## setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
setwd('~/Dropbox (University of Oregon)/skyislands')

library(gplots)
library(bipartite)
library(RColorBrewer)
library(tidyverse)

setwd("analysis/parasites")
source("src/misc.R")
load(file="saved/spec_weights.Rdata")

## screened
rownames(spec.net) <- NULL
spec.screened <- spec.net[spec.net$WeightsSp ==1,]

## no positives for Nosema bombi, only 1 for ceranae

parasites <- c("AscosphaeraSpp", "ApicystisSpp",
               "CrithidiaExpoeki", "CrithidiaMellificae",
               "CrithidiaBombi", "CrithidiaSpp")
parasite.cols <- c( "SpCrithidiaPresence",
                   paste0("Sp", parasites))
spec.screened <- spec.screened[, c("GenusSpecies", "Genus",
                                   "MtRange",
                                   "SampleRound", "Year",
                                   "SpScreened",
                                   parasite.cols)]


## sum over sample rounds, years
sum.genus.screened <- spec.screened  %>%
    group_by(MtRange, Genus) %>%
    summarise(
        SpCrithidiaPresence= sum(SpCrithidiaPresence),
        SpApicystisSpp = sum(SpApicystisSpp),
        SpAscosphaeraSpp= sum(SpAscosphaeraSpp),
        SpCrithidiaExpoeki = sum(SpCrithidiaExpoeki),
        SpCrithidiaBombi = sum(SpCrithidiaBombi),
        SpCrithidiaMellificae = sum(SpCrithidiaMellificae),
        SpCrithidiaSpp = sum(SpCrithidiaSpp),
        SpScreened = sum(SpScreened)
    )

sum.screened <- spec.screened  %>%
    group_by(MtRange, GenusSpecies, Genus) %>%
    summarise(
        SpCrithidiaPresence= sum(SpCrithidiaPresence),
        SpApicystisSpp = sum(SpApicystisSpp),
        SpAscosphaeraSpp= sum(SpAscosphaeraSpp),
        SpCrithidiaExpoeki = sum(SpCrithidiaExpoeki),
        SpCrithidiaBombi = sum(SpCrithidiaBombi),
        SpCrithidiaMellificae = sum(SpCrithidiaMellificae),
        SpCrithidiaSpp = sum(SpCrithidiaSpp),
        SpScreened = sum(SpScreened)
    )

## sum by species
sp.sum <- spec.screened  %>%
    group_by(GenusSpecies) %>%
    summarise(
        TotalScreened = sum(SpScreened)
    )
## sum by specues
gen.sum <- spec.screened  %>%
    group_by(Genus) %>%
    summarise(
        TotalScreened = sum(SpScreened)
    )
## sum by site
site.sum <- spec.screened  %>%
    group_by(MtRange) %>%
    summarise(
        TotalScreened = sum(SpScreened)
    )

## create proportions
sum.screened[, parasite.cols ] <- sum.screened[, parasite.cols
                                               ]/sum.screened$SpScreened
sum.genus.screened[, parasite.cols ] <-
    sum.genus.screened[, parasite.cols ]/sum.genus.screened$SpScreened

## add N to species/genus
sum.screened$GenusSpecies <- paste0(sum.screened$GenusSpecies, " (",
                                    sp.sum$TotalScreened[match(
                                        sum.screened$GenusSpecies,
                                        sp.sum$GenusSpecies)],
                                    ")")

sum.genus.screened$Genus <- paste0(
    sum.genus.screened$Genus, " (",
    gen.sum$TotalScreened[match(
        sum.genus.screened$Genus,
        gen.sum$Genus)],
    ")")


## add N to sites
sum.screened$MtRange <- paste0(sum.screened$MtRange, " (",
                               site.sum$TotalScreened[match(
                                   sum.screened$MtRange,
                                   site.sum$MtRange)],
                               ")")

sum.genus.screened$MtRange <- paste0(
    sum.genus.screened$MtRange, " (",
    site.sum$TotalScreened[match(
        sum.genus.screened$MtRange,
        site.sum$MtRange)],
    ")")

sum.screened$SpScreened <- NULL
sum.genus.screened$SpScreened <- NULL

## make community matrices 
makeParMat <- function(parasite, screened, sp.col="GenusSpecies"){
    colnames(screened)[colnames(screened) == parasite] <-
        "parasite"
    colnames(screened)[colnames(screened) == sp.col] <- "SpCat"
    par.mat <- screened %>%
        select(SpCat, MtRange, parasite) %>%
        pivot_wider(names_from = SpCat, values_from =
                                            parasite)
    par.mat <- as.data.frame(par.mat)
    rownames(par.mat) <- par.mat$MtRange
    par.mat$MtRange <- NULL
    par.mat <- as.matrix(par.mat)
    return(par.mat)
}

par.mats.species <- lapply(parasite.cols, makeParMat, sum.screened)
names(par.mats.species) <- parasite.cols
par.mats.genus <- lapply(parasite.cols, makeParMat,
                         sum.genus.screened, sp.col="Genus")
names(par.mats.genus) <- parasite.cols

## heat maps of # of infected individuals by species
plotParasiteMap <- function(){
    colfunc <- colorRampPalette(c("grey", "red"))
    if(parasite == "SpAscosphaeraSpp"){
        par(oma=c(15,4,3,10), mar=c(1,2,2,1),
            mgp=c(1.5,0.5,0))
    } else{
        par(oma=c(10,4,3,8), mar=c(1,2,2,1),
            mgp=c(1.5,0.5,0))
    }
    heatmap.2(bipartite::empty(par.mats.species[[parasite]]),
              trace="none",
              col=colfunc,
              breaks=seq(0, 1, 0.1))
    mtext(parasite, 3, line=2.5)
}
for(parasite in parasite.cols){
    pdf.f(plotParasiteMap,
          file=sprintf("figures/heatmaps/%s.pdf", parasite),
          width=10, height=7)
}

## heat maps of # of infected individuals by genus
plotParasiteMapGenus <- function(){
    colfunc <- colorRampPalette(c("grey", "red"))
    par(oma=c(11,4,3,8), mar=c(4,2,2,1),
        mgp=c(1.5,0.5,0))
    heatmap.2(bipartite::empty(par.mats.genus[[parasite]]),
              trace="none",
              col=colfunc,
              breaks=seq(0, 1, 0.1))
    mtext(parasite, 3, line=2.5)
}
for(parasite in parasite.cols){
    pdf.f(plotParasiteMapGenus,
          file=sprintf("figures/heatmaps/genus_%s.pdf", parasite),
          width=10, height=7)
}

## ***************************************************************
## barplots - by genus
## ***************************************************************
## sum over genera
sum.genus.screened <- spec.screened  %>%
    group_by(Genus) %>%
    summarise(
        ## SpCrithidiaPresence= sum(SpCrithidiaPresence),
        SpCrithidiaExpoeki = sum(SpCrithidiaExpoeki),
        SpCrithidiaBombi = sum(SpCrithidiaBombi),
        SpCrithidiaMellificae = sum(SpCrithidiaMellificae),
        SpCrithidiaSpp = sum(SpCrithidiaSpp),
        SpApicystisSpp = sum(SpApicystisSpp),
        ## SpNosemaBombi= sum(SpNosemaBombi),
        ## SpNosemaCeranae = sum(SpNosemaCeranae),
        SpAscosphaeraSpp= sum(SpAscosphaeraSpp),
        SpScreened = sum(SpScreened)
    )

sum.genus.screened[, parasite.cols[-1] ] <-
    sum.genus.screened[, parasite.cols[-1] ]/
    sum.genus.screened$SpScreened

## sum.genus.screened[, parasite.cols ] <-
##     sum.genus.screened[, parasite.cols ]/
##     sum.genus.screened$SpScreened


sum.genus.screened$SpScreened <- NULL

long.sum <- sum.genus.screened %>%
    pivot_longer(
        cols = starts_with("Sp"),
        names_to = "Parasite",
        values_drop_na = TRUE
    )
long.sum <- long.sum[long.sum$Genus != "Dufourea",]

library(RColorBrewer)
pp <- ggplot(long.sum, aes(Genus, value, fill=Parasite)) +
    labs(y = "Proportion tested positive") +
    scale_fill_brewer(palette = "Dark2",
                      name = "Parasite",
                      labels = c("Apicystis spp.",
                                 "Ascosphaera spp.",
                                 "Crithidia bombi",
                                 "Crithidia expoeki",
                                 "Crithidia mellificae",
                                 "Crithidia spp."))
pp <- pp + geom_bar(stat = "identity",
                    position = 'dodge')

ggsave(pp, file="figures/parasite_barplots.pdf",
       height=2.5, width=7)

## ***************************************************************
## boxplots by species
## ***************************************************************

long.sum <- sum.screened %>%
    pivot_longer(
        cols = starts_with("Sp"),
        names_to = "Parasite",
        values_drop_na = TRUE
    )
long.sum <- long.sum[long.sum$Genus != "Dufourea",]
long.sum <- long.sum[!long.sum$Parasite %in%
                     c("SpCrithidiaExpoeki", "SpCrithidiaMellificae",
                       "SpCrithidiaBombi", "SpCrithidiaSpp"),]

pp2 <- ggplot(long.sum, aes(Genus, value, fill=Parasite)) +
    geom_boxplot() + 
    labs(y = "Proportion tested positive") +
    scale_fill_brewer(palette = "Dark2",
                      name = "Parasite",
                      labels = c("Apicystis spp.",
                                 "Ascosphaera spp.",
                                 "Crithidia spp.")) 
ggsave(pp2, file="figures/parasite_boxplots.pdf",
       height=2.5, width=7)
