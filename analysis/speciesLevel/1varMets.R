## setwd('~/Dropbox/skyIslands')
setwd('analysis/speciesLevel')
## net.type <- "YrSR"
net.type <- "Yr"
species <- c("Plant", "Pollinator")
## species <- c("Pollinator", "Parasite")
source('src/initialize.R')
source('src/calcVar.R')
source('src/plotting.R')

metrics <- c("species.strength",
             "d",
             "betweenness",
             "closeness",
             "rare.degree",
             "niche.overlap")

## formulas for the relationship between the var/mean species-level
## network metrics and the number of occurrences
if(species[1] == "Plant"){
    mod.FUN <- lmer
    formulas.pv <-lapply(paste(metrics, "pv", sep="."),
                         function(y) {
                             as.formula(paste(y, "~",
                             "scale(maxN) + (1|Year.Site) + (1|GenusSpecies)"))
                         })

    formulas.mean <-lapply(metrics,
                           function(y) {
                               as.formula(paste(y, "~",
                     "scale(maxN) +  (1|Site) + (1|GenusSpecies)"))
                           })
} else if(species[2] == "Parasite"){
    mod.FUN <- lm
    formulas.pv <-lapply(paste(metrics, "pv", sep="."),
                         function(y) {
                             as.formula(paste(y, "~",
                                              "scale(maxN)"))
                         })
    formulas.mean <-lapply(metrics,
                           function(y) {
                               as.formula(paste(y, "~",
                                                "scale(maxN) + GenusSpecies"))
                           })
}

## *****************************************************************
## mean/variability between sites within a year
## *****************************************************************

by.yr <- split(sp.lev, sp.lev$Year)
by.sp.geo <- lapply(by.yr, function(x) {split(x, x$GenusSpecies)})
by.sp.geo <- unlist(by.sp.geo, recursive=FALSE)

geo.met <- calcVar(by.sp.geo, PV)
sp.lev$maxN <- geo.met$maxN[match(sp.lev$GenusSpecies,
                                  geo.met$GenusSpecies)]

## higher level (pollinators/parasites)
mods.pv.higher <- lapply(formulas.pv, function(x){
    try(mod.FUN(x, data=geo.met[geo.met$speciesType == "higher.level",],
                REML = FALSE), silent=TRUE)
})


mods.mean.higher <- lapply(formulas.mean, function(x){
    mod.FUN(x, data=sp.lev[sp.lev$speciesType == "higher.level",],
            REML = FALSE)
})

names(mods.pv.higher) <- names(mods.mean.higher) <- metrics

mods.pv.higher <- mods.pv.higher[!sapply(mods.pv.higher, function(x)
    inherits(x, "try-error"))]

## results
## lapply(mods.pv.higher, summary)
## check sig levels with method other than wald CI
lapply(mods.pv.higher, anova)

##  pollinator
## most variable in rare.degree, + with occurrence

##  parasite
## nothing sig

## lapply(mods.mean.higher, summary)
## check sig levels with method other than wald CI
lapply(mods.mean.higher, anova)

##  pollinator
## all sig with occurrence

##  parasite
## species strength, niche overlap, closeness + with occurrence

## lower
mods.pv.lower <- lapply(formulas.pv, function(x){
    try(mod.FUN(x, data=geo.met[geo.met$speciesType == "lower.level",],
            REML = FALSE), silent=TRUE)
})
mods.mean.lower <- lapply(formulas.mean, function(x){
    mod.FUN(x, data=sp.lev[sp.lev$speciesType == "lower.level",])
})
names(mods.pv.lower) <- names(mods.mean.lower) <- metrics

mods.pv.lower <- mods.pv.lower[!sapply(mods.pv.lower, function(x)
    inherits(x, "try-error"))]

## results
## lapply(mods.pv.lower, summary)
## check sig levels with method other than wald CI
lapply(mods.pv.lower, anova)

## plant
## no sig relationships with occurrence

## pollinator
## rare.degree +  occurrence

## lapply(mods.mean.lower, summary)
## check sig levels with method other than wald CI
lapply(mods.mean.lower, anova)

## plant
## species strength
## d
## niche overlap

## pollinator
## none sig

plotAllGeo()


## *****************************************************************
## mean/variability between years within a site
## *****************************************************************
## if(species[1] == "Plant"){
##     by.site <- split(sp.lev, sp.lev$Site)
##     by.sp.site <- lapply(by.site, function(x) {split(x, x$GenusSpecies)})
##     by.sp.site <- unlist(by.sp.site, recursive=FALSE)

##     yr.met <- calcVar(by.sp.site, PV)

##     ## pollinators
##     mods.pv.higher <- lapply(formulas.pv, function(x){
##         lmer(x, data=yr.met[yr.met$speciesType == "higher.level",])
##     })
##     mods.mean.higher <- lapply(formulas.mean, function(x){
##         lmer(x, data=yr.met[yr.met$speciesType == "higher.level",])
##     })

##     names(mods.pv.higher) <- names(mods.mean.higher) <- metrics
##     ## results
##     lapply(mods.pv.higher, summary)
##     ## check sig levels with method other than wald CI
##     lapply(mods.pv.higher, anova)
##     ## most variable in betweenness + number of years detected at a site
##     ## most variable in degree + number of years detected at a site

##     lapply(mods.mean.higher, summary)
##     ## check sig levels with method other than wald CI
##     lapply(mods.mean.higher, anova)
##     ##  all metrics + number of years detected at a site


##     ## lowers
##     mods.pv.lower <- lapply(formulas.pv, function(x){
##         mod.FUN(x, data=yr.met[yr.met$speciesType == "lower.level",])
##     })
##     mods.mean.lower <- lapply(formulas.mean, function(x){
##         mod.FUN(x, data=yr.met[yr.met$speciesType == "lower.level",])
##     })
##     names(mods.pv.lower) <- names(mods.mean.lower) <- metrics
##     ## results
##     lapply(mods.pv.lower, summary)
##     ## check sig levels with method other than wald CI
##     lapply(mods.pv.lower, anova)
##     ## most variable in species strength + number of years detected at a site
##     ## most variable in betweenness  + number of years detected at a site

##     lapply(mods.mean.lower, summary)
##     ## check sig levels with method other than wald CI
##     lapply(mods.mean.lower, anova)
##     ## betweenness
##     ## closeness
##     ## niche overlap

##     plotAllYr()
## }

## ## variability for between years maybe doesn't make any sense with
## ## only three years of data.
