#Network building for 16s and rbcL FFAR 2019
setwd('~/Dropbox (University of Oregon)')
#dir.bombus <- '/Volumes/bombus/rhayes/Dropbox (University of Oregon)'

#setwd(dir.bombus)
rm(list=ls())
## setwd('~')
setwd('skyIslands_saved')
library(bipartite)
library(parallel)
source('../skyIslands/dataPrep/src/nets_init.R')

ndim.func <- 10
options(cores=1)

## ********************************************************
## pollinator abundance by site to multiple interactions by
## ********************************************************
pol.abund <- read.csv(file.path(mainDir, "spstats_net.csv"))

## prep for making a vector of bee abund at each site
pol.abund.site <- split(pol.abund, pol.abund$Site)

## ********************************************************
## RBCL
## ********************************************************

## list of rbcl data columns
rbcl_names <- grep('RBCL',names(spec.net))

# note check out RBCL:NA

## split out specimens that were sequenced
seq.df <- spec.net[apply(spec.net[, rbcl_names], 1, function(x) !all(is.na(x))),]

## makes community matrix/network of plants by individuals
indivNet_rbcl <- illumSplit(seq.df, "Site", rbcl_names)
## Format: each element of this list is a large matrix for that site


pol.sp.rbcl <- lapply(indivNet_rbcl,
                function(x){
                spec.net$GenusSpecies[match(colnames(x), spec.net$UniqueID)]
})


pol.abund.site.rbcl <- pol.abund.site[names(indivNet_rbcl)]

## species level network
spNet_rbcl <-  mapply(makeSpNet,
    indiv.nets = indivNet_rbcl,
    sp.names = pol.sp.rbcl,
    sp.abund=pol.abund.site.rbcl,
    SIMPLIFY = FALSE)


## save these matrices for use elsewhere
save(indivNet_rbcl, spNet_rbcl, pol.sp.rbcl,
     file=file.path(save.dir, 'rbclNets.RData'))

## ********************************************************

## calculate summary statistics for each "species" (individual) in the
## networks

## set m to the min number of traits. The Functin's default is to use
## each trait as a dimensions, which get's stuck due to the high
## number of zeros, and "traits"
## see https://stat.ethz.ch/pipermail/r-sig-ecology/2016-January/005264.html

nindiv.rbcl <- sapply(indivNet_rbcl, ncol) > 1

nsp.rbcl <- sapply(spNet_rbcl, ncol) > 1

indivNetSums_rbcl  <-  netSums(indivNet_rbcl[nindiv.rbcl],
                               binary=FALSE,
                               spec=spec.net, FUN=mclapply,
                               m=ndim.func)

indivNetSums_rbclBinary  <-  netSums(indivNet_rbcl[nindiv.rbcl],
                                     binary=TRUE,
                                     spec=spec.net,
                                     FUN=mclapply,
                                     m=ndim.func)

spNetSums_rbcl  <-  netSums(spNet_rbcl[nsp.rbcl],
                            binary=FALSE,
                            spec=spec.net,
                            type="sp")

spNetSums_rbclBinary  <-  netSums(spNet_rbcl[nsp.rbcl],
                                  binary=TRUE,
                                  spec=spec.net,
                                  type="sp")


save(indivNetSums_rbcl, indivNetSums_rbclBinary,
     spNetSums_rbcl,
     spNetSums_rbclBinary,
     file=file.path(save.dir, 'NetSums_rbcl.RData'))

## ********************************************************
## 16s
## ********************************************************
#list of column names for 16s
names_16s  <-  grep('16s',names(spec.net))

#split out specimens that were sequenced
seq.micro <- spec.net[apply(spec.net[, names_16s], 1,
                        function(x) !all(is.na(x))),]

indivNet_micro <- illumSplit(seq.micro,"Site", names_16s)
#Format: each element of this list is a large matrix for that site

pol.sp.micro <- lapply(indivNet_micro,
                function(x){
                    spec.net$GenusSpecies[match(colnames(x),
                                            spec.net$UniqueID)]
})

## species level
pol.abund.site.micro <- pol.abund.site[names(indivNet_micro)]

## species level network
spNet_micro <-  mapply(makeSpNet,
    indiv.nets = indivNet_micro,
    sp.names = pol.sp.micro,
    sp.abund=pol.abund.site.micro,
    SIMPLIFY = FALSE)


#save these matrices for use elsewhere
save(indivNet_micro, spNet_micro, pol.sp.micro,
     file = file.path(save.dir, 'microNets.RData'))


## calculate network summary statistics

indivNetSums_micro  <-  netSums(indivNet_micro,
                                binary=FALSE,
                                spec=spec.net,
                                FUN=mclapply,
                                m=ndim.func)

indivNetSums_microBinary  <-  netSums(indivNet_micro,
                                      binary=TRUE,
                                      spec=spec.net,
                                      FUN=mclapply,
                                      m=ndim.func)

# spNetSums_micro  <-  netSums(spNet_micro, binary=FALSE,
#                              spec=spec.net,
#                              FUN=mclapply,
#                              type="sp")
# 
# spNetSums_microBinary  <-  netSums(spNet_micro, binary=TRUE,
#                                    spec=spec.net,
#                                    FUN=mclapply,
#                                    type="sp")


save(indivNetSums_micro, indivNetSums_microBinary,
     #spNetSums_micro, spNetSums_microBinary,
     file =file.path(save.dir, 'NetSums_micro.RData'))

## ***********************************************************************
## Parasite networks
## ***********************************************************************


# parasite networks not working right now!
# 
# 
# #get columns of parasite presence to feed into illumSplit
# 
# names.Para <- c("AscosphaeraSpp","ApicystisSpp", "AspergillusSpp",
#                 "CrithidiaBombi", "CrithidiaExpoeki",
#                 "NosemaCeranae", "NosemaBombi" )
# 
# ## no       "CrithidiaSpp",
# 
# spec.net[, names.Para] <- apply(spec.net[, names.Para], 2, as.numeric)
# seq.para <- spec.net[apply(spec.net[, names.Para], 1, function(x) !all(is.na(x))),]
# 
# indivNet_para  <-  illumSplit(seq.para,"Site",names.Para)
# 
# 
# pol.sp.para <- lapply(indivNet_para,
#                 function(x){
#                     spec.net$GenusSpecies[match(colnames(x),
#                                             spec.net$UniqueID)]
# })
# 
# ## species level
# pol.abund.site.para <- pol.abund.site[names(indivNet_para)]
# 
# ## species level network
# spNet_para <-  mapply(makeSpNet,
#     indiv.nets = indivNet_para,
#     sp.names = pol.sp.para,
#     sp.abund=pol.abund.site.para,
#     SIMPLIFY = FALSE)
# 
# 
# #save networks themselves
# save(indivNet_para, spNet_para, pol.abund.site.para,
#      file =file.path(save.dir, 'paraNets.RData'))
# 
# #calculate network summary statistics
# 
# indivNetSums_para  <-  netSums(indivNet_para, spec=spec.net)
# 
# # spNetSums_para  <-  netSums(spNet_para, spec=spec.net,
# #                             type="sp", FUN=lapply)
# 
# 
# save(indivNetSums_para, #spNetSums_para,
#      file = file.path(save.dir, 'NetSums_para.RData'))

## ***********************************************************************
## combined dataframe of all the network metrics
## ***************************************************************
load( file =file.path(save.dir, 'NetSums_micro.RData'))
load(   file=file.path(save.dir, 'NetSums_rbcl.RData'))
#load(   file=file.path(save.dir, 'NetSums_para.RData'))


#net.metrics <- colnames(indivNetSums_para)[!
#    colnames(indivNetSums_para) %in% c("Site", "GenusSpecies", "UniqueID")]

net.metrics <- colnames(indivNetSums_micro)[!
   colnames(indivNetSums_micro) %in% c("Site", "GenusSpecies", "UniqueID")]

makeNewColnames <- function(netSums, name){
    colnames(netSums)[colnames(netSums) %in% net.metrics] <-
        paste(name,  colnames(netSums)[colnames(netSums) %in% net.metrics],
              sep="_")
    return(colnames(netSums))
}

#colnames(indivNetSums_para) <- makeNewColnames(indivNetSums_para, "Parasite")
colnames(indivNetSums_rbcl) <- makeNewColnames(indivNetSums_rbcl, "RBCL")
colnames(indivNetSums_micro) <- makeNewColnames(indivNetSums_micro, "Micro")

#sp level not working rn
# colnames(spNetSums_para) <- makeNewColnames(spNetSums_para, "Parasite")
# colnames(spNetSums_rbcl) <- makeNewColnames(spNetSums_rbcl, "RBCL")
# colnames(spNetSums_micro) <- makeNewColnames(spNetSums_micro, "Micro")


## combine individual level metrics
all.indiv.mets <- merge(indivNetSums_micro, indivNetSums_rbcl, all=TRUE)
#all.indiv.mets <- merge(all.indiv.mets, indivNetSums_para, all=TRUE)

all.indiv.mets$Site <- spec.net$Site[match(all.indiv.mets$UniqueID,
                                       spec.net$UniqueID)]
all.indiv.mets$GenusSpecies <- spec.net$GenusSpecies[match(all.indiv.mets$UniqueID,
                                       spec.net$UniqueID)]

## combine species level metrics -- not working
# all.sp.mets <- merge(spNetSums_micro, spNetSums_rbcl, all=TRUE)
# all.sp.mets <- merge(all.sp.mets, spNetSums_para, all=TRUE)


## do the same thing for binary networks
colnames(indivNetSums_rbclBinary) <-
    makeNewColnames(indivNetSums_rbclBinary, "RBCL")
colnames(indivNetSums_microBinary) <- makeNewColnames(indivNetSums_microBinary,
                                                      "Micro")

#not working right now

# colnames(spNetSums_rbclBinary) <- makeNewColnames(spNetSums_rbclBinary,
#                                                   "RBCL")
# colnames(spNetSums_microBinary) <- makeNewColnames(spNetSums_microBinary,
#                                                    "Micro")


## combine individual level metrics
all.indiv.mets.binary <- merge(indivNetSums_microBinary,
                               indivNetSums_rbclBinary, all=TRUE)
#all.indiv.mets.binary <- merge(all.indiv.mets.binary,
#                               indivNetSums_para, all=TRUE)
all.indiv.mets.binary$Site <- spec.net$Site[match(all.indiv.mets.binary$UniqueID,
                                       spec.net$UniqueID)]
all.indiv.mets.binary$GenusSpecies <- spec.net$GenusSpecies[
                                       match(all.indiv.mets.binary$UniqueID,
                                             spec.net$UniqueID)]

all.indiv.mets.binary <- merge(all.indiv.mets.binary, spec.net)

## combine species level metrics
# all.sp.mets.binary <- merge(spNetSums_microBinary,
#                             spNetSums_rbclBinary,
#                             all=TRUE)
# all.sp.mets.binary <- merge(all.sp.mets.binary, spNetSums_para, all=TRUE)


## merge with specimen data at individual level
#load('../skyIslands/data/spec_traits.Rdata')

## names.Para <- c("Ascosphaera","Apicystis",
##                 "CrithidiaSpp", "CrithidiaBombi", "CrithidiaExpoeki",
##                 "NosemaCeranae", "NosemaBombi" )

#spec.net[, names.Para] <- NULL
all.indiv.mets <- merge(all.indiv.mets, spec.net)
## all.indiv.mets.binary <- merge(all.indiv.mets.binary, spec.net)


## species level merge
#sp.by.site <- read.csv("../skyIslands/data/spstats.csv")

#all.sp.mets <- merge(all.sp.mets, sp.by.site)
#all.sp.mets.binary <- merge(all.sp.mets.binary, sp.by.site)


save(all.sp.mets.binary, all.sp.mets,
     all.indiv.mets.binary,
     all.indiv.mets,
     file = file.path(save.dir, 'allNetSums.RData'))



