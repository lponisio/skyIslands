rm(list=ls())
library(vegan)
library(fields)
library(lme4)
setwd("~/Dropbox/SkyIslands/analysis/distanceDecay")
source('src/misc.R')

spec <-
  read.csv("../data/spec.csv")


##************************************************************************
## d' by occurance
##************************************************************************

prep.comm <- aggregate(spec$GenSp,
                       list(site=spec$PlantGenSp,
                            sp=spec$GenSp), length)
comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x)

d <- specieslevel(comm, index="d")

nest <- specieslevel(comm, index="nestedrank")

betw <- specieslevel(comm, index="betweenness")

occ <- function(type, colname, famcol, metric, d.scores){
  bysite <- aggregate(spec$Site, list(sp=spec[,colname],
                                      site=spec$Site), length)
  sp.bysite <- samp2site.spp(bysite$site, bysite$sp, bysite$x)
  sp.bysite[sp.bysite > 1] <- 1
  occ.sp <- apply(sp.bysite, 2, sum)
  glm.sp <- glm(cbind(occ.sp, 5-occ.sp) ~ d.scores,
                family="binomial")
  
  fams <- spec[,famcol][match(names(occ.sp), spec[,colname])]
  if(length(unique(fams)) <= 12){
    cols <- brewer.pal(length(unique(fams)), "Paired")
  }else{
    cols <- rainbow(length(unique(fams)))
  }
  names(cols) <- unique(fams)
  all.cols <- cols[match(fams, names(cols))]
  
  f <- function(){
    boxplot(d.scores ~occ.sp, ylab=metric,
            xlab="Occurance")
  }
  path <- '~/Dropbox/SkyIslands/analysis/figures/occurence' 
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste(type, "occ", metric, "box", sep="_"))),
        width=5, height=5)

  f2 <- function(){
    plot(occ.sp ~ d.scores, xlab=metric,
         ylab="Occurence", pch=16, col=all.cols, xlim=range(d.scores),
         ylim=c(1,5))
    #abline(coef=glm.sp$coeff)
    legend("topright", legend=unique(fams), col=cols, pch=16, cex=0.5)
  }
  pdf.f(f2, file= file.path(path, sprintf("%s.pdf",
              paste(type, "occ", metric, "plot", sep="_"))),
        width=5, height=5)
  return(summary(glm.sp))
}

occ.pol <- occ(type="pol", colname="GenSp", d.scores=
               d$'higher level'$d, famcol="Family", metric="Specialization")

occ.plant <- occ(type="plant", colname="PlantGenSp",
                 d.scores= d$'lower level'$d, famcol="PlantFamily",
                 metric="Specialization")

## occ.pol.nest <- occ(type="pol", colname="GenSp", d.scores=
##                nest$'higher level'$nestedrank, "Family", "Nested Rank")

## occ.plant.nest <- occ(type="plant", colname="PlantGenSp",
##                  d.scores= nest$'lower level'$nestedrank,
##                  "PlantFamily", "Nested Rank")

## occ.pol.betw <- occ(type="pol", colname="GenSp", d.scores=
##                betw$'higher level'$weighted.betweeness, "Family", "Betweenness")

## occ.plant.betw <- occ(type="plant", colname="PlantGenSp",
##                  d.scores= betw$'lower level'$weighted.betweeness,
##                  "PlantFamily", "Betweenness")






