
plotVarMean <- function(metrics, level, mods, mean.pv, mets){
    panels <- vector("list", length(metrics))
    for(i in 1:length(metrics)){
        ee <- Effect("maxN", mods[[metrics[i]]])

        this.dat <- mets[mets$speciesType == level,]
        colnames(this.dat)[colnames(this.dat) ==
                           paste(metrics[i], mean.pv, sep=".")] <- "y"
        panels[[i]] <- ggplot(as.data.frame(ee),
                              aes(maxN, fit))+
            geom_line()+
            geom_point(aes(x=maxN, y=y), data=this.dat) +
            ## colour=NA suppresses edges of the ribbon
            geom_ribbon(colour=NA,alpha=0.1,
                        aes(ymin=lower,ymax=upper))+
            labs(x="Occurrence", y=metrics[i])
    }

    do.call(grid.arrange, panels)
}


plotGeoVarHigher <- function(){
    plotVarMean(metrics= metrics, level="higher.level",
                mods=mods.pv.higher, mean.pv="pv", mets=geo.met)
}


plotGeoVarLower <- function(){
    plotVarMean(metrics= metrics, level="lower.level",
                mods=mods.pv.lower, mean.pv="pv", mets=geo.met)
}


plotGeoMeanHigher <- function(){
    plotVarMean(metrics= metrics, level="higher.level",
                mods=mods.mean.higher, mean.pv="mean", mets=geo.met)
}


plotGeoMeanLower <- function(){
    plotVarMean(metrics= metrics, level="lower.level",
                mods=mods.mean.lower, mean.pv="mean", mets=geo.met)
}


plotAllGeo <- function(){
    pdf.f(plotGeoVarHigher,  file=sprintf("figures/GeoVar%s%s.pdf", net.type,
                                       paste(species[2], collapse="")),
          height=9, width=8)

    pdf.f(plotGeoVarLower,  file=sprintf("figures/GeoVar%s%s.pdf", net.type,
                                         paste(species[1], collapse="")),
          height=9, width=8)


    pdf.f(plotGeoMeanHigher,  file=sprintf("figures/GeoMean%s%s.pdf", net.type,
                                        paste(species[2], collapse="")),
          height=9, width=8)

    pdf.f(plotGeoMeanLower,  file=sprintf("figures/GeoMean%s%s.pdf", net.type,
                                          paste(species[1], collapse="")),
          height=9, width=8)

}



plotYrVarHigher <- function(){
    plotVarMean(metrics= metrics, level="higher.level",
                mods=mods.pv.higher, mean.pv="pv", mets=yr.met)
}


plotYrVarLower <- function(){
    plotVarMean(metrics= metrics, level="lower.level",
                mods=mods.pv.lower, mean.pv="pv", mets=yr.met)
}


plotYrMeanHigher <- function(){
    plotVarMean(metrics= metrics, level="higher.level",
                mods=mods.mean.higher, mean.pv="mean", mets=yr.met)
}


plotYrMeanLower <- function(){
    plotVarMean(metrics= metrics, level="lower.level",
                mods=mods.mean.lower, mean.pv="mean", mets=yr.met)
}


plotAllYr <- function(){
    pdf.f(plotYrVarHigher,  file=sprintf("figures/YrVar%s%s.pdf", net.type,
                                       paste(species[2], collapse="")),
          height=9, width=8)

    pdf.f(plotYrVarLower,  file=sprintf("figures/YrVar%s%s.pdf", net.type,
                                         paste(species[1], collapse="")),
          height=9, width=8)


    pdf.f(plotYrMeanHigher,  file=sprintf("figures/YrMean%s%s.pdf", net.type,
                                        paste(species[2], collapse="")),
          height=9, width=8)

    pdf.f(plotYrMeanLower,  file=sprintf("figures/YrMean%s%s.pdf", net.type,
                                          paste(species[1], collapse="")),
          height=9, width=8)

}
