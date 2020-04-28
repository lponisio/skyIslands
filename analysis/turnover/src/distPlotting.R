library(viridis)

plotDists <- function(){

    panels <- vector("list", length(yvars))

    plotGeoTurnover <- function(){
        for(i in 1:length(yvars)){
            ee <- geo.mods[[yvars[i]]]$eff

            this.beta <- beta.same.year
            colnames(this.beta)[colnames(this.beta) == yvars[i]] <-
                "y"
            panels[[i]] <- ggplot(as.data.frame(ee),
                                  aes(GeoDist, fit))+
                geom_line()+
                geom_point(aes(x=GeoDist, y=y), data=this.beta) +
            ## colour=NA suppresses edges of the ribbon
            geom_ribbon(colour=NA,alpha=0.1,
                        aes(ymin=lower,ymax=upper))+
                labs(x="Geographic Distance", y=ylabs[i])
        }

        do.call(grid.arrange, panels)
    }

    pdf.f(plotGeoTurnover,  file=sprintf("figures/GeoBeta%s%s.pdf", net.type,
                                         paste(species, collapse="")
                                         ),
          height=9, width=8)

    plotYrTurnover <- function(){
        for(i in 1:length(yvars)){
            this.beta <- beta.same.site
            colnames(this.beta)[colnames(this.beta) == yvars[i]] <- "y"

            panels[[i]] <- ggplot(this.beta) +
                geom_boxplot(aes(factor(YrDist), y=y)) +
                scale_x_discrete(labels = c("1", "5", "6")) +
                labs(x="Years between surveys", y=ylabs[i]) +
                lims(y=c(0,1))
        }

        do.call(grid.arrange, panels)
    }

    pdf.f(plotYrTurnover,  file=sprintf("figures/YrBeta%s%s.pdf", net.type,
                                        paste(species, collapse="")
                                        ),
          height=9, width=8)

    if(nets.by.SR){
        plotSRTurnover <- function(){
            for(i in 1:length(yvars)){
                this.beta <- beta.same.site.year
                colnames(this.beta)[colnames(this.beta) == yvars[i]] <- "y"
                panels[[i]] <- ggplot(this.beta,
                                      aes(x=SRDist, y=y)) +
                    geom_point() + geom_smooth(method="lm") +
                    labs(x="Days", y=ylabs[i]) +
                    lims(y=c(0,1))
            }

            do.call(grid.arrange, panels)
        }

        pdf.f(plotSRTurnover,  file=sprintf("figures/SRBeta%s%s.pdf", net.type,
                                            paste(species, collapse="")
                                            ),
              height=9, width=8)
    }
}

plotPCAs <- function(){
    panels <- vector("list", length(xvars))
    plotPCATurnover <- function(){
        for(i in 1:length(xvars)){
            ee <- pca.mods[[xvars[i]]]$eff
            this.beta <-  pcas.beta
            colnames(this.beta)[colnames(this.beta) == xvars[i]] <-
                "x"
            panels[[i]] <- ggplot(as.data.frame(ee),
                                  aes(x, fit))+
                geom_line()+
                geom_point(aes(x=x, y=diffPca,
                               color=GenusSpecies),
                           data=this.beta) +
                 scale_color_viridis(discrete = TRUE, option = "D") +
                geom_ribbon(colour=NA,alpha=0.1,
                            aes(ymin=lower,ymax=upper))+
                labs(y="Change in network position", x=xlabs[i])  +
                theme(legend.position = "none")

        }

        do.call(grid.arrange, panels)
    }
    pdf.f(plotPCATurnover,  file=sprintf("figures/pcaDiff%s%s%s.pdf",
                                         net.type,
                                         paste(species, collapse=""),
                                         gsub("[.]", "", sp.level)
                                         ),
          height=8, width=8)
}
