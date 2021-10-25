setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/multilayer')

library(ggplot2)
library(igraph)

load('../../data/spec.Rdata')
load('../../data/nets.Rdata')

same.year <- split(nets.graph.uw, years)

## **************************************************************
## 2012 multi layer network
## **************************************************************

multiplex.network <- list()
multiplex.network[[1]] <-  same.year$'2012'$JC.2012.1
multiplex.network[[2]] <-  same.year$'2012'$SC.2012.1
multiplex.network[[3]] <-  same.year$'2012'$MM.2012.1
multiplex.network[[4]] <-  same.year$'2012'$PL.2012.1
multiplex.network[[5]] <-  same.year$'2012'$CH.2012.1
names(multiplex.network) <- c("JC", "SC", "MM", "PL", "CH")

alpha.vals <- seq(from=0,to=1,by=0.25) ## distinct values of alpha
alpha.vals <- alpha.vals[alpha.vals!=0]
alpha.vals <- c(alpha.vals,seq(from=2,to=5,by=1))
l <- length(alpha.vals)
## number of distinct values

scale <- 20									# node scale for graph plots
formats <- c(								# file format of the plots
    "PDF",
    "PNG"
)




number.layers <- length(multiplex.network)
number.nodes <- length(unique(
    unlist(sapply(multiplex.network,function (x) names(V(x))))))

                                        #setup personal opinion
personal.opinion <- array(0.9, c(number.layers, number.nodes))

                                        # init centrality tables

                                        # process our centrality measure
cat("  Processing the opinion centrality\n",sep="")
for(i in 1:l)
{	cat("    for alpha=",alpha.vals[i]," (",i,"/",l,")\n",sep="")
    alpha <- matrix(alpha.vals[i],nrow=number.nodes, ncol=number.layers)
                                        #print(alpha)

####### process opinion centrality measure
    elapsed.time <- system.time(
        centrality <- process.opinion.centrality(network=multiplex.network, alpha, budget=1/(number.nodes)^3, personal.opinion)
    )
                                        # normalize
    centrality <- centrality - min(centrality)
    centrality <- centrality / sum(centrality)

    elapsed.times[network.name,i] <- elapsed.time["elapsed"]
    cat("Elapsed times:\n");print(elapsed.times)
    write.csv2(elapsed.times, file=time.data.file)
                                        #elapsed.times <- read.csv2(file=time.data.file,header=TRUE,row.names=1,check.names=FALSE)
                                        #net.prop <- read.csv2(file=netprop.file,header=TRUE,row.names=1,check.names=FALSE)
    plot.time.perf(elapsed.times, net.prop, plot.file=time.plot.file, dispersion=FALSE, formats)

    opinion.centralities[i,] <- t(centrality)
    stdev <- sd(opinion.centralities[i,])
    if(stdev==0)
        cat("....WARNING: stdev=",stdev," (value=",opinion.centralities[i,1],")\n",sep="")
    else
        cat("....stdev=",stdev,"\n",sep="")
}

                                        # record our measure as a table
out.file <- paste(net.plot.folder,"opinion-centrality.csv",sep="")
col.node <- 1:number.nodes
centr <- cbind(col.node, t(opinion.centralities))
colnames(centr) <- c("Node", paste("alpha=",alpha.vals,sep=""))
write.csv2(centr, file=out.file, row.names=FALSE)

                                        # compare each measure to our own
cat("  Compare the opinion measure to the others\n")
                                        # we consider all possible values of alpha
for(i in 1:l)
{	cat("    Processing alpha=",alpha.vals[i],"\n",sep="")

                                        # produce the opinion centrality histogram for the considered value of alpha
    cat("      Generate histogram for the opinion centrality with alpha=",alpha.vals[i],"\n",sep="")
    measure.histo(vals=opinion.centralities[i,], alpha=alpha.vals[i], folder=net.plot.folder, formats=formats)

                                        # process each alternate measure individually
    for(measure in measures)
    {	# check if the considered measure was processed
        if(any(is.na(other.centralities[,measure])))
            cat("      Measure ",measure," was not processed for this dataset\n",sep="")
        else
        {	# process the rank correlation with our measure
            correlation.values[i,measure] <- cor(opinion.centralities[i,], other.centralities[,measure], method="spearman")

                                        # plot ranking differences
            cat("      Generate line plot representing ranking differences with measure ",measure,"\n",sep="")
            rank.diff.lineplot(ref.vals=other.centralities[,measure], comp.vals=opinion.centralities[i,], ref.measure=measure, alpha=alpha.vals[i], folder=net.plot.folder, formats=formats)

                                        # ranking differences as a barplot
            cat("      Generate barplot representing ranking differences with measure ",measure,"\n",sep="")
            rank.diff.barplot(ref.vals=other.centralities[,measure], comp.vals=opinion.centralities[i,], ref.measure=measure, alpha=alpha.vals[i], folder=net.plot.folder, formats=formats)

                                        # store rank difference for k most central nodes (according to alt. measure)
            if(i==1) # just for the first alpha value, since it does not affect ranks
            {	ref.rk <- rank(other.centralities[,measure],ties.method="min")
                comp.rk <- rank(opinion.centralities[i,],ties.method="min")
                diff <- comp.rk - ref.rk
					#print(diff)
                idx <- order(other.centralities[,measure], decreasing=TRUE)[1:knodes]
					#print(idx)
                all.rank.diff[[measure]][network.name,1:knodes] <- diff[idx]
            }

                                        # plot the aggregated network with each existing measure as the size, and the opinion measure as the color
            cat("      Generate a plot representing the (aggregated) graph and measure ",measure,"\n",sep="")
            graph.plot(g=aggregated.network, ref.vals=other.centralities[,measure], comp.vals=opinion.centralities[i,], ref.measure=measure, alpha=alpha.vals[i], folder=net.plot.folder, layout=lay, scale=scale, formats=formats)
        }
    }
}

                                        # generate plots comparing existing centrality measures
cat("  Generate plots comparing existing centrality measures\n")
for(m1 in 1:(length(measures)-1))
{	for(m2 in (m1+1):length(measures))
            rank.diff.barplot(ref.vals=other.centralities[,measures[m1]], comp.vals=other.centralities[,measures[m2]], ref.measure=measures[m1], comp.measure=measures[m2], alpha=NA, folder=net.plot.folder, formats=formats)
}

                                        # plot the correlation between our measure and the other ones
cat("  Plot the correlations between the opinion measure and the other measures\n")
for(measure in measures)
{	# check if the considered measure was processed
    if(any(is.na(other.centralities[,measure])) | length(unique(other.centralities[,measure]))==1)
        cat("    WARNING: could not find the values for measure ",measure,", so no correlation plot\n",sep="")

                                        # check if our own centrality could be processed
    else if(all(is.na(correlation.values[,measure])))
        cat("    WARNING: Opinion centrality could not be processed, so no correlation plot for ",measure,"\n",sep="")
    else
    {	cat("    With measure ",measure,"\n",sep="")
        corr.plot(cor.vals=correlation.values[,measure], alpha.vals, measure, folder=net.plot.folder, formats=formats)
    }
}

                                        # record the correlations between our measure and the other ones
corr.folder <- paste(net.plot.folder,"/corr_plots/",sep="")
dir.create(corr.folder,showWarnings=FALSE,recursive=TRUE)
out.file <- paste(corr.folder,"all-correlations.csv",sep="")
write.csv2(correlation.values, file=out.file)

                                        # draw all measure correlations in the same plot
if(!all(is.na(correlation.values)))
    corr.plot.all(cor.vals=correlation.values, alpha.vals, measures, folder=net.plot.folder, formats=formats)

                                        # complete and record the overall correlation table
all.corr.table[network.name,] <- apply(correlation.values,2,mean)
write.csv2(all.corr.table, file=all.corr.file)
cat("Overall correlation table:\n");print(all.corr.table)

                                        # process and record correlation matrix for opinion measure only
cat("  Record the correlations for the opinion measure only (in function of alpha)\n")
opinion.correlation <- matrix(NA,nrow=l,ncol=l)
colnames(opinion.correlation) <- alpha.vals
rownames(opinion.correlation) <- alpha.vals
for(i in 1:l)
{	for(j in 1:l)
            opinion.correlation[i,j] <- cor(opinion.centralities[i,], opinion.centralities[j,], method="spearman")
}
correlation.plot(corr.mat=opinion.correlation, folder=net.plot.folder, formats=formats)
}

