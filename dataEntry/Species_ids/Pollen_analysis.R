source('~/Documents/Hedgerow Networks/Code/Functions/samp2site_spp.R', chdir = TRUE)
source('~/Documents/Hedgerow Networks/Code/Final_net_null_with_ISA_H2_div.R', chdir = TRUE)
bfly <- read.csv("~/Dropbox/SI_data_entry/butterfly_labels3.csv")
pollen <- read.csv("~/Documents/Mentoring/Becky Wong/Data.csv")

############# pollen data prep #############

## merge the pollen data by species so each row is a species and each column is a plant species

bysp.pollen <- sapply(3:90, function(i) tapply(pollen[,i], pollen$gen_sp, sum))
colnames(bysp.pollen) <- colnames(pollen[3:90])

abund.bysp.pollen <- bysp.pollen ##keep track of the abundance for later

bysp.pollen[bysp.pollen > 0] <- 1 ## for richness replace all counts with 1

############# butterfly data prep #############
##merge the observation data by species so each row is a species and each column is a plant species
ag.bfly <- aggregate(bfly$finalPlantID, list(gen_sp = bfly$genus_species, plant_sp = bfly$finalPlantID), length)

abund.ag.bfly <- ag.bfly ##keep track of the abundance for later

bysp.obs.abund <- samp2site.spp(abund.ag.bfly$gen_sp, abund.ag.bfly$plant_sp, abund.ag.bfly$x)

ag.bfly$x[ag.bfly$x > 0] <- 1 ## for richness replace all counts with 1

bysp.obs <- samp2site.spp(ag.bfly$gen_sp, ag.bfly$plant_sp, ag.bfly$x)


## calculate the richness for each butterfly species
pollen.rich <- apply(bysp.pollen,1,sum)
obs.rich <- apply(bysp.obs,1,sum)

##take the subset of species that are in both obs data and pollen data 
obs.rich.abrev <- obs.rich[names(obs.rich) %in% names(pollen.rich)]
pollen.rich.abrev <- pollen.rich[names(pollen.rich) %in% names(obs.rich)]

## paried t-test on richness 
t.test(log(obs.rich.abrev), log(pollen.rich.abrev), paired=TRUE)

#	Paired t-test

#data:  log(obs.rich.abrev) and log(pollen.rich.abrev) 
#t = -7.0421, df = 12, p-value = 1.352e-05
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -1.236735 -0.652277 
#sample estimates:
#mean of the differences 
 #            -0.9445061 


boxplot(obs.rich.abrev, pollen.rich.abrev, names=c("Observation","Pollen"), col= c("darkolivegreen", "magenta4"), log="y", ylab="ln Richness", cex.lab=1.3, cex.axis=1.5)

quartz(width=5,height=5)
par(xpd=NA, mar=c(8.5,5,0,1))

order.pol <- order(pollen.rich.abrev)
plot(obs.rich.abrev[order.pol], col="darkolivegreen", pch=16, cex=2, ylim=c(0,25), xaxt="n", xlab="", ylab="Richness")

labs <- c("Poanes zabulon", "Vanessa atlanta","Phyciodes tharos","Pontia protodice","Limenitis bredowii","Papilio polyxenes", "Polygonia gracilis","Vanessa cardui","Eurema nicippe","Vanessa carye","Speyeria atlantis","Euptoieta Claudia","Colias eurytheme")
axis(side=1, at=1:13, labels=FALSE)
text(1:13, y=rep(-2.5,13), labels=labs, srt=45, adj=1)

points(pollen.rich.abrev[order.pol], col="magenta4", pch=16, cex=2)
legend("topleft", legend=c("Observation","Pollen"),pch=16, cex=1.5, col= c("darkolivegreen", "magenta4"), bty="n")

#####specialization analysis#######
library(bipartite) 

pollen.d <- specieslevel(abund.bysp.pollen, index="d")
obs.d <- specieslevel(bysp.obs.abund, index="d")

##take the subset of species that are in both obs data and pollen data 
obs.d.abrev <- obs.d[[2]][rownames(obs.d$'lower trophic level') %in% rownames(pollen.d$'lower trophic level'),]
pollen.d.abrev <- pollen.d[[2]][rownames(pollen.d$'lower trophic level') %in% rownames(obs.d$'lower trophic level'),]

t.test(obs.d.abrev, pollen.d.abrev, paired=TRUE)

boxplot(obs.d.abrev, pollen.d.abrev, names=c("Observation","Pollen"), col= c("darkolivegreen", "magenta4"), ylab="Specialization", cex.lab=1.3, cex.axis=1.5)


quartz(width=5,height=5)

par(xpd=NA, mar=c(8.5,5,0,1))
order.pol <- order(pollen.d.abrev)
plot(obs.d.abrev[order.pol], col="darkolivegreen", pch=16, cex=2, ylim=c(0,1), xlab="", ylab="Specialization", xaxt="n")

axis(side=1, at=1:13, labels=FALSE)
text(1:13, y=rep(-0.1,13), labels=labs, srt=45, adj=1)

points(pollen.d.abrev[order.pol], col="magenta4", pch=16, cex=2)
legend("topleft", legend=c("Observation","Pollen"),pch=16, cex=1.5, col= c("darkolivegreen", "magenta4"), bty="n")

#	Paired t-test

#data:  obs.d.abrev and pollen.d.abrev 
#t = 0.5591, df = 12, p-value = 0.5864
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
# -0.1095336  0.1851532 
#sample estimates:
#mean of the differences 
 #            0.03780979 


###########network fig
library(bipartite)
bysp.pollen.abund.abrev <- abund.bysp.pollen[rownames(abund.bysp.pollen) %in% rownames(bysp.obs.abund),]
bysp.pollen.abund.abrev <- empty(bysp.pollen.abund.abrev)
pollen.stats <- network.metrics(bysp.pollen.abund.abrev, 10)

bysp.obs.abund.abrev <- bysp.obs.abund[rownames(bysp.obs.abund) %in% rownames(abund.bysp.pollen),]
bysp.obs.abund.abrev <- empty(bysp.obs.abund.abrev)
obs.stats <- network.metrics(bysp.obs.abund.abrev, 10)

library(sna)
	gplot(bysp.pollen.abund.abrev,
			gmode="twomode",
			label= c(rownames(bysp.pollen.abund.abrev),colnames(bysp.pollen.abund.abrev)),
			vertex.col=c(rep("goldenrod1",nrow(bysp.pollen.abund.abrev)), rep("forestgreen",ncol(bysp.pollen.abund.abrev))),
			label.cex = 0.4,
			usearrows="FALSE",
			mode="fruchtermanreingold",
			main="Pollen",
			object.scale = 0.03)
			
			
gplot(bysp.obs.abund.abrev,
			gmode="twomode",
			label= c(rownames(bysp.obs.abund.abrev),colnames(bysp.obs.abund.abrev)),
			vertex.col=c(rep("goldenrod1",nrow(bysp.obs.abund.abrev)), rep("forestgreen",ncol(bysp.obs.abund.abrev))),
			label.cex = 0.4,
			usearrows="FALSE",
			mode="fruchtermanreingold",
			main="Observation",
			object.scale = 0.03)
			