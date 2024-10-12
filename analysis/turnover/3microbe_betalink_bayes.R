## setwd("~/Dropbox/skyIslands/")
rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)

setwd("skyIslands/analysis/turnover")
#source("src/initialize.R")
source("src/chao.R")
source("src/betaNet.R")
source("src/writeResultsTable.R")
library(ggplot2)
library(lme4)
library(lmerTest)
library(igraph)
library(ggpubr)
library(emmeans)
library(bipartite)
library(dplyr)
library(fields)
library(brms)
library(tidybayes)
library(gridExtra)
library(bayesplot)
library(glmmTMB)
library(performance)

load("../../data/networks/microNets.RData")

load("../../data/spec_RBCL_16s.RData")

source("src/writeResultsTable.R")
source("src/networkTurnover.R")
##############################################################################


## Trying for just one to test code

#if(obligate.microbe.network==TRUE){

## now do the same for the three groups of obligate, parasitic, transient


## WIP


## prep obligate network
## error  object 'CH' not found
only_obligate_network <- prep_obligate_network(raw_network=spNet_micro)

CH <- only_obligate_network$CH
CH <- CH[,colnames(CH)!=""]
HM <- only_obligate_network$HM
JC <- only_obligate_network$JC
MM <- only_obligate_network$MM 
MM <- MM[,colnames(MM)!=""]
PL <- only_obligate_network$PL
PL <- PL[,colnames(PL)!=""]
RP <- only_obligate_network$RP
SC <- only_obligate_network$SC 
SM <- only_obligate_network$SM


lower.order <- "Microbes"
higher.order <- "Pollinators"


obligate_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, RP, SM, SC),
                                          partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)

#View(obligate_poll_betalink)

colnames(obligate_poll_betalink) <- c("Site1",
                                      "Site2",
                                      "DissimilaritySpeciesComposition",
                                      "OnlySharedLinks",
                                      "WholeNetworkLinks",
                                      "SpeciesTurnoverLinks",
                                      paste("TurnoverAbsence",lower.order,sep=""),
                                      paste("TurnoverAbsence",higher.order,sep=""),
                                      "TurnoverAbsenceBoth")


geo <- unique(spec.net[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site

## add column for geographic distance between sites
obligate_poll_betalink$GeoDist <- apply(obligate_poll_betalink, 1, function(x){
  geo.dist[x["Site1"],  x["Site2"]]
})

dir.create("figures", showWarnings = FALSE)
dir.create("figures/obligate_microbe_poll", showWarnings = FALSE)

View(obligate_poll_betalink)

## prep transient now
only_transient_network <- prep_transient_network(raw_network=spNet_micro)

#### species level networks
CH <- only_transient_network$CH
CH <- CH[,colnames(CH)!=""]
HM <- only_transient_network$HM
JC <- only_transient_network$JC
MM <- only_transient_network$MM 
MM <- MM[,colnames(MM)!=""]
PL <- only_transient_network$PL
PL <- PL[,colnames(PL)!=""]
RP <- only_transient_network$RP
SC <- only_transient_network$SC 
SM <- only_transient_network$SM


lower.order <- "Microbes"
higher.order <- "Pollinators"


transient_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, RP, SM, SC),
                                           partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)

#View(transient_poll_betalink)

colnames(transient_poll_betalink) <- c("Site1",
                                       "Site2",
                                       "DissimilaritySpeciesComposition",
                                       "OnlySharedLinks",
                                       "WholeNetworkLinks",
                                       "SpeciesTurnoverLinks",
                                       paste("TurnoverAbsence",lower.order,sep=""),
                                       paste("TurnoverAbsence",higher.order,sep=""),
                                       "TurnoverAbsenceBoth")


geo <- unique(spec.net[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site

## add column for geographic distance between sites
transient_poll_betalink$GeoDist <- apply(transient_poll_betalink, 1, function(x){
  geo.dist[x["Site1"],  x["Site2"]]
})

dir.create("figures", showWarnings = FALSE)
dir.create("figures/transient_microbe_poll", showWarnings = FALSE)

View(transient_poll_betalink)
## A. rewiring

rewiring.obligate.mod <- run_network_turnover_mod(this_component="OnlySharedLinks",
                                                          this_network=obligate_poll_betalink)

rewiring.transient.mod <- run_network_turnover_mod(this_component="OnlySharedLinks",
                                                          this_network=transient_poll_betalink)


## this is working, need to make function that plots both rewiring together
ob.rewiring.plot <- plot_network_turnover_mod_single(mod1=rewiring.obligate.mod, 
                                                          this.network=obligate_poll_betalink,
                                                          network_type='Obligate',
                                                          this.effect="GeoDist",
                                                          this.resp="OnlySharedLinks",
                                                          label="Rewiring"
                                                          )

trans.rewiring.plot <- plot_network_turnover_mod_single(mod1=rewiring.transient.mod, 
                                                     this.network=transient_poll_betalink,
                                                     network_type='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="OnlySharedLinks",
                                                     label="Rewiring"
)

### working up to here, now need to make function to add both plots together!

### WIP





species.turnover <- calculate_and_plot_betalinkr("DissimilaritySpeciesComposition",
                                                 obligate_poll_betalink,
                                                 "Dissimilarity: \nSpecies Composition",
                                                 network_type="Obligate")
ggsave(species.turnover[[2]], file="figures/obligate_microbe_poll/DissimilaritySpeciesTurnover.pdf", height=6, width=6)

interaction.turnover <- calculate_and_plot_betalinkr("WholeNetworkLinks",
                                                     obligate_poll_betalink,
                                                     "Dissimilarity: \nInteraction Turnover",
                                                     network_type="Obligate")
ggsave(interaction.turnover[[2]], file="figures/obligate_microbe_poll/DissimilarityInteractionTurnover.pdf", height=6, width=6)

int.turnover.rewiring <- calculate_and_plot_betalinkr("OnlySharedLinks",
                                                      obligate_poll_betalink,
                                                      "Rewiring",
                                                      network_type="Obligate")
ggsave(int.turnover.rewiring[[2]], file="figures/obligate_microbe_poll/InteractionDissimilarityRewiring.pdf", height=6, width=6)

int.turnover.species.turnover <- calculate_and_plot_betalinkr("SpeciesTurnoverLinks",
                                                              obligate_poll_betalink,
                                                              "Interaction Turnover: \nSpecies Turnover",
                                                              network_type="Obligate")
ggsave(int.turnover.species.turnover[[2]], file="figures/obligate_microbe_poll/InteractionTurnoverSpeciesComp.pdf", height=6, width=6)

sp.turnover.microbes <- calculate_and_plot_betalinkr("TurnoverAbsenceMicrobes",
                                                     obligate_poll_betalink,
                                                     "Microbe-driven Turnover",
                                                     network_type="Obligate")
ggsave(sp.turnover.microbes[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceMicrobes.pdf", height=6, width=6)

sp.turnover.bees <- calculate_and_plot_betalinkr("TurnoverAbsencePollinators",
                                                 obligate_poll_betalink,
                                                 "Host-driven Turnover",
                                                 network_type="Obligate")
ggsave(sp.turnover.bees[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceBees.pdf", height=6, width=6)

sp.turnover.both <- calculate_and_plot_betalinkr("TurnoverAbsenceBoth",
                                                 obligate_poll_betalink,
                                                 "Complete Turnover",
                                                 network_type="Obligate")
ggsave(sp.turnover.both[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceBoth.pdf", height=6, width=6)

}

## make panel figure
# panelA <- species.turnover[[2]] + labs(tag="A.")
# panelB <- interaction.turnover[[2]] + labs(tag="B.")
# panelC <- sp.turnover.bees[[2]] + labs(tag="C.")
# panelD <- int.turnover.species.turnover[[2]] + labs(tag="D.")
# panelE <- sp.turnover.microbes[[2]] + labs(tag="E.")
# panelF <- int.turnover.rewiring[[2]] + labs(tag="F.")

## updated panel fig
## rewiring
dir.create("tables/obligate_microbe_poll", showWarnings = FALSE)

panelA <- int.turnover.rewiring[[2]] + labs(tag="A.")
write.table(int.turnover.rewiring[[1]],
            file=sprintf("tables/obligate_microbe_poll/rewiring.txt"), sep="&")
write.csv(int.turnover.rewiring[[1]],
          file=sprintf("tables/obligate_microbe_poll/rewiring.csv"))

panelB <- sp.turnover.bees[[2]] + labs(tag="B.")
write.table(sp.turnover.bees[[1]],
            file=sprintf("tables/obligate_microbe_poll/host-driven.txt"), sep="&")
write.csv(sp.turnover.bees[[1]],
          file=sprintf("tables/obligate_microbe_poll/host-driven.csv"))

panelC <- sp.turnover.microbes[[2]] + labs(tag="C.")
write.table(sp.turnover.microbes[[1]],
            file=sprintf("tables/obligate_microbe_poll/microbe-driven.txt"), sep="&")
write.csv(sp.turnover.microbes[[1]],
          file=sprintf("tables/obligate_microbe_poll/microbe-driven.csv"))

panelD <- sp.turnover.both[[2]] + labs(tag="D.")
write.table(sp.turnover.both[[1]],
            file=sprintf("tables/obligate_microbe_poll/complete-turnover.txt"), sep="&")
write.csv(sp.turnover.both[[1]],
          file=sprintf("tables/obligate_microbe_poll/complete-turnover.csv"))

pdf("figures/obligate_microbe_poll/grid_obligate_microbes.pdf", width = 8.5, height = 8.5) # Open a new pdf file
grid.arrange(panelA,
             panelB,
             panelC,
             panelD,
             ncol=2) # Write the grid.arrange in the file
dev.off()


if(transient.microbe.network==TRUE){

## now do for non obligate
dir.create("figures", showWarnings = FALSE)
dir.create("figures/transient_microbe_poll", showWarnings = FALSE)

site_list <- names(spNet_micro)

## obligate symbionts
bee.obligates <- "Lactobacillus|Bifidobacterium|Snodgrassella|Gilliamella|Frischella|Bartonella|Commensalibacter"



only_transient_network <- list()

for (x in site_list){

  trans_rows <- rownames(spNet_micro[[x]])
  
  trans_rows_to_keep <- !grep(bee.obligates, trans_rows)
  
  trans_new_net <- spNet_micro[[x]][!trans_rows_to_keep,]
  
  new_name <- x
  
  only_transient_network[[new_name]] <- trans_new_net
  
}

#### species level networks
CH <- only_transient_network$CH
CH <- CH[,colnames(CH)!=""]
HM <- only_transient_network$HM
JC <- only_transient_network$JC
MM <- only_transient_network$MM 
MM <- MM[,colnames(MM)!=""]
PL <- only_transient_network$PL
PL <- PL[,colnames(PL)!=""]
RP <- only_transient_network$RP
SC <- only_transient_network$SC 
SM <- only_transient_network$SM


lower.order <- "Microbes"
higher.order <- "Pollinators"


transient_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, RP, SM, SC),
                                          partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)

#View(transient_poll_betalink)

colnames(transient_poll_betalink) <- c("Site1",
                                      "Site2",
                                      "DissimilaritySpeciesComposition",
                                      "OnlySharedLinks",
                                      "WholeNetworkLinks",
                                      "SpeciesTurnoverLinks",
                                      paste("TurnoverAbsence",lower.order,sep=""),
                                      paste("TurnoverAbsence",higher.order,sep=""),
                                      "TurnoverAbsenceBoth")


geo <- unique(spec.net[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site

## add column for geographic distance between sites
transient_poll_betalink$GeoDist <- apply(transient_poll_betalink, 1, function(x){
  geo.dist[x["Site1"],  x["Site2"]]
})

dir.create("figures", showWarnings = FALSE)
dir.create("figures/transient_microbe_poll", showWarnings = FALSE)

species.turnover <- calculate_and_plot_betalinkr("DissimilaritySpeciesComposition",
                                                 transient_poll_betalink,
                                                 "Dissimilarity: \nSpecies Composition",
                                                 network_type = 'Transient')
ggsave(species.turnover[[2]], file="figures/transient_microbe_poll/DissimilaritySpeciesTurnover.pdf", height=6, width=6)

interaction.turnover <- calculate_and_plot_betalinkr("WholeNetworkLinks",
                                                     transient_poll_betalink,
                                                     "Dissimilarity: \nInteraction Turnover",
                                                     network_type = 'Transient')
ggsave(interaction.turnover[[2]], file="figures/transient_microbe_poll/DissimilarityInteractionTurnover.pdf", height=6, width=6)

int.turnover.rewiring <- calculate_and_plot_betalinkr("OnlySharedLinks",
                                                      transient_poll_betalink,
                                                      "Rewiring",
                                                      network_type = 'Transient')
ggsave(int.turnover.rewiring[[2]], file="figures/transient_microbe_poll/InteractionDissimilarityRewiring.pdf", height=6, width=6)

int.turnover.species.turnover <- calculate_and_plot_betalinkr("SpeciesTurnoverLinks",
                                                              transient_poll_betalink,
                                                              "Interaction Turnover: \nSpecies Turnover",
                                                              network_type = 'Transient')
ggsave(int.turnover.species.turnover[[2]], file="figures/transient_microbe_poll/InteractionTurnoverSpeciesComp.pdf", height=6, width=6)

sp.turnover.microbes <- calculate_and_plot_betalinkr("TurnoverAbsenceMicrobes",
                                                     transient_poll_betalink,
                                                     "Microbe-driven Turnover",
                                                     network_type = 'Transient')
ggsave(sp.turnover.microbes[[2]], file="figures/transient_microbe_poll/SpeciesTurnoverAbsenceMicrobes.pdf", height=6, width=6)

sp.turnover.bees <- calculate_and_plot_betalinkr("TurnoverAbsencePollinators",
                                                 transient_poll_betalink,
                                                 "Host-driven Turnover",
                                                 network_type = 'Transient')
ggsave(sp.turnover.bees[[2]], file="figures/transient_microbe_poll/SpeciesTurnoverAbsenceBees.pdf", height=6, width=6)

sp.turnover.both <- calculate_and_plot_betalinkr("TurnoverAbsenceBoth",
                                                 transient_poll_betalink,
                                                 "Complete Turnover",
                                                 network_type = 'Transient')
ggsave(sp.turnover.both[[2]], file="figures/transient_microbe_poll/SpeciesTurnoverAbsenceBoth.pdf", height=6, width=6)
}

## make panel figure
# panelA <- species.turnover[[2]] + labs(tag="A.")
# panelB <- interaction.turnover[[2]] + labs(tag="B.")
# panelC <- sp.turnover.bees[[2]] + labs(tag="C.")
# panelD <- int.turnover.species.turnover[[2]] + labs(tag="D.")
# panelE <- sp.turnover.microbes[[2]] + labs(tag="E.")
# panelF <- int.turnover.rewiring[[2]] + labs(tag="F.")

dir.create("tables/transient_microbe_poll", showWarnings = FALSE)

panelA <- int.turnover.rewiring[[2]] + labs(tag="A.")
write.table(int.turnover.rewiring[[1]],
            file=sprintf("tables/transient_microbe_poll/rewiring.txt"), sep="&")
write.csv(int.turnover.rewiring[[1]],
          file=sprintf("tables/transient_microbe_poll/rewiring.csv"))

panelB <- sp.turnover.bees[[2]] + labs(tag="B.")
write.table(sp.turnover.bees[[1]],
            file=sprintf("tables/transient_microbe_poll/host-driven.txt"), sep="&")
write.csv(sp.turnover.bees[[1]],
          file=sprintf("tables/transient_microbe_poll/host-driven.csv"))

panelC <- sp.turnover.microbes[[2]] + labs(tag="C.")
write.table(sp.turnover.microbes[[1]],
            file=sprintf("tables/transient_microbe_poll/microbe-driven.txt"), sep="&")
write.csv(sp.turnover.microbes[[1]],
          file=sprintf("tables/transient_microbe_poll/microbe-driven.csv"))

panelD <- sp.turnover.both[[2]] + labs(tag="D.")
write.table(sp.turnover.both[[1]],
            file=sprintf("tables/transient_microbe_poll/complete-turnover.txt"), sep="&")
write.csv(sp.turnover.both[[1]],
          file=sprintf("tables/transient_microbe_poll/complete-turnover.csv"))

pdf("figures/transient_microbe_poll/grid_transient_microbes.pdf", width = 8.5, height = 8.5) # Open a new pdf file
grid.arrange(panelA,
             panelB,
             panelC,
             panelD,
             ncol=2) # Write the grid.arrange in the file
dev.off()


