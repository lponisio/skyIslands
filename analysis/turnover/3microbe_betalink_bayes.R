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
library(performance)

load("../../data/networks/microNets.RData")

load("../../data/spec_RBCL_16s.RData")


whole.microbe.network = FALSE

obligate.microbe.network = FALSE

transient.microbe.network = TRUE
#reworking script to run brms models

##adapted from Lauren's 1betalink in skyIslands folder



#### species level networks
CH <- spNet_micro$CH
HM <- spNet_micro$HM
JC <- spNet_micro$JC
MM <- spNet_micro$MM 
PL <- spNet_micro$PL
RP <- spNet_micro$RP
SC <- spNet_micro$SC 
SM <- spNet_micro$SM

lower.order <- "Microbes"
higher.order <- "Pollinators"


microbe_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, RP, SM, SC),
                                         partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)

#View(microbe_poll_betalink)

colnames(microbe_poll_betalink) <- c("Site1",
                                     "Site2",
                                     "DissimilaritySpeciesComposition", #beta_S, the dissimilarity in species composition
                                     "OnlySharedLinks", #beta_OS, the dissimilarity (component) explained by “rewiring” among shared species
                                     "WholeNetworkLinks", #beta_WN, the dissimilarity between the two networks
                                     "SpeciesTurnoverLinks", #beta_ST, the dissimilarity (component) explained by difference in species community composition
                                     paste("TurnoverAbsence",lower.order,sep=""), #additive dissimilarity component (beta_ST) due to absence of resource species (lower, ST.l)
                                     paste("TurnoverAbsence",higher.order,sep=""), #dissimilarity due to absence of consumer species (higher, ST.h)
                                     "TurnoverAbsenceBoth") #dissimilarity due to absence of both






###will need to update LP's function networkBetaDiversity because most of the packages
### are no longer compatible :( 

geo <- unique(spec.net[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site

## add column for geographic distance between sites
microbe_poll_betalink$GeoDist <- apply(microbe_poll_betalink, 1, function(x){
  geo.dist[x["Site1"],  x["Site2"]]
})





##############################################################################

# Define a function called 'calculate_and_plot_betalinkr' that takes three arguments: this_component (which dissimilarity component should you make plots for?), 
#                                                                                     this_network(which network object will be used for input), and label (y axis label).


calculate_and_plot_betalinkr <- function(this_component, this_network, label){
  
  # Assign the value of 'this_component' to a new variable 'y'.
  y <- this_component
  
  #brms won't let me run the model unless there is a column called this_component, so I copied the variable of interest into a new
  # column called this_component.... janky? maybe
  this_network$this_component <- this_network[[this_component]]
  
  #browser()
  # Define a formula for brms with this_component the response variable,
  # 'GeoDist' as a fixed effect, and 'Site1' and 'Site2' as random effects.
  forms <- bf(formula(this_component~GeoDist + (1|Site1) + (1|Site2)))
  
  # Fit model
  mod1 <-  brm(forms, this_network,
              cores=1,
              iter = 10000,
              chains = 1,
              thin=1,
              init=0,
              control = list(adapt_delta = 0.99),
              save_pars = save_pars(all = TRUE))
  
  mod_summary <- write.summ.table(mod1)
  
  diag <- plot(check_model(mod1, panel = TRUE))
  
  newdata.mod <- tidyr::crossing(GeoDist = seq(min(this_network$GeoDist),
                                            max(this_network$GeoDist),
                                            length.out=10),
                                      Site1="CH",
                                      Site2="JC")
  
  #predictions
  pred_mod <- mod1 %>%
    epred_draws(newdata = newdata.mod ,
                resp = this_component,
                allow_new_levels = TRUE)
  
  
  fig <- ggplot(pred_mod, aes(x = GeoDist, y =.epred)) +
    stat_lineribbon(show.legend = FALSE) +
    scale_fill_brewer(palette = "Oranges") +
    labs(x = "Geographic Distance", y = label,
         fill = "Credible Interval") +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none") +
    theme_classic() +
    geom_point(data=this_network,
               aes(y=this_component, x=GeoDist), cex=2, alpha=0.5) + ylim(0,1)

  # Return a list containing the model summary[1] and the generated 'turnover.plot'[2].
  return(list(mod_summary, fig, diag, mod1))
  
  
}

################################################################################

if(whole.microbe.network){

dir.create("figures", showWarnings = FALSE)
dir.create("saved", showWarnings = FALSE)
dir.create("figures/diagnostic_plots", showWarnings = FALSE)
dir.create("figures/microbe_poll", showWarnings = FALSE)

species.turnover <- calculate_and_plot_betalinkr("DissimilaritySpeciesComposition",
                                                 microbe_poll_betalink,
                                                 "Dissimilarity: \nSpecies Composition")
ggsave(species.turnover[[2]], file="figures/microbe_poll/DissimilaritySpeciesTurnover.pdf", height=4, width=6)
ggsave(species.turnover[[3]], file="figures/diagnostic_plots/DissimilaritySpeciesTurnover.pdf", height=8, width=11)
write.csv(species.turnover[[1]], file="tables/DissimilaritySpeciesTurnover.csv")
save(species.turnover, file="saved/DissimilaritySpeciesTurnover.Rdata")

interaction.turnover <- calculate_and_plot_betalinkr("WholeNetworkLinks",
                                                     microbe_poll_betalink,
                                                     "Dissimilarity: \nInteraction Turnover")
ggsave(interaction.turnover[[2]], file="figures/microbe_poll/DissimilarityInteractionTurnover.pdf", height=4, width=6)
ggsave(interaction.turnover[[3]], file="figures/diagnostic_plots/DissimilarityInteractionTurnover.pdf", height=8, width=11)
write.csv(interaction.turnover[[1]], file="tables/DissimilarityInteractionTurnover.csv")
save(interaction.turnover, file="saved/DissimilarityInteractionTurnover.Rdata")


int.turnover.rewiring <- calculate_and_plot_betalinkr("OnlySharedLinks",
                                                      microbe_poll_betalink,
                                                      "Interaction Turnover: \nRewiring")
ggsave(int.turnover.rewiring[[2]], file="figures/microbe_poll/InteractionDissimilarityRewiring.pdf", height=4, width=6)
ggsave(int.turnover.rewiring[[3]], file="figures/diagnostic_plots/InteractionDissimilarityRewiring.pdf", height=8, width=11)
write.csv(int.turnover.rewiring[[1]], file="tables/InteractionDissimilarityRewiring.csv")
save(int.turnover.rewiring, file="saved/InteractionDissimilarityRewiring.Rdata")


int.turnover.species.turnover <- calculate_and_plot_betalinkr("SpeciesTurnoverLinks",
                                                              microbe_poll_betalink,
                                                              "Interaction Turnover: \nSpecies Turnover")
ggsave(int.turnover.species.turnover[[2]], file="figures/microbe_poll/InteractionTurnoverSpeciesComp.pdf", height=4, width=6)
ggsave(int.turnover.species.turnover[[3]], file="figures/diagnostic_plots/InteractionTurnoverSpeciesComp.pdf", height=8, width=11)
write.csv(int.turnover.species.turnover[[1]], file="tables/InteractionTurnoverSpeciesComp.csv")
save(int.turnover.species.turnover, file="saved/InteractionTurnoverSpeciesComp.Rdata")


sp.turnover.microbes <- calculate_and_plot_betalinkr("TurnoverAbsenceMicrobes",
                                                     microbe_poll_betalink,
                                                     "Species Turnover: \nAbsence of Microbes")
ggsave(sp.turnover.microbes[[2]], file="figures/microbe_poll/SpeciesTurnoverAbsenceMicrobes.pdf", height=4, width=6)
ggsave(sp.turnover.microbes[[3]], file="figures/diagnostic_plots/SpeciesTurnoverAbsenceMicrobes.pdf", height=8, width=11)
write.csv(sp.turnover.microbes[[1]], file="tables/SpeciesTurnoverAbsenceMicrobes.csv")
save(sp.turnover.microbes, file="saved/SpeciesTurnoverAbsenceMicrobes.Rdata")


sp.turnover.bees <- calculate_and_plot_betalinkr("TurnoverAbsencePollinators",
                                                 microbe_poll_betalink,
                                                 "Species Turnover: \nAbsence of Bees")
ggsave(sp.turnover.bees[[2]], file="figures/microbe_poll/SpeciesTurnoverAbsenceBees.pdf", height=4, width=6)
ggsave(sp.turnover.bees[[3]], file="figures/diagnostic_plots/SpeciesTurnoverAbsenceBees.pdf", height=8, width=11)
write.csv(sp.turnover.bees[[1]], file="tables/SpeciesTurnoverAbsenceBees.csv")
save(sp.turnover.bees, file="saved/SpeciesTurnoverAbsenceBees.Rdata")


sp.turnover.both <- calculate_and_plot_betalinkr("TurnoverAbsenceBoth",
                                                 microbe_poll_betalink,
                                                 "Species Turnover: \nAbsence of Both")
ggsave(sp.turnover.both[[2]], file="figures/microbe_poll/SpeciesTurnoverAbsenceBoth.pdf", height=4, width=6)
ggsave(sp.turnover.both[[3]], file="figures/diagnostic_plots/SpeciesTurnoverAbsenceBoth.pdf", height=8, width=11)
write.csv(sp.turnover.both[[1]], file="tables/SpeciesTurnoverAbsenceBoth.csv")
save(sp.turnover.both, file="saved/SpeciesTurnoverAbsenceBoth.Rdata")

# diss.whole.network.replace <- calculate_and_plot_betalinkr(microbe_poll_betalink$WholeNetworkReplaced, microbe_poll_betalink, "Dissimilarity: Whole Network Replacement")
# ggsave(diss.whole.network.replace[[2]], file="figures/microbe_poll/DissimilarityWholeNetReplace.pdf", height=4, width=6)
#
# diss.only.shared.replace <- calculate_and_plot_betalinkr(microbe_poll_betalink$OnlySharedReplaced, microbe_poll_betalink, "Dissimilarity: Only Shared Species Replacement")
# ggsave(diss.only.shared.replace[[2]], file="figures/microbe_poll/DissimilarityOnlySharedReplace.pdf", height=4, width=6)
#
# diss.whole.network.rich.dif <- calculate_and_plot_betalinkr(microbe_poll_betalink$WholeNetworkRichnessDifference, microbe_poll_betalink, "Dissimilarity: Whole Network Richness Difference")
# ggsave(diss.whole.network.rich.dif[[2]], file="figures/microbe_poll/DissimilarityWholeNetRichDif.pdf", height=4, width=6)
#
# diss.only.shared.rich.dif <- calculate_and_plot_betalinkr(microbe_poll_betalink$OnlySharedRichnessDifference, microbe_poll_betalink, "Dissimilarity: Only Shared Richness Difference")
# ggsave(diss.only.shared.rich.dif[[2]], file="figures/microbe_poll/DissimilarityOnlySharedRichDif.pdf", height=4, width=6)

## make panel figure
panelA <- species.turnover[[2]] + labs(tag="A.")
panelB <- interaction.turnover[[2]] + labs(tag="B.")
panelC <- sp.turnover.bees[[2]] + labs(tag="C.")
panelD <- int.turnover.species.turnover[[2]] + labs(tag="D.")
panelE <- sp.turnover.microbes[[2]] + labs(tag="E.")
panelF <- int.turnover.rewiring[[2]] + labs(tag="F.")


grid.arrange(panelA,
             panelB,
             panelC,
             panelD,
             panelE,
             panelF,
             ncol=2)}
#################################################################

if(obligate.microbe.network==TRUE){

## now do the same for the three groups of obligate, parasitic, transient

site_list <- names(spNet_micro)

## obligate symbionts
these_obligates <- c("Lactobacillus",
                     "Bifidobacterium",
                     "Snodgrassella",
                     "Gilliamella",
                     "Frischella",
                     "Bartonella",
                     "Commensalibacter")


only_obligate_network <- list()

for (x in site_list){
  
  obligates_rows <- rownames(spNet_micro[[x]])
  
  ob_rows_to_keep <- grep(paste(these_obligates, collapse = "|"), obligates_rows)
  
  ob_new_net <- spNet_micro[[x]][ob_rows_to_keep,]
  
  new_name <- x
  
  only_obligate_network[[new_name]] <- ob_new_net
  
}

#### species level networks
CH <- only_obligate_network$CH
HM <- only_obligate_network$HM
JC <- only_obligate_network$JC
MM <- only_obligate_network$MM 
PL <- only_obligate_network$PL
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

species.turnover <- calculate_and_plot_betalinkr("DissimilaritySpeciesComposition",
                                                 obligate_poll_betalink,
                                                 "Dissimilarity: \nSpecies Composition")
ggsave(species.turnover[[2]], file="figures/obligate_microbe_poll/DissimilaritySpeciesTurnover.pdf", height=6, width=6)

interaction.turnover <- calculate_and_plot_betalinkr("WholeNetworkLinks",
                                                     obligate_poll_betalink,
                                                     "Dissimilarity: \nInteraction Turnover")
ggsave(interaction.turnover[[2]], file="figures/obligate_microbe_poll/DissimilarityInteractionTurnover.pdf", height=6, width=6)

int.turnover.rewiring <- calculate_and_plot_betalinkr("OnlySharedLinks",
                                                      obligate_poll_betalink,
                                                      "Interaction Turnover: \nRewiring")
ggsave(int.turnover.rewiring[[2]], file="figures/obligate_microbe_poll/InteractionDissimilarityRewiring.pdf", height=6, width=6)

int.turnover.species.turnover <- calculate_and_plot_betalinkr("SpeciesTurnoverLinks",
                                                              obligate_poll_betalink,
                                                              "Interaction Turnover: \nSpecies Turnover")
ggsave(int.turnover.species.turnover[[2]], file="figures/obligate_microbe_poll/InteractionTurnoverSpeciesComp.pdf", height=6, width=6)

sp.turnover.microbes <- calculate_and_plot_betalinkr("TurnoverAbsenceMicrobes",
                                                     obligate_poll_betalink,
                                                     "Species Turnover: \nAbsence of Microbes")
ggsave(sp.turnover.microbes[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceMicrobes.pdf", height=6, width=6)

sp.turnover.bees <- calculate_and_plot_betalinkr("TurnoverAbsencePollinators",
                                                 obligate_poll_betalink,
                                                 "Species Turnover: \nAbsence of Bees")
ggsave(sp.turnover.bees[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceBees.pdf", height=6, width=6)

sp.turnover.both <- calculate_and_plot_betalinkr("TurnoverAbsenceBoth",
                                                 obligate_poll_betalink,
                                                 "Species Turnover: \nAbsence of Both")
ggsave(sp.turnover.both[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceBoth.pdf", height=6, width=6)

## make panel figure
panelA <- species.turnover[[2]] + labs(tag="A.")
panelB <- interaction.turnover[[2]] + labs(tag="B.")
panelC <- sp.turnover.bees[[2]] + labs(tag="C.")
panelD <- int.turnover.species.turnover[[2]] + labs(tag="D.")
panelE <- sp.turnover.microbes[[2]] + labs(tag="E.")
panelF <- int.turnover.rewiring[[2]] + labs(tag="F.")


grid.arrange(panelA,
             panelB,
             panelC,
             panelD,
             panelE,
             panelF,
             ncol=2)
}

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
HM <- only_transient_network$HM
JC <- only_transient_network$JC
MM <- only_transient_network$MM 
PL <- only_transient_network$PL
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
                                                 "Dissimilarity: \nSpecies Composition")
ggsave(species.turnover[[2]], file="figures/transient_microbe_poll/DissimilaritySpeciesTurnover.pdf", height=6, width=6)

interaction.turnover <- calculate_and_plot_betalinkr("WholeNetworkLinks",
                                                     transient_poll_betalink,
                                                     "Dissimilarity: \nInteraction Turnover")
ggsave(interaction.turnover[[2]], file="figures/transient_microbe_poll/DissimilarityInteractionTurnover.pdf", height=6, width=6)

int.turnover.rewiring <- calculate_and_plot_betalinkr("OnlySharedLinks",
                                                      transient_poll_betalink,
                                                      "Interaction Turnover: \nRewiring")
ggsave(int.turnover.rewiring[[2]], file="figures/transient_microbe_poll/InteractionDissimilarityRewiring.pdf", height=6, width=6)

int.turnover.species.turnover <- calculate_and_plot_betalinkr("SpeciesTurnoverLinks",
                                                              transient_poll_betalink,
                                                              "Interaction Turnover: \nSpecies Turnover")
ggsave(int.turnover.species.turnover[[2]], file="figures/transient_microbe_poll/InteractionTurnoverSpeciesComp.pdf", height=6, width=6)

sp.turnover.microbes <- calculate_and_plot_betalinkr("TurnoverAbsenceMicrobes",
                                                     transient_poll_betalink,
                                                     "Species Turnover: \nAbsence of Microbes")
ggsave(sp.turnover.microbes[[2]], file="figures/transient_microbe_poll/SpeciesTurnoverAbsenceMicrobes.pdf", height=6, width=6)

sp.turnover.bees <- calculate_and_plot_betalinkr("TurnoverAbsencePollinators",
                                                 transient_poll_betalink,
                                                 "Species Turnover: \nAbsence of Bees")
ggsave(sp.turnover.bees[[2]], file="figures/transient_microbe_poll/SpeciesTurnoverAbsenceBees.pdf", height=6, width=6)

sp.turnover.both <- calculate_and_plot_betalinkr("TurnoverAbsenceBoth",
                                                 transient_poll_betalink,
                                                 "Species Turnover: \nAbsence of Both")
ggsave(sp.turnover.both[[2]], file="figures/transient_microbe_poll/SpeciesTurnoverAbsenceBoth.pdf", height=6, width=6)

## make panel figure
panelA <- species.turnover[[2]] + labs(tag="A.")
panelB <- interaction.turnover[[2]] + labs(tag="B.")
panelC <- sp.turnover.bees[[2]] + labs(tag="C.")
panelD <- int.turnover.species.turnover[[2]] + labs(tag="D.")
panelE <- sp.turnover.microbes[[2]] + labs(tag="E.")
panelF <- int.turnover.rewiring[[2]] + labs(tag="F.")


grid.arrange(panelA,
             panelB,
             panelC,
             panelD,
             panelE,
             panelF,
             ncol=2)

}