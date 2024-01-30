## setwd("~/Dropbox/skyIslands/")
setwd('~/Dropbox (University of Oregon)/skyIslands/') ## Rebecca wd

rm(list=ls())
setwd("analysis/turnover")
#source("src/initialize.R")
source("src/chao.R")
source("src/betaNet.R")
library(ggplot2)
library(lme4)
library(lmerTest)
library(igraph)
library(ggpubr)
library(emmeans)
library(bipartite)
library(dplyr)
library(fields)

library(gridExtra)

load("C:/Users/rah10/Dropbox (University of Oregon)/skyIslands/data/networks/microNets.RData")

load("C:/Users/rah10/Dropbox (University of Oregon)/skyIslands/data/spec_RBCL_16s.RData")



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
  
  # Define a formula for linear mixed-effects modeling with 'y' as the response variable,
  # 'GeoDist' as a fixed effect, and 'Site1' and 'Site2' as random effects.
  forms <- formula(y~GeoDist + (1|Site1) + (1|Site2))
  
  # Fit a linear mixed-effects model using the 'lmer' function, providing the formula, data, and specifying REML to be FALSE.
  mod1 <- do.call(lmer,
                  list(formula=forms,
                       data=this_network,
                       REML = FALSE))
  
  # Summarize the model 'mod1'.
  mod_summary <- summary(mod1)
  
  # Compute the significance of the model using ANOVA.
  model.sig <- anova(mod1)
  
  # Check if the p-value from the ANOVA is less than or equal to 0.05 and set the 'ribbon_col' accordingly.
  if(model.sig$`Pr(>F)` <= 0.05){
    ribbon_col = '#FFA500'
  } else {
    ribbon_col = 'grey80'
  }
  
  # Create a reference grid 'gr1' based on the model, keeping 'GeoDist' as a covariate.
  gr1 <- ref_grid(mod1, cov.keep= c('GeoDist'))
  
  # Compute estimated marginal means (emmeans) for 'GeoDist' at a 95% confidence level and store it in 'emm1'.
  emm1 <- data.frame(emmeans(gr1, spec= c("GeoDist"), level= 0.95))
  
  # Determine the upper bound based on 'emm1' and set it to 1.05 if the upper confidence limit is greater than 1, or 1 otherwise.
  if(max(emm1$upper.CL) > 1){
    upper.bound = 1.05
  } else {
    upper.bound = 1
  }
  
  # Determine the lower bound based on 'emm1' and set it to -0.10 if the lower confidence limit is less than 0, or 0 otherwise.
  if(min(emm1$lower.CL) < 0){
    lower.bound = -0.15
  } else {
    lower.bound = 0
  }
  # Create a ggplot for plotting the data in 'this_network' with 'GeoDist' on the x-axis and 'this_component' on the y-axis.
  
  turnover.plot <- ggplot(this_network, 
                          aes(x=GeoDist, y=this_component)) +
    # Add a ribbon to the plot based on 'emm1' for visualizing confidence intervals.
    geom_ribbon(data= emm1, aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= ribbon_col, alpha=0.5) +
    # Add a line to the plot based on 'emm1' for the estimated marginal means.
    geom_line(data= emm1, aes(y= emmean)) +
    geom_point() + 
    theme_classic() + 
    labs(x='Geographic Distance (km)', y=label) +
    # Set the y-axis limits based on the computed lower and upper bounds.
    scale_y_continuous(limits=c(lower.bound,upper.bound))
  # Return a list containing the model summary[1] and the generated 'turnover.plot'[2].
  return(list(mod_summary, turnover.plot))
  
  
}

################################################################################



dir.create("figures", showWarnings = FALSE)
dir.create("figures/microbe_poll", showWarnings = FALSE)

species.turnover <- calculate_and_plot_betalinkr(microbe_poll_betalink$DissimilaritySpeciesComposition,
                                                 microbe_poll_betalink,
                                                 "Dissimilarity: \nSpecies Composition")
ggsave(species.turnover[[2]], file="figures/microbe_poll/DissimilaritySpeciesTurnover.pdf", height=4, width=6)

interaction.turnover <- calculate_and_plot_betalinkr(microbe_poll_betalink$WholeNetworkLinks,
                                                     microbe_poll_betalink,
                                                     "Dissimilarity: \nInteraction Turnover")
ggsave(interaction.turnover[[2]], file="figures/microbe_poll/DissimilarityInteractionTurnover.pdf", height=4, width=6)

int.turnover.rewiring <- calculate_and_plot_betalinkr(microbe_poll_betalink$OnlySharedLinks,
                                                      microbe_poll_betalink,
                                                      "Interaction Turnover: \nRewiring")
ggsave(int.turnover.rewiring[[2]], file="figures/microbe_poll/InteractionDissimilarityRewiring.pdf", height=4, width=6)

int.turnover.species.turnover <- calculate_and_plot_betalinkr(microbe_poll_betalink$SpeciesTurnoverLinks,
                                                              microbe_poll_betalink,
                                                              "Interaction Turnover: \nSpecies Turnover")
ggsave(int.turnover.species.turnover[[2]], file="figures/microbe_poll/InteractionTurnoverSpeciesComp.pdf", height=4, width=6)

sp.turnover.microbes <- calculate_and_plot_betalinkr(microbe_poll_betalink$TurnoverAbsenceMicrobes,
                                                     microbe_poll_betalink,
                                                     "Species Turnover: \nAbsence of Microbes")
ggsave(sp.turnover.microbes[[2]], file="figures/microbe_poll/SpeciesTurnoverAbsenceMicrobes.pdf", height=4, width=6)

sp.turnover.bees <- calculate_and_plot_betalinkr(microbe_poll_betalink$TurnoverAbsencePollinators,
                                                 microbe_poll_betalink,
                                                 "Species Turnover: \nAbsence of Bees")
ggsave(sp.turnover.bees[[2]], file="figures/microbe_poll/SpeciesTurnoverAbsenceBees.pdf", height=4, width=6)

sp.turnover.both <- calculate_and_plot_betalinkr(microbe_poll_betalink$TurnoverAbsenceBoth,
                                                 microbe_poll_betalink,
                                                 "Species Turnover: \nAbsence of Both")
ggsave(sp.turnover.both[[2]], file="figures/microbe_poll/SpeciesTurnoverAbsenceBoth.pdf", height=4, width=6)

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
             ncol=2)
#################################################################



## now do the same for the three groups of obligate, parasitic, transient

site_list <- names(spNet_micro)

## obligate symbionts
these_obligates <- c("Rosenbergiella", "Pseudomonas", "Gilliamella",
                                "Lactobacillus", "Caulobacter", "Snodgrassella",
                                "Acinetobacter", "Corynebacterium", "Sphingomonas",
                                "Commensalibacter", "Methylobacterium",
                                "Massilia","Stenotrophomonas", "Bifidobacterium", "Frischella", "Bartonella")

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

species.turnover <- calculate_and_plot_betalinkr(obligate_poll_betalink$DissimilaritySpeciesComposition,
                                                 obligate_poll_betalink,
                                                 "Dissimilarity: \nSpecies Composition")
ggsave(species.turnover[[2]], file="figures/obligate_microbe_poll/DissimilaritySpeciesTurnover.pdf", height=4, width=6)

interaction.turnover <- calculate_and_plot_betalinkr(obligate_poll_betalink$WholeNetworkLinks,
                                                     obligate_poll_betalink,
                                                     "Dissimilarity: \nInteraction Turnover")
ggsave(interaction.turnover[[2]], file="figures/obligate_microbe_poll/DissimilarityInteractionTurnover.pdf", height=4, width=6)

int.turnover.rewiring <- calculate_and_plot_betalinkr(obligate_poll_betalink$OnlySharedLinks,
                                                      obligate_poll_betalink,
                                                      "Interaction Turnover: \nRewiring")
ggsave(int.turnover.rewiring[[2]], file="figures/obligate_microbe_poll/InteractionDissimilarityRewiring.pdf", height=4, width=6)

int.turnover.species.turnover <- calculate_and_plot_betalinkr(obligate_poll_betalink$SpeciesTurnoverLinks,
                                                              obligate_poll_betalink,
                                                              "Interaction Turnover: \nSpecies Turnover")
ggsave(int.turnover.species.turnover[[2]], file="figures/obligate_microbe_poll/InteractionTurnoverSpeciesComp.pdf", height=4, width=6)

sp.turnover.microbes <- calculate_and_plot_betalinkr(obligate_poll_betalink$TurnoverAbsenceMicrobes,
                                                     obligate_poll_betalink,
                                                     "Species Turnover: \nAbsence of Microbes")
ggsave(sp.turnover.microbes[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceMicrobes.pdf", height=4, width=6)

sp.turnover.bees <- calculate_and_plot_betalinkr(obligate_poll_betalink$TurnoverAbsencePollinators,
                                                 obligate_poll_betalink,
                                                 "Species Turnover: \nAbsence of Bees")
ggsave(sp.turnover.bees[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceBees.pdf", height=4, width=6)

sp.turnover.both <- calculate_and_plot_betalinkr(obligate_poll_betalink$TurnoverAbsenceBoth,
                                                 obligate_poll_betalink,
                                                 "Species Turnover: \nAbsence of Both")
ggsave(sp.turnover.both[[2]], file="figures/obligate_microbe_poll/SpeciesTurnoverAbsenceBoth.pdf", height=4, width=6)

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

