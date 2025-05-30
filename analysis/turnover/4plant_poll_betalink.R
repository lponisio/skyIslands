## setwd("~/Dropbox/skyIslands/")
setwd('~/Dropbox (University of Oregon)/skyIslands/') ## Rebecca wd

rm(list=ls())
setwd("analysis/turnover")
source("src/initialize.R")
source("src/chao.R")
source("src/betaNet.R")
library(ggplot2)
library(lme4)
library(lmerTest)
library(igraph)
library(ggpubr)
library(emmeans)

load("C:/Users/rah10/Dropbox (University of Oregon)/skyIslands/data/networks/rbclNets.RData")


##adapted from Lauren's 1betalink in skyIslands folder



####################### must use functions from bipartite bc betalink package is depreciated

CH <- spNet_rbcl$CH
HM <- spNet_rbcl$HM
JC <- spNet_rbcl$JC
MM <- spNet_rbcl$MM 
PL <- spNet_rbcl$PL
SC <- spNet_rbcl$SC 
SM <- spNet_rbcl$SM

lower.order <- "Plants"
higher.order <- "Pollinators"


plant_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC),
                                         partitioning="commondenom", binary=TRUE, distofempty='zero', partition.st=TRUE, partition.rr=TRUE)

#View(plant_poll_betalink)

colnames(plant_poll_betalink) <- c("Site1",
                                     "Site2",
                                     "DissimilaritySpeciesComposition",
                                     "OnlySharedLinks",
                                     "WholeNetworkLinks",
                                     "SpeciesTurnoverLinks",
                                     paste("TurnoverAbsence",lower.order,sep=""),
                                     paste("TurnoverAbsence",higher.order,sep=""),
                                     "TurnoverAbsenceBoth",
                                     "WholeNetworkReplaced",
                                     "OnlySharedReplaced",
                                     "WholeNetworkRichnessDifference",
                                     "OnlySharedRichnessDifference"
)






###will need to update LP's function networkBetaDiversity because most of the packages
### are no longer compatible :( 

geo <- unique(spec.net[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site

## add column for geographic distance between sites
plant_poll_betalink$GeoDist <- apply(plant_poll_betalink, 1, function(x){
  geo.dist[x["Site1"],  x["Site2"]]
})





###############################################################################

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
    lower.bound = -0.10
  } else {
    lower.bound = 0
  }
  # Create a ggplot for plotting the data in 'this_network' with 'GeoDist' on the x-axis and 'this_component' on the y-axis.
  
  turnover.plot <- ggplot(this_network, 
                          aes(x=GeoDist, y=this_component)) +
    geom_point() + 
    # Add a ribbon to the plot based on 'emm1' for visualizing confidence intervals.
    geom_ribbon(data= emm1, aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= ribbon_col, alpha=0.5) +
    # Add a line to the plot based on 'emm1' for the estimated marginal means.
    geom_line(data= emm1, aes(y= emmean)) +
    theme_classic() + 
    labs(x='Geographic Distance (km)', y=label) +
    # Set the y-axis limits based on the computed lower and upper bounds.
    scale_y_continuous(limits=c(lower.bound,upper.bound))
  # Return a list containing the model summary[1] and the generated 'turnover.plot'[2].
  return(list(mod_summary, turnover.plot))
  
  
}

################################################################################


dir.create("figures/plant_poll", showWarnings = FALSE)

species.turnover <- calculate_and_plot_betalinkr(plant_poll_betalink$DissimilaritySpeciesComposition, plant_poll_betalink, "Species Turnover")
ggsave(species.turnover[[2]], file="figures/plant_poll/DissimilaritySpeciesTurnover.pdf", height=4, width=6)

interaction.turnover <- calculate_and_plot_betalinkr(plant_poll_betalink$WholeNetworkLinks, plant_poll_betalink, "Interaction Turnover")
ggsave(interaction.turnover[[2]], file="figures/plant_poll/DissimilarityInteractionTurnover.pdf", height=4, width=6)

int.turnover.rewiring <- calculate_and_plot_betalinkr(plant_poll_betalink$OnlySharedLinks, plant_poll_betalink, "Interaction Turnover: Rewiring")
ggsave(int.turnover.rewiring[[2]], file="figures/plant_poll/InteractionDissimilarityRewiring.pdf", height=4, width=6)

int.turnover.species.turnover <- calculate_and_plot_betalinkr(plant_poll_betalink$SpeciesTurnoverLinks, plant_poll_betalink, "Interaction Turnover: Species Turnover")
ggsave(int.turnover.species.turnover[[2]], file="figures/plant_poll/InteractionTurnoverSpeciesComp.pdf", height=4, width=6)

sp.turnover.plants <- calculate_and_plot_betalinkr(plant_poll_betalink$TurnoverAbsencePlants, plant_poll_betalink, "Species Turnover: plants")
ggsave(sp.turnover.plants[[2]], file="figures/plant_poll/SpeciesTurnoverAbsencePlants.pdf", height=4, width=6)

sp.turnover.bees <- calculate_and_plot_betalinkr(plant_poll_betalink$TurnoverAbsencePollinators, plant_poll_betalink, "Species Turnover: Bees")
ggsave(sp.turnover.bees[[2]], file="figures/plant_poll/SpeciesTurnoverAbsenceBees.pdf", height=4, width=6)

sp.turnover.both <- calculate_and_plot_betalinkr(plant_poll_betalink$TurnoverAbsenceBoth, plant_poll_betalink, "Species Turnover: Both")
ggsave(sp.turnover.both[[2]], file="figures/plant_poll/SpeciesTurnoverAbsenceBoth.pdf", height=4, width=6)

diss.whole.network.replace <- calculate_and_plot_betalinkr(plant_poll_betalink$WholeNetworkReplaced, plant_poll_betalink, "Dissimilarity: Whole Network Replacement")
ggsave(diss.whole.network.replace[[2]], file="figures/plant_poll/DissimilarityWholeNetReplace.pdf", height=4, width=6)

diss.only.shared.replace <- calculate_and_plot_betalinkr(plant_poll_betalink$OnlySharedReplaced, plant_poll_betalink, "Dissimilarity: Only Shared Species Replacement")
ggsave(diss.only.shared.replace[[2]], file="figures/plant_poll/DissimilarityOnlySharedReplace.pdf", height=4, width=6)

diss.whole.network.rich.dif <- calculate_and_plot_betalinkr(plant_poll_betalink$WholeNetworkRichnessDifference, plant_poll_betalink, "Dissimilarity: Whole Network Richness Difference")
ggsave(diss.whole.network.rich.dif[[2]], file="figures/plant_poll/DissimilarityWholeNetRichDif.pdf", height=4, width=6)

diss.only.shared.rich.dif <- calculate_and_plot_betalinkr(plant_poll_betalink$OnlySharedRichnessDifference, plant_poll_betalink, "Dissimilarity: Only Shared Richness Difference")
ggsave(diss.only.shared.rich.dif[[2]], file="figures/plant_poll/DissimilarityOnlySharedRichDif.pdf", height=4, width=6)



# ### facet plots
# bee_microbe_turnover <- ggarrange(spec.turnover.plot,
#                                   pollinator.turnover.plot,
#                                   microbe.turnover.plot,
#                                   interaction.turnover.plot,
#                                   int.turnover.spcomp.plot,
#                                   int.turnover.rewiring.plot,
#                                   ncol=3, nrow=2)
# 
# bee_microbe_turnover
# 
# dir.create("figures")
# 
# ggsave(bee_microbe_turnover, file="figures/poll_microbe_betaComponents.pdf",
#        height=8, width=11)


########## combined plots

# Assign the value of 'this_component' to a new variable 'y'.
y <- plant_poll_betalink$OnlySharedLinks
y2 <- plant_poll_betalink$SpeciesTurnoverLinks

# Define a formula for linear mixed-effects modeling with 'y' as the response variable,
# 'GeoDist' as a fixed effect, and 'Site1' and 'Site2' as random effects.
forms <- formula(y~GeoDist + (1|Site1) + (1|Site2))
forms2 <- formula(y2~GeoDist + (1|Site1) + (1|Site2))

# Fit a linear mixed-effects model using the 'lmer' function, providing the formula, data, and specifying REML to be FALSE.
mod2 <- do.call(lmer,
                list(formula=forms2,
                     data=plant_poll_betalink,
                     REML = FALSE))
mod1 <- do.call(lmer,
                list(formula=forms,
                     data=plant_poll_betalink,
                     REML = FALSE))

# Summarize the model 'mod1'.
mod_summary1 <- summary(mod1)
mod_summary2 <- summary(mod2)

# Compute the significance of the model using ANOVA.
model.sig1 <- anova(mod1)
model.sig2 <- anova(mod2)

# Check if the p-value from the ANOVA is less than or equal to 0.05 and set the 'ribbon_col' accordingly.
if(model.sig2$`Pr(>F)` <= 0.05){
  ribbon_col2 = '#FFA500'
} else {
  ribbon_col2 = 'grey80'
}

if(model.sig1$`Pr(>F)` <= 0.05){
  ribbon_col1 = '#FFA500'
} else {
  ribbon_col1 = 'grey80'
}

# Create a reference grid 'gr1' based on the model, keeping 'GeoDist' as a covariate.
gr2 <- ref_grid(mod2, cov.keep= c('GeoDist'))
gr1 <- ref_grid(mod1, cov.keep= c('GeoDist'))

# Compute estimated marginal means (emmeans) for 'GeoDist' at a 95% confidence level and store it in 'emm1'.
emm2 <- data.frame(emmeans(gr2, spec= c("GeoDist"), level= 0.95))
emm1 <- data.frame(emmeans(gr1, spec= c("GeoDist"), level= 0.95))

# Determine the upper bound based on 'emm1' and set it to 1.05 if the upper confidence limit is greater than 1, or 1 otherwise.
if(max(emm1$upper.CL) > 1){
  upper.bound1 = 1.05
} else {
  upper.bound1 = 1
}

if(max(emm2$upper.CL) > 1){
  upper.bound2 = 1.05
} else {
  upper.bound2 = 1
}


# Determine the lower bound based on 'emm1' and set it to -0.10 if the lower confidence limit is less than 0, or 0 otherwise.
if(min(emm1$lower.CL) < 0){
  lower.bound1 = -0.10
} else {
  lower.bound1 = 0
}

if(min(emm2$lower.CL) < 0){
  lower.bound2 = -0.10
} else {
  lower.bound2 = 0
}
# Create a ggplot for plotting the data in 'this_network' with 'GeoDist' on the x-axis and 'this_component' on the y-axis.

turnover.plot <- ggplot(plant_poll_betalink, 
                        aes(x=GeoDist, y=OnlySharedLinks)) +
  geom_point() + 
  # Add a ribbon to the plot based on 'emm1' for visualizing confidence intervals.
  geom_ribbon(data= emm1, aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= ribbon_col1, alpha=0.5) +
  # Add a line to the plot based on 'emm1' for the estimated marginal means.
  geom_line(data= emm1, aes(y= emmean)) +
  theme_classic() + 
  geom_point(aes(x=GeoDist, y=SpeciesTurnoverLinks)) +
  # Add a ribbon to the plot based on 'emm1' for visualizing confidence intervals.
  geom_ribbon(data= emm2, aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= ribbon_col2, alpha=0.5) +
  geom_line(data= emm2, aes(y= emmean)) +
  labs(x='Geographic Distance (km)', y="Network Dissimilarity") +
  # Set the y-axis limits based on the computed lower and upper bounds.
  scale_y_continuous(limits=c(lower.bound2,upper.bound))

turnover.plot
