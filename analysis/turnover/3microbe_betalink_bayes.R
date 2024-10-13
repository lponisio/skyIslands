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

## prep obligate and transient network
only_obligate_network <- prep_obligate_network(raw_network=spNet_micro)

only_transient_network <- prep_transient_network(raw_network=spNet_micro)

## run network betalinkr prep script
source("src/betalinkrPrep.R")


## Run mods
run.mods=FALSE
if (run.mods==TRUE){
## Rewiring
rewiring.obligate.mod <- run_network_turnover_mod(this_component="OnlySharedLinks",
                                                  this_network=obligate_poll_betalink)

rewiring.transient.mod <- run_network_turnover_mod(this_component="OnlySharedLinks",
                                                   this_network=transient_poll_betalink)

## Host-driven turnover
host.driven.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsencePollinators",
                                                     this_network=obligate_poll_betalink)

host.driven.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsencePollinators",
                                                      this_network=transient_poll_betalink)

## Microbe-driven turnover
microbe.driven.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceMicrobes",
                                                        this_network=obligate_poll_betalink)

microbe.driven.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceMicrobes",
                                                         this_network=transient_poll_betalink)

## Complete turnover
complete.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceBoth",
                                                  this_network=obligate_poll_betalink)

complete.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceBoth",
                                                   this_network=transient_poll_betalink)

## save out models
save(rewiring.obligate.mod,
       rewiring.transient.mod,
       host.driven.obligate.mod,
       host.driven.transient.mod,
       microbe.driven.obligate.mod,
       microbe.driven.transient.mod,
       complete.obligate.mod,
       complete.transient.mod,
       file="../microbiome/saved/turnover_mods.Rdata")
} else {
  load("../microbiome/saved/turnover_mods.Rdata")
}


## A. rewiring

rewiring.plot <- plot_network_turnover_mod_compare(mod1=rewiring.obligate.mod,
                                              mod2=rewiring.transient.mod,
                                              this.network1=obligate_poll_betalink,
                                              this.network2=transient_poll_betalink,
                                              network_type1='Obligate',
                                              network_type2='Transient',
                                              this.effect="GeoDist",
                                              this.resp="OnlySharedLinks",
                                              label="Rewiring")
rewiring.plot[[1]]

panelA <- rewiring.plot[[1]] + labs(tag="A.")
rewiring.table <- rewiring.plot[[2]]

## B. Host-driven turnover

host.driven.plot <- plot_network_turnover_mod_compare(mod1=host.driven.obligate.mod,
                                                   mod2=host.driven.transient.mod,
                                                   this.network1=obligate_poll_betalink,
                                                   this.network2=transient_poll_betalink,
                                                   network_type1='Obligate',
                                                   network_type2='Transient',
                                                   this.effect="GeoDist",
                                                   this.resp="TurnoverAbsencePollinators",
                                                   label="Host-Driven Turnover")
host.driven.plot[[1]]
panelB <- host.driven.plot[[1]] + labs(tag="B.")
host.table <- host.driven.plot[[2]]

## C. Microbe-driven turnover

microbe.driven.plot <- plot_network_turnover_mod_compare(mod1=microbe.driven.obligate.mod,
                                                      mod2=microbe.driven.transient.mod,
                                                      this.network1=obligate_poll_betalink,
                                                      this.network2=transient_poll_betalink,
                                                      network_type1='Obligate',
                                                      network_type2='Transient',
                                                      this.effect="GeoDist",
                                                      this.resp="TurnoverAbsenceMicrobes",
                                                      label="Microbe-Driven Turnover")
microbe.driven.plot[[1]]
panelC <- microbe.driven.plot[[1]] + labs(tag="C.")
microbe.table <- microbe.driven.plot[[2]]

## D. Complete turnover

complete.plot <- plot_network_turnover_mod_compare(mod1=complete.obligate.mod,
                                                         mod2=complete.transient.mod,
                                                         this.network1=obligate_poll_betalink,
                                                         this.network2=transient_poll_betalink,
                                                         network_type1='Obligate',
                                                         network_type2='Transient',
                                                         this.effect="GeoDist",
                                                         this.resp="TurnoverAbsenceBoth",
                                                         label="Complete Turnover")

complete.plot[[1]]
panelD <- complete.plot[[1]] + labs(tag="D.")
complete.table <- complete.plot[[2]]

## Make panel figs and save out
pdf("../microbiome/figures/final/turnover_combined.pdf", width = 8.5, height = 8.5) # Open a new pdf file
grid.arrange(panelA,
             panelB,
             panelC,
             panelD,
             ncol=2) 
dev.off()

## Combine results tables and save out
combined.table <- bind_rows(rewiring.table,
                            host.table,
                            microbe.table,
                            complete.table)

write.csv(combined.table,
          file=sprintf("../microbiome/saved/tables/turnover.csv"))
