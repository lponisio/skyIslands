
prep_obligate_network <- function(raw_network,
                                  specimens_data){
  site_list <- names(raw_network)

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
  
  obligates_rows <- rownames(raw_network[[x]])
  
  ob_rows_to_keep <- grep(paste(these_obligates, collapse = "|"), obligates_rows)
  
  ob_new_net <- spNet_micro[[x]][ob_rows_to_keep,]
  
  new_name <- x
  
  only_obligate_network[[new_name]] <- ob_new_net
  
}

#### species level networks
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


geo <- unique(specimens_data[, c("Site", "Lat", "Long")])
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

}


run_network_turnover_mod <- function(this_component, this_network){
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
               iter = 20000,
               chains = 1,
               thin=1,
               init=0,
               control = list(adapt_delta = 0.99),
               save_pars = save_pars(all = TRUE))
  #browser()
  mod1
}


plot_network_turnover_mod_single <- function(mod1, 
                                      network_type,
                                      this.effect,
                                      this.resp,
                                      ){
  

  mod_summary <- write.summ.table(mod1)
  model_geodist <- mod_summary[rownames(mod_summary) == "GeoDist",]
  
  if(network_type == "Obligate") {
    point_color <- "darkgreen"
    if(model_geodist$Pgt0 >= 0.95){
      ribbon_color <- "darkgreen"
    } else if (model_geodist$Pgt0 <= 0.05) {
      ribbon_color <- "darkgreen"
    } else {ribbon_color <- NA}
  }
  
  if(network_type == "Transient") {
    point_color <- "darkorange"
    if(model_geodist$Pgt0 >= 0.95){
      ribbon_color <- "darkorange"
    } else if (model_geodist$Pgt0 <= 0.05) {
      ribbon_color <- "darkorange"
    } else {ribbon_color <- NA}
  }
  
  # Extract the data from conditional_effects
  cond_effects_data <- conditional_effects(mod1, effects = this.effect, resp = this.resp, plot = FALSE)
  plot_data <- cond_effects_data[[paste(this.resp,".",this.resp, "_",this.effect, sep='')]]
  
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill = ribbon_color,
                color = ribbon_color, linetype='dotted') +
    geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=ribbon_color, 
                color = ribbon_color, linetype='dashed') +
    geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=ribbon_color,
                color = ribbon_color, linetype='solid') +
    #Add line for the estimates
    geom_line(data = plot_data, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data, color = point_color, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add points for original data
    geom_point(data = point.data, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = point_color, alpha = 0.6,color="black", pch=21, cex=3) +
    # Labels and theme
    labs(x = xlabel, y = ylabel) +
    scale_x_continuous(breaks = axis.breaks, labels = axis.labs) +
    theme_classic()  +
    labs(x = "Geographic Distance (km)", y = label,
         fill = "Credible Interval") +
    theme_classic() +
    geom_point(data=this_network,
               aes(y=this_component, x=GeoDist), fill=point_color, color="black", pch=21, cex=2, alpha=0.9) + ylim(0,1) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  plot_obj
}