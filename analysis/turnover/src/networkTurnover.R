

## WIP
# trying to break down function to get around issue of not finding CH
prep_obligate_network <- function(raw_network=spNet_micro){
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
    #browser()
    obligates_rows <- rownames(raw_network[[x]])
    
    ob_rows_to_keep <- grep(paste(these_obligates, collapse = "|"), obligates_rows)
    
    ob_new_net <- raw_network[[x]][ob_rows_to_keep,]
    print(dim(ob_new_net))
    new_name <- x
    
    only_obligate_network[[new_name]] <- ob_new_net
    
  }
  only_obligate_network
  #### species level networks
  
}





prep_transient_network <- function(raw_network=spNet_micro){
  site_list <- names(raw_network)
  
  ## obligate symbionts
  bee.obligates <- "Lactobacillus|Bifidobacterium|Snodgrassella|Gilliamella|Frischella|Bartonella|Commensalibacter"
  
  
  only_transient_network <- list()
  #browser()
  for (x in site_list){
    
    trans_rows <- rownames(raw_network[[x]])
    
    trans_rows_to_drop <- !grepl(bee.obligates, trans_rows)
    
    trans_new_net <- raw_network[[x]][trans_rows_to_drop,]
    print(dim(trans_new_net))
    new_name <- x
    
    only_transient_network[[new_name]] <- trans_new_net
    #browser()
  }
  only_transient_network
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

obligate_betalinkr <- function(only_obligate_network=only_obligate_network){
  
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
  
  return(obligate_poll_betalink)
}

transient_betalinkr <- function(only_transient_network=only_transient_network){
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
  
  return(transient_poll_betalink)
}



plot_network_turnover_mod_single <- function(mod1, 
                                      this.network,
                                      network_type,
                                      this.effect,
                                      this.resp,
                                      label
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
  plot_data <- cond_effects_data[[this.effect]]
  #browser()
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill = ribbon_color,
                color = point_color, linetype='dotted') +
    geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=ribbon_color, 
                color = point_color, linetype='dashed') +
    geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=ribbon_color,
                color = point_color, linetype='solid') +
    #Add line for the estimates
    geom_line(data = plot_data, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data, color = point_color, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add points for original data
    geom_point(data = this.network, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = point_color, alpha = 0.9,color="black", pch=21, cex=3) +
    # Labels and theme
    theme_classic()  +
    labs(x = "Geographic Distance (km)", y = label,
         fill = "Credible Interval") +
    theme_classic() +
    ylim(0,1) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  return(list(plot_obj, model_geodist))
}

plot_network_turnover_mod_compare <- function(mod1,
                                              mod2,
                                              this.network1,
                                              this.network2,
                                              network_type1,
                                              network_type2,
                                              this.effect,
                                              this.resp,
                                              label
){
  
  # prep cond effects mod 1
  
  mod_summary1 <- write.summ.table(mod1)
  model_geodist1 <- mod_summary1[rownames(mod_summary1) == "GeoDist",]
  
  if(network_type1 == "Obligate") {
    point_color1 <- "darkgreen"
    if(model_geodist1$Pgt0 >= 0.95){
      ribbon_color1 <- "darkgreen"
    } else if (model_geodist1$Pgt0 <= 0.05) {
      ribbon_color1 <- "darkgreen"
    } else {ribbon_color1 <- NA}
  }
  
  if(network_type1 == "Transient") {
    point_color1 <- "darkorange"
    if(model_geodist1$Pgt0 >= 0.95){
      ribbon_color1 <- "darkorange"
    } else if (model_geodist1$Pgt0 <= 0.05) {
      ribbon_color1 <- "darkorange"
    } else {ribbon_color1 <- NA}
  }
  
  # Extract the data from conditional_effects
  cond_effects_data1 <- conditional_effects(mod1, effects = this.effect, resp = this.resp, plot = FALSE)
  plot_data1 <- cond_effects_data1[[this.effect]]
  #browser()
  
  # now repeat for mod 2
  mod_summary2 <- write.summ.table(mod2)
  model_geodist2 <- mod_summary2[rownames(mod_summary2) == "GeoDist",]
  
  if(network_type2 == "Obligate") {
    point_color2 <- "darkgreen"
    if(model_geodist2$Pgt0 >= 0.95){
      ribbon_color2 <- "darkgreen"
    } else if (model_geodist2$Pgt0 <= 0.05) {
      ribbon_color2 <- "darkgreen"
    } else {ribbon_color2 <- NA}
  }
  
  if(network_type2 == "Transient") {
    point_color2 <- "darkorange"
    if(model_geodist2$Pgt0 >= 0.95){
      ribbon_color2 <- "darkorange"
    } else if (model_geodist2$Pgt0 <= 0.05) {
      ribbon_color2 <- "darkorange"
    } else {ribbon_color2 <- NA}
  }
  
  # Extract the data from conditional_effects
  cond_effects_data2 <- conditional_effects(mod2, effects = this.effect, resp = this.resp, plot = FALSE)
  plot_data2 <- cond_effects_data2[[this.effect]]
  
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data1, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill = ribbon_color1,
                color = point_color1, linetype='dotted') +
    geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=ribbon_color1, 
                color = point_color1, linetype='dashed') +
    geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=ribbon_color1,
                color = point_color1, linetype='solid') +
      # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(data=plot_data2, aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill = ribbon_color2,
                color = point_color2, linetype='dotted') +
    geom_ribbon(data=plot_data2, aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=ribbon_color2, 
                color = point_color2, linetype='dashed') +
    geom_ribbon(data=plot_data2, aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=ribbon_color2,
                color = point_color2, linetype='solid') +
    #Add line for the estimates
    geom_line(data = plot_data1, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data1, color = point_color1, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add points for original data
    geom_point(data = this.network1, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = point_color1, alpha = 0.9,color="black", pch=21, cex=3) +
    #Add line for the estimates
    geom_line(data = plot_data2, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data2, color = point_color2, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add points for original data
    geom_point(data = this.network2, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = point_color2, alpha = 0.9,color="black", pch=21, cex=3) +
    theme_classic()  +
    labs(x = "Geographic Distance (km)", y = label,
         fill = "Credible Interval") +
    theme_classic() +
    ylim(-0.15,1) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  combined_mods <- bind_rows(model_geodist1, model_geodist2)
  
  rownames(combined_mods) <- c(paste(network_type1,label), paste(network_type2,label))
  #TODO add together model geodist objs for output
  return(list(plot_obj, combined_mods))
}
