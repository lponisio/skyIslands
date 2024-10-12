

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
