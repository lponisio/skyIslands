plot_model_condeff_compare <- function(model.a=fit.microbe.bombus,
                                       model.b=fit.microbe.melissodes,
                                       this.effect='BeeDiversity',
                                       this.resp.a='PDobligate', ## TODO: potentially update to both be PD obligate log?
                                       this.resp.b='PDobligatelog',
                                       point.data.a=bombus.obligate,
                                       point.data.b=melissodes.obligate,
                                       axis.breaks=axis.bee.div,
                                       axis.labs=labs.bee.div,
                                       xlabel="Bee Diversity",
                                       ylabel="Obligate Microbe PD",
                                       mod1color='navy',
                                       mod2color='coral',
                                       fill.a=FALSE,
                                       fill.b=FALSE
){
  
  ## PD obligate ~ Bee Diversity
  
  # Extract the data from conditional_effects
  cond_effects_data_a <- conditional_effects(model.a, effects = this.effect, resp = this.resp.a, plot = FALSE)
  plot_data_a <- cond_effects_data_a[[paste(this.resp.a,".",this.resp.a, "_",this.effect, sep='')]]
  #browser()
  
  ## fixing col names TODO: should fix the col names in prep
  if ('PD.obligate' %in% colnames(plot_data_a)){
    plot_data_a$PDobligate <- plot_data_a$PD.obligate
  }
  if ('PD.obligate.log' %in% colnames(plot_data_a)){
    plot_data_a$PDobligatelog <- plot_data_a$PD.obligate.log
  }
  if ('PD.transient' %in% colnames(plot_data_a)){
    plot_data_a$PDtransient <- plot_data_a$PD.transient
  }
  if ('PD.transient.log' %in% colnames(plot_data_a)){
    plot_data_a$PDtransientlog <- plot_data_a$PD.transient.log
  }
  
  # Extract the data from conditional_effects
  cond_effects_data_b <- conditional_effects(model.b, effects = this.effect, resp = this.resp.b, plot = FALSE)
  plot_data_b <- cond_effects_data_b[[paste(this.resp.b, ".", this.resp.b,"_",this.effect, sep='')]]
  
  ## fixing col names TODO: should fix the col names in prep
  if ('PD.obligate' %in% colnames(plot_data_b)){
    plot_data_b$PDobligate <- plot_data_b$PD.obligate
  }
  if ('PD.obligate.log' %in% colnames(plot_data_b)){
    plot_data_b$PDobligatelog <- plot_data_b$PD.obligate.log
  }
  if ('PD.transient' %in% colnames(plot_data_b)){
    plot_data_b$PDtransient <- plot_data_b$PD.transient
  }
  if ('PD.transient.log' %in% colnames(plot_data_b)){
    plot_data_b$PDtransientlog <- plot_data_b$PD.transient.log
  }
  
  ## model a fill
  if (fill.a==TRUE){
    mod1fill=mod1color
  } else {mod1fill=NA}
  
  ## model a fill
  if (fill.b==TRUE){
    mod2fill=mod2color
  } else {mod2fill=NA}
  
  
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data_a, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill=mod1fill,
                color = mod1color, linetype='dotted') +
    geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=mod1fill, 
                color = mod1color, linetype='dashed') +
    geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=mod1fill,
                color = mod1color, linetype='solid') +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(data=plot_data_b, aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill=mod2fill, color = mod2color, linetype='dotted') +
    geom_ribbon(data=plot_data_b, aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                                      ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, fill=mod2fill, color = mod2color, linetype='dashed') +
    geom_ribbon(data=plot_data_b, aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                                      ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, fill=mod2fill, color = mod2color, linetype='solid') +
    #Add line for the estimates
    geom_line(data = plot_data_a, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data_a, color = mod1color, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data_b, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add line for the estimates
    geom_line(data=plot_data_b, linewidth=2, aes(x = .data[[this.effect]],  y = .data$estimate__), color = mod2color) +
    # Add points for original data
    geom_point(data = point.data.a, aes(x = .data[[this.effect]], y = .data[[this.resp.a]]),
               fill = mod1color, alpha = 0.9,color="black", pch=21, cex=3) +
    # Add points for original data
    geom_point(data = point.data.b, aes(x = .data[[this.effect]], y = .data[[this.resp.b]]),
               fill = mod2color, alpha = 0.9, color="black", pch=21, cex=3) +
    # Labels and theme
    labs(x = xlabel, y = ylabel) +
    scale_x_continuous(breaks = axis.breaks, labels = axis.labs) +
    theme_classic() +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  plot_obj
  #browser()
}

plot_model_condeff_single <- function(model=fit.microbe.bombus,
                                       this.effect='BeeDiversity',
                                       this.resp='PDobligate', ## TODO: potentially update to both be PD obligate log?
                                       point.data=bombus.obligate,
                                       axis.breaks=axis.bee.div,
                                       axis.labs=labs.bee.div,
                                       xlabel="Bee Diversity",
                                       ylabel="Obligate Microbe PD",
                                       mod1color='navy'
){
  
  ## PD obligate ~ Bee Diversity
  
  # Extract the data from conditional_effects
  cond_effects_data <- conditional_effects(model, effects = this.effect, resp = this.resp, plot = FALSE)
  plot_data <- cond_effects_data[[paste(this.resp,".",this.resp, "_",this.effect, sep='')]]
  
  ## fixing col names TODO: should fix the col names in prep
  if ('PD.obligate' %in% colnames(plot_data)){
    plot_data$PDobligate <- plot_data$PD.obligate
  }
  if ('PD.obligate.log' %in% colnames(plot_data)){
    plot_data$PDobligatelog <- plot_data$PD.obligate.log
  }
  if ('PD.transient' %in% colnames(plot_data)){
    plot_data$PDtransient <- plot_data$PD.transient
  }
  if ('PD.transient.log' %in% colnames(plot_data)){
    plot_data$PDtransientlog <- plot_data$PD.transient.log
  }
  
  
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill=mod1color,
                color = mod1color, linetype='dotted') +
    geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=mod1color, 
                color = mod1color, linetype='dashed') +
    geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=mod1color,
                color = mod1color, linetype='solid') +
    #Add line for the estimates
    geom_line(data = plot_data, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data, color = mod1color, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add points for original data
    geom_point(data = point.data, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = mod1color, alpha = 0.9,color="black", pch=21, cex=3) +
    # Labels and theme
    labs(x = xlabel, y = ylabel) +
    scale_x_continuous(breaks = axis.breaks, labels = axis.labs) +
    theme_classic() +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  plot_obj
}

