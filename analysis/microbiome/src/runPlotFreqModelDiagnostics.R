
## funciton to run frequentist version of brms models and plot diagnostics

#currently works with all families included in ?family, as well as zero inflated neg binomial, negbinomial, Gamma, inverse gaussian, poisson, quasibinomial
# will need to add other families as needed 



run_plot_freq_model_diagnostics <- function(this_formula, #brms model formula
                            this_data, #data frame, subsetted to correct weights!
                            num_chains=1, 
                            num_iter=10000, 
                            this_family #model family
                            ){
  if (this_family == 'gaussian'){
  #run model
  this_model_output <- brms::brm(this_formula,
                         data = this_data, 
                         chains = num_chains, 
                         iter = num_iter, family=this_family)

  
  this_model_output
  # return a list of single plots
  diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))} else if (this_family=='negbinomial') {
    
    this_model_output <- glmer.nb(this_formula, data=this_data)
    diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
  } else if (this_family=='zero_inflated_negbinomial') {
    
    this_model_output <- glmmadmb(this_formula, data=this_data, zeroInflation=TRUE, family='nbinom')
    diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
  }
  else if (this_family=='hurdle_gamma') {
    
    this_model_output <- glmer(this_formula, data=this_data, family=Gamma)
    diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
  } else if (this_family=='negbinomial') {
    this_model_output <- glmmadmb(this_formula, data=this_data, zeroInflation=FALSE, family='nbinom')
    diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
  } else if (this_family == 'Gamma'){
    #run model
    this_model_output <- glmer(this_formula, data=this_data, family=Gamma(link=log))
   
    # return a list of single plots
    diagnostic.plots <- plot(check_model(this_model_output, 
                                         panel = TRUE))
  } else if (this_family == 'inverse.gaussian'){
    #run model
    this_model_output <- brms::brm(this_formula,
                                   data = this_data, 
                                   chains = num_chains, 
                                   iter = num_iter, family=this_family)
    # return a list of single plots
    diagnostic.plots <- plot(check_model(this_model_output, 
                                         panel = TRUE))
  } else if (this_family == 'poisson'){
    #run model
    this_model_output <- brms::brm(this_formula,
                                   data = this_data, 
                                   chains = num_chains, 
                                   iter = num_iter, family=this_family)
    # return a list of single plots
    diagnostic.plots <- plot(check_model(this_model_output, 
                                         panel = TRUE), type = "discrete_dots")
  } 
  else if (this_family == 'quasibinomial'){
    #run model
    this_model_output <- brms::brm(this_formula,
                                   data = this_data, 
                                   chains = num_chains, 
                                   iter = num_iter, family=this_family)
    # return a list of single plots
    diagnostic.plots <- plot(check_model(this_model_output, 
                                         panel = TRUE))
  } 
  
  diagnostic.plots
}
