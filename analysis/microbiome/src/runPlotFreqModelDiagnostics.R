
## funciton to run frequentist version of brms models and plot diagnostics

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
    
    this_model_output <- MASS::zeroinfl(this_formula, data=this_data)
    diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
  }
  
  diagnostic.plots
}
