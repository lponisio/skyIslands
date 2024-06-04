## This function creates and runs the models for the parasites. 
runCombinedParasiteModels <- function(spec.data,## data
                                      species.group,## genus of bee group
                                      parasites,  ## name of the parasites for the model
                                      xvars,## explanatory variables for the parasite model
                                      ncores, ## number of cores
                                      iter = 10^4,
                                      chains = 1,
                                      thin=1,
                                      init=0, data2 = NULL,
                                      SEM = TRUE,
                                      neg.binomial = FALSE){
  ## Create a list with the formulas for the different parasites models
  bf.parasite.formulas <- vector(mode="list",
                                 length=length(parasites))
  names(bf.parasite.formulas) <- parasites
  ## Create the models of the parasites using the variables provided in xvars
  if(!neg.binomial){
    print("binomial")
    for(parasite in parasites){
      formula.parasite  <- as.formula(paste(
        paste(parasite, "| weights(WeightsPar)"),
        paste(xvars,
              collapse=" + "),
        sep=" ~ "))
      bf.parasite.formulas[[parasite]] <-  bf(formula.parasite,
                                              family="bernoulli")
    }
  } else{
    print("negbinomial")
    for(parasite in parasites){
      formula.parasite  <- as.formula(paste(
        paste(paste0("Sp", parasite), "| weights(WeightsSp)"),
        paste(c(xvars,"offset(SpScreened)"),
              collapse=" + "),
        sep=" ~ "))
      bf.parasite.formulas[[parasite]] <-  bf(formula.parasite,
                                              family="negbinomial") 
    }}

  if(SEM){
    ## When there are two parasites or 1 parasite create a parasite model for each. 
    ## Select the bee abundance based on the species.group.
    if(length(parasites) == 2){
      ## Bombus
      if(species.group == "bombus"){
        print("Bombus")
        bform <- bf.fabund + bf.fdiv +
          bf.bombusabund +
          bf.bdiv  +    
          bf.parasite.formulas[[1]]+
          bf.parasite.formulas[[2]] +
          set_rescor(FALSE)
      } ## Apis
      else if (species.group == "apis"){
        print("Apis")
        bform <- bf.fabund + bf.fdiv +
          bf.HBabund +
          bf.bdiv  +    
          bf.parasite.formulas[[1]]+
          bf.parasite.formulas[[2]] +
          set_rescor(FALSE)
      } ## Other bees
      else if (species.group != "bombus" & species.group != "apis"){
        print("Other")
        bform <- bf.fabund + bf.fdiv +
          bf.babund +
          bf.bdiv  +    
          bf.parasite.formulas[[1]]+
          bf.parasite.formulas[[2]] +
          set_rescor(FALSE)
      }
    }else  if(length(parasites) == 1){
      ## Bombus
      if(species.group == "bombus"){
        print("Bombus")
        bform <- bf.fabund + bf.fdiv +
          bf.bombusabund +
          bf.bdiv  +    
          bf.parasite.formulas[[1]]+
          set_rescor(FALSE)
      } ## Apis
      else if (species.group == "apis"){
        print("Apis")
        bform <- bf.fabund + bf.fdiv +
          bf.HBabund +
          bf.bdiv  +    
          bf.parasite.formulas[[1]]+
          set_rescor(FALSE)
      } ## Other bees
      else if (species.group != "bombus" & species.group != "apis"){
        print("Other")
        bform <- bf.fabund + bf.fdiv +
          bf.babund 
        bf.bdiv  +    
          bf.parasite.formulas[[1]]+
          set_rescor(FALSE)
      }
    }
  } else {
    ## When there are two parasites or 1 parasite create a parasite model for each. 
    ## Select the bee abundance based on the species.group.
    if(length(parasites) == 2){
      bform <- 
        bf.parasite.formulas[[1]]+
        bf.parasite.formulas[[2]] +
        set_rescor(FALSE)      
    }else  if(length(parasites) == 1){
      bform <- 
        bf.parasite.formulas[[1]] +
        set_rescor(FALSE)
    }
  }
  ## Fit brms model to the complete model

  fit.parasite <- brm(bform, spec.data,
                      cores= ncores,
                      iter = iter,
                      chains = chains,
                      thin=thin,
                      init=init,
                      control = list(adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE), data2 = data2)
  ## Create a table with the results. 
  write.ms.table(fit.parasite,
                 sprintf("parasitism_%s_%s",
                         species.group, paste(parasites, collapse="")))
  ## Calculate r2 values 
  r2 <- bayes_R2(fit.parasite)
  print(round(r2, 2))
  ## Save the model results as a rdata file
  save(fit.parasite, spec.data, r2,
       file=sprintf("saved/parasiteFit_%s_%s.Rdata",
                    species.group, paste(parasites, collapse="")))
  ## Plot the residuals
  plot.res(fit.parasite,  sprintf("%s_%s",
                                  species.group, paste(parasites, collapse="")))

  return(list(fit=fit.parasite, r2=r2))
}
