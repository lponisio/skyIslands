## This function is to run the multilevel model using brms.
## The function needs a data set and to specify in the parameters which bee species
## groups and parasite to include as well as the x variables that we are using
runParasiteModels <- function(spec.data,
                              species.group, parasite,
                              xvars){
## Create the parasite formula. y is the parasite weight, x variables are specified 
## in 1amplification_dilution.R
  
  ## formula.parasite  <- as.formula(paste(
  ##   paste(parasite, "| weights(WeightsPar)"),
  ##   paste(xvars,
  ##         collapse=" + "),
  ##   sep=" ~ "))

## Drop the specimens that were no screened for parasites
  data.par <- spec.data[spec.data$WeightsPar == 1,]
  
  formula.parasite  <- as.formula(paste(
    parasite,
    paste(xvars,
          collapse=" + "),
    sep=" ~ "))
## Test the formula with GLMM
  mod.test <- glmer(formula.parasite, family="binomial",
                    data=data.par)
## Check for variance inflation factor to determine if there is colinearity 
## between variables  
  print(vif(mod.test))
  print(summary(mod.test))
## Brms parasite formula
  bf.parasite <- bf(formula.parasite, family="bernoulli")
## Fit the brms formula to the model.  
  fit.parasite <- brm(bf.parasite, data.par,
                      cores=ncores,
                      iter = 10^4,
                      chains = 2,
                      thin=1,
                      init=0,
                      control = list(adapt_delta = 0.99))
## Create a csv table   
  write.ms.table(fit.parasite,
                 sprintf("parasitism_%s_%s",
                         species.group, parasite))
## Calculate the r2 value. This determines how much of the variation is explained by
## out x variables
  r2 <- bayes_R2(fit.parasite)
  print(r2)
## Create a r file with the output
  save(fit.parasite, data.par, r2,
       file=sprintf("saved/parasiteFit_%s_%s.Rdata",
                    species.group, parasite))
## The function returns the results of the brms model 
  return(fit.parasite)
}
