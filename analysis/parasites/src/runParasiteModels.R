
runParasiteModels <- function(spec.data,
                              species.group, parasite,
                              xvars){

  ## formula.parasite  <- as.formula(paste(
  ##   paste(parasite, "| weights(WeightsPar)"),
  ##   paste(xvars,
  ##         collapse=" + "),
  ##   sep=" ~ "))


  data.par <- spec.data[spec.data$WeightsPar == 1,]
  
  formula.parasite  <- as.formula(paste(
    parasite,
    paste(xvars,
          collapse=" + "),
    sep=" ~ "))

  mod.test <- glmer(formula.parasite, family="binomial",
                    data=data.par)
  print(vif(mod.test))
  print(summary(mod.test))

  bf.parasite <- bf(formula.parasite, family="bernoulli")
  
  fit.parasite <- brm(bf.parasite, data.par,
                      cores=ncores,
                      iter = 10^4,
                      chains = 2,
                      thin=1,
                      ## init=0,
                      control = list(adapt_delta = 0.99))
  
  write.ms.table(fit.parasite,
                 sprintf("parasitism_%s_%s",
                         species.group, parasite))

  r2 <- loo_R2(fit.parasite)
  print(r2)

  save(fit.parasite, data.par, r2,
       file=sprintf("saved/parasiteFit_%s_%s.Rdata",
                    species.group, parasite))
  return(fit.parasite)
}
