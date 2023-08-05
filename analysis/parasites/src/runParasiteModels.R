
runParasiteModels <- function(spec.data,
                              species.group, parasite,
                              xvars){

    formula.parasite  <- as.formula(paste(
        paste(parasite, "| weights(WeightsPar)"),
        paste(xvars,
              collapse=" + "),
        sep=" ~ "))


    formula.parasite.lmer  <- as.formula(paste(
        parasite,
        paste(xvars,
              collapse=" + "),
        sep=" ~ "))

  mod.test <- glmer(formula.parasite.lmer, family="binomial",
                    data=spec.data)
  print(vif(mod.test))

    bf.parasite <- bf(formula.parasite, family="bernoulli")
    bform.parasite <- bf.fabund + bf.fdiv +
        bf.babund + bf.bombusabund + bf.HBabund +
        bf.bdiv + bf.parasite +
        set_rescor(FALSE)

    fit.parasite <- brm(bform.parasite, spec.data,
                        cores=ncores,
                        iter = 10^4,
                        chains = 1,
                        thin=1,
                        init=0,
                        control = list(adapt_delta = 0.99))
    write.ms.table(fit.parasite,
                   sprintf("parasitism_%s_%s",
                           species.group, parasite))
    save(fit.parasite, spec.data,
         file=sprintf("saved/parasiteFit_%s_%s.Rdata",
                      species.group, parasite))
    return(fit.parasite)
}
