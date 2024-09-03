## This function creates and runs the models for the parasites. 
runCombinedParasiteModels <- function(spec.data,## data
                                      species.group,## genus of bee group
                                      parasites =
                                          c("CrithidiaPresence",
                                            "ApicystisSpp"),
                                      ## name of the parasites for the model
                                      xvars,## explanatory variables for the parasite model
                                      ncores=ncores, ## number of cores
                                      iter = 10^4,
                                      chains = 1,
                                      thin=1,
                                      init=0,
                                      data2 = NULL, ## Data for the Phylogeny
                                      SEM = TRUE,
                                      neg.binomial = FALSE,
                                      site.lat,
                                      xvar.name,
                                      ...){
    ## Create a list with the formulas for the different parasites models
    bf.parasite.formulas <- vector(mode="list",
                                   length=length(parasites))
    names(bf.parasite.formulas) <- parasites
    ## Create the models of the parasites using the variables provided in xvars
    if(!neg.binomial){
        print("binomial")
        for(parasite in parasites){
            formula.parasite  <- as.formula(paste(
                paste(parasite, "| subset(WeightsPar) + trials(1)"),
                paste(xvars,
                      collapse=" + "),
                sep=" ~ "))
            bf.parasite.formulas[[parasite]] <-  bf(formula.parasite,
                                                    family="beta_binomial")

            if(species.group == "bombus"){
                freq.formula <- as.formula(paste(
                    parasite,
                    paste(xvars[-length(xvars)], ## last xvar must by the phylogeny!!!!
                          collapse=" + "),
                    sep=" ~ "))
            } else{
                freq.formula <- as.formula(paste(
                    parasite,
                    paste(xvars, 
                          collapse=" + "),
                    sep=" ~ "))
            }
            
            run_plot_freq_model_diagnostics(
                freq.formula,
                this_data=spec.data[spec.data$WeightsPar == 1,],
                this_family="bernoulli",
                site.lat=paste(site.lat, xvar.name, sep="_"),
                species.group=species.group)

            freq.mod <- glmer(freq.formula, family="binomial",
                              data=spec.data[spec.data$WeightsPar ==
                                             1,],
                              glmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 100000)))
            print(summary(freq.mod))
            print(vif(freq.mod))

        }
    } else{
        print("negbinomial")
        for(parasite in parasites){
            formula.parasite  <- as.formula(paste(
                paste(paste0("Sp", parasite), "|  subset(WeightsSp)"),
                paste(c(xvars,"offset(SpScreened)"),
                      collapse=" + "),
                sep=" ~ "))
            bf.parasite.formulas[[parasite]] <-  bf(formula.parasite,
                                                    family="negbinomial")
            freq.formula <- as.formula(paste(
                paste(paste0("Sp", parasite)),
                paste(xvars[-length(xvars)], ## last xvar must by the phylogeny!!!!
                      collapse=" + "),
                sep=" ~ "))
            
            run_plot_freq_model_diagnostics(
                freq.formula,
                this_data=spec.data[spec.data$WeightsSp == 1,],
                this_family="negbinomial",
                site.lat=paste(site.lat, xvar.name, sep="_"),
                offset="SpScreened",
                species.group=species.group)

            freq.mod <- glmer.nb(freq.formula,
                                 data=spec.data[spec.data$WeightsSp ==
                                                1,],
                                 glmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 100000)))
            print(summary(freq.mod))
            print(vif(freq.mod))
        }}

    if(SEM){
        ## When there are two parasites or 1 parasite create a parasite model for each. 
        ## Select the bee abundance based on the species.group.
        if(length(parasites) == 2){
            ## Bombus or apis
            if(species.group == "bombus" | species.group == "apis"){
                print("Bombus")
                bform <- bf.fabund + bf.fdiv +
                    bf.bombusabund + bf.babund + bf.HBabund +
                    bf.bdiv  +    
                    bf.parasite.formulas[[1]]+
                    bf.parasite.formulas[[2]] +
                    set_rescor(FALSE)
            } ## Other bees
            else if (species.group != "bombus" & species.group != "apis"){
                print("Other")
                ## only all bee abundance
                bform <- bf.fabund + bf.fdiv +
                    bf.babund  +
                    bf.bdiv  +    
                    bf.parasite.formulas[[1]]+
                    bf.parasite.formulas[[2]] +
                    set_rescor(FALSE)
            }
        }else  if(length(parasites) == 1){
            ## Bombus or apis
            if(species.group == "bombus" | species.group == "apis"){
                print("Bombus")
                bform <- bf.fabund + bf.fdiv +
                    bf.bombusabund + bf.babund + bf.HBabund +
                    bf.bdiv  +    
                    bf.parasite.formulas[[1]]+
                    set_rescor(FALSE)
            } ## Other bees
            else if (species.group != "bombus" & species.group != "apis"){
                print("Other")
                ## only all bee abundance
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
                        control = list(adapt_delta = 0.999999999999,  max_treedepth=12),
                        save_pars = save_pars(all = TRUE),
                        data2 = data2,
                        drop_unused_levels = TRUE,
                        ...)
    ## Create a table with the results. 
    write.ms.table(fit.parasite,
                   sprintf("parasitism_%s_%s_%s_%s",
                           species.group, paste(parasites,
                                                collapse=""), site.lat,
                           xvar.name))
    ## Calculate r2 values 
    r2 <- bayes_R2(fit.parasite)
    print(round(r2, 2))
    ## Save the model results as a rdata file
    save(fit.parasite, spec.data, r2,
         file=sprintf("saved/parasiteFit_%s_%s_%s_%s.Rdata",
                      species.group, paste(parasites, collapse=""),
                      site.lat,
                      xvar.name))
    ## Plot the residuals
    plot.res(fit.parasite,  sprintf("%s_%s_%s_%s",
                                    species.group, paste(parasites,
                                                         collapse=""),
                                    site.lat, xvar.name))

    loo.crithidia <- loo(fit.parasite, resp="CrithidiaPresence")
    loo.apicystis <- loo(fit.parasite, resp="ApicystisSpp")

    return(list(fit=fit.parasite, loo.crithidia= loo.crithidia,
                loo.apicystis=loo.apicystis,
                r2=r2))

}
