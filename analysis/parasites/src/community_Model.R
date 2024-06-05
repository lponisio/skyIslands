## This function creates a community model and calculates it R2 and diagnostics plots
runCommunityModels <- function(bform.model, ## brms formulas
                               spec.data,## data
                                      ncores, ## number of cores
                                      iter = 10^4,
                                      chains = 1,
                                      thin=1,
                                      init=0, data2 = NULL){
  
  fit.community <- brm(bform.community, spec.data,
                       cores = ncores,
                       iter = iter,
                       chains = chains,
                       thin= thin,
                       control = list(adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
  write.ms.table(fit.community,
                 sprintf("parasitism_%s_%s",
                         species.group="all", parasite="none"))
  
  r2loo <- loo_R2(fit.community)
  r2 <- rstantools::bayes_R2(fit.community)
  save(fit.community, spec.net, r2, spec.orig,
       file="saved/communityFit.Rdata")
  plot.res(fit.community,  sprintf("%s",
                                  "community"))
  
  return(list(fit= fit.community, r2=r2))
  
}