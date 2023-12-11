standardizeVars <- function(spec.data, vars, key){
  ##  center all of the x variables, need to use unique values to avoid
  ##  repetition by the number of specimens
  unique.site.vals <-  unique(spec.data[,c("Site", key, vars)])
  unique.site.vals[, vars] <- apply(unique.site.vals[, vars], 2, standardize)
  print("Dimensions of the data before merging the standardize data")
  print(dim(spec.data))
  spec.data[, vars] <- NULL
  spec.data <- merge(spec.data, unique.site.vals, all.x=TRUE)
  print("Dimensions of the data after merging the standardize data")
  print(dim(spec.data))
  layout(matrix(1:(2*round(length(vars)/2)), nrow=2))
  for(var in vars){
    hist(unique.site.vals[, var], main=var)
  }
                
  return(spec.data)
}

prepParasiteWeights <- function(){
  ## create a dummy varaible "WeightPar" for the parasite data. The
  ## intention is to keep stan from dropping data for site-level models,
  ## but weight is 0 for parasite models.
  spec.net$WeightsPar <- 1
  spec.net$WeightsPar[spec.net$Apidae == 0 |
                        is.na(spec.net$Apidae)] <- 0
  ## stan drops all NA data, so can set AnyParasite to 0 with WeightsPar
  ## to keep it in the models
  spec.net$ParasitePresence[is.na(spec.net$ParasitePresence)] <- 0
  spec.net$CrithidiaBombi[is.na(spec.net$CrithidiaBombi)] <- 0
  spec.net$CrithidiaPresence[is.na(spec.net$CrithidiaPresence)] <- 0
  spec.net$ApicystisSpp[is.na(spec.net$ApicystisSpp)] <- 0
  spec.net$Year <- as.factor(spec.net$Year)
  return(spec.net)
}


prepDataSEM <-
  function(spec.data,#individual level specimen data
           variables.to.log = NULL, #variables to be logged
           variables.to.log.1 = NULL, #variables to be logged + 1
           vars_yearsr = NULL,#variables to standardize at year site sampling round level 
           vars_sp = NULL) #variables to standardize at year site sampling round at the species level
  {
    ## Function for making the SEM weights and standarizing variables.
    spec.data <- spec.data[order(spec.data$Site), ]
    
    ## create a dummy variable "Weight" to deal with the data sets being at
    ## different levels to get around the issue of having to pass in one
    ## data set into brms
    spec.data$YearSR <-
      paste(spec.data$Year, spec.data$SampleRound, sep = ";")
    spec.data$YearSRGenusSpecies <-
      paste(spec.data$YearSR, spec.data$GenusSpecies, sep = ";")
    
    ## will need to modify when we have multiple years
    print("Number of unique site, year, sampling round combinations")
    print(length(unique(paste(spec.data$Site, spec.data$YearSR))))
    spec.data <- makeDataMultiLevel(spec.data, "Site", "YearSR")
    print("Number of individuals with Weights == 1, should be the same as above")
    print(sum(spec.data$Weights))
    
    if(!is.null(variables.to.log)){
      spec.data[, variables.to.log] <-
        log(spec.data[, variables.to.log])
    }
    if(!is.null(variables.to.log.1)){
      spec.data[, variables.to.log.1] <-
        log(spec.data[, variables.to.log.1]+ 1)
    }
    ##  center all of the x variables, need to use unique values to avoid
    ##  repetition by the number of specimens
    
    if(!is.null(vars_yearsr)){
      print("Standardizing variables with year, sampling round, site combinations")
      spec.data <- standardizeVars(spec.data, vars_yearsr, "YearSR")
    }
    if(!is.null(vars_sp)){
      print("Standardizing variables with year, sampling round,site, individual species")
      spec.data <-
        standardizeVars(spec.data, vars_sp, "YearSRGenusSpecies")
    }
    
    ## create a dumby varaible "WeightPar" for the parasite data. The
    ## original intention was to keep stan from dropping data for
    ## site-level models, but weight is 0 for parasite models.
    print("Number of successful parasite screenings")
    print(sum(spec.data$Apidae, na.rm = TRUE))
    spec.data <- prepParasiteWeights()
    print("Number of of individuals with WeightsPar == 1, should be the same as above")
    print(sum(spec.data$WeightsPar))
    return(spec.data)
  }