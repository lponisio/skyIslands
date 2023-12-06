standardizeVars <- function(spec.data, vars, key){
  ##  center all of the x variables, need to use unique values to avoid
  ##  repetition by the number of specimens

  unique.site.vals <- spec.data[spec.data$Weights == 1,
                               c("Site", key, vars)]
  unique.site.vals[, vars] <- apply(unique.site.vals[, vars], 2, standardize)
  print(dim(spec.data))
  spec.data[, vars] <- NULL
  spec.data <- merge(spec.data, unique.site.vals, all.x=TRUE)
  print(dim(spec.data))
  #layout(matrix(1:(2*round(length(vars)/2)), nrow=2))
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
                        is.na(spec.net$Apidae) |
                        spec.net$Genus != "Bombus"] <- 0
  spec.net$WeightsMicrobe <- 1
  spec.net$WeightsMicrobe[spec.net$Apidae == 0 |
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