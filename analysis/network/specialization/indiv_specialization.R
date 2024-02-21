
### Purpose: Compare individual specialization within and between species using rbcL networks
           # Assess pollen fidelity using the rbcL networks

# Libraries 
  
library(vegan)
library(bipartite) 
library(tnet) 
library(lavaan)
library(tidyverse)


# Load individual networks _________________________________________________________________________

### Run scripts 4 + 6 + 7 under skyIslands/DataPrep to get networks

### Load rbcL networks for each site ______________

    # indivNet_rbcl should have the networks grouped by site already (script 7)
        ## I concatenated  UniqueIDs to include GenusSpecies info before building the networks for grep() later - not committed to script 7

# See format
print(indivNet_rbcl[["PL"]]) # pollinators are columns, plants are rows. Transpose:

# transpose and remove first 2 rows
indivNet.rbcl <- lapply(indivNet_rbcl,function(x){
  t(x)[-1,]
})
# print transposed data frame
indivNet.rbcl[[1]]
indivNet_rbcl[[1]]


# Subset pollinator data from each site - had to pare it down a lot because of error later on. Only including indivs with rbcL data. 
  # Will need to look into this and properly address later when we have all of the data.

poldata.rbcl <- subset(pol.abund, GenusSpecies  %in% c("Bombus.huntii",
                                               "Bombus.flavifrons",
                                               "Bombus.centralis",
                                               "Bombus.bifarius",
                                               "Apis.mellifera",
                                               "Anthophora.montana",
                                               "Agapostemon.angelicus") & Year == '2018' & SpParasitismRate != 'NA')

no.pollinator.rbcl<- aggregate(poldata.rbcl$Abundance, by=list(Category=poldata.rbcl$Site,poldata.rbcl$GenusSpecies), FUN=sum)
no.pollinator.rbcl$Site <- no.pollinator.rbcl$Category
no.pollinator.rbcl$GenusSpecies <- no.pollinator.rbcl$Group.2
no.pollinator.rbcl$Abundance <- no.pollinator.rbcl$x
no.pollinator.rbcl <- subset(no.pollinator.rbcl, select = -c(1:3))

 #list of unique pollinator species sampled at each site + column with the # of individuals of that species found at that site (abundance)
   
pollinators.rbcl<-list(subset(no.pollinator.rbcl,Site=="CH"),subset(no.pollinator.rbcl,Site=="HM"),subset(no.pollinator.rbcl,Site=="JC"),subset(no.pollinator.rbcl,Site=="MM"),
                  subset(no.pollinator.rbcl,Site=="PL"),#subset(no.pollinator.rbcl,Site=="RP"),
                  subset(no.pollinator.rbcl,Site=="SC"),subset(no.pollinator.rbcl,Site=="SM"))#,
                  #subset(no.pollinator.rbcl,Site=="VC")) 
                                                                      
# Subset pollinator species to have at least 5? (Tur et al., capped at 5 I think) individuals 
sp_names.rbcl<-list()
for(i in 1:length(pollinators.rbcl)){
sp_names.rbcl[[i]]<-as.character(subset(pollinators.rbcl[[i]],pollinators.rbcl[[i]]$Abundance>=5)$GenusSpecies) # 3 individuals right now, may change
}


print(sp_names.rbcl)

# Pollinator species subnetworks 

subnetworks.rbcl<-list()
for(n in 1:length(indivNet.rbcl)){
  m<-list()
  for (i in 1:length(sp_names.rbcl[[n]])){     
    m[[i]]<-indivNet.rbcl[[n]][grep(sp_names.rbcl[[n]][i],rownames(indivNet.rbcl[[n]])),]
  }
  subnetworks.rbcl[[n]]<-m
}

# Check how some look
subnetworks.rbcl[[1]]
subnetworks.rbcl[[5]]

subnetworks.rbcl = subnetworks.rbcl[-1] # CH site is acting strange, remove for this test

sp_names.rbcl = sp_names.rbcl[-1] # CH site is acting strange, remove for this test

# Individual Specialization ____________________________________________________________________

# Based on Roughgarden and Tur et al. papers

# Total niche width = within-individual component + between-individual component (TNW = WIC + BIC)

# Relative degree of individual specialization, or proportion of TNW explained by WIC: WIC/TNW
    #### A WIC/TNW of 1 means that individuals share the same niche breadth as their species. 
    #### A WIC/TNW of < 1 means there is individual specialization, with higher specialization as it approaches 0.

indices.rbcl<-list()
for (n in 1:length(subnetworks.rbcl)){
  results<-matrix(data=NA,nrow=length(subnetworks.rbcl[[n]]),ncol=4)
  rownames(results)<-sp_names.rbcl[[n]]
  colnames(results)<-c("WIC","TNW","BIC","WIC/TNW")
  for(i in 1:length(subnetworks.rbcl[[n]])){
# WIC
  results[i,1]<-sum((rowSums(subnetworks.rbcl[[n]][[i]])/sum(rowSums(subnetworks.rbcl[[n]][[i]]))) 
                        * 
                      (diversity(subnetworks.rbcl[[n]][[i]],"shannon"))
                    ) 

#### This includes calculating proportion pollen type per individual, which is what our network has. Fix!!  
# TNW   
  results[i,2]<- -sum((colSums(subnetworks.rbcl[[n]][[i]])/sum(rowSums(subnetworks.rbcl[[n]][[i]]))) 
                        * 
                      log((colSums(subnetworks.rbcl[[n]][[i]])/sum(rowSums(subnetworks.rbcl[[n]][[i]]))))
                      ) 
# BIC 
  results[i,3]<-results[i,2]-results[i,1] 
  
# WIC/TNW    
  results[i,4]<-results[i,1]/results[i,2]
  
  }
  indices.rbcl[[n]]<-round(results,2)
}

indices.rbcl # List with WIC/TNW values for species at all sites


#### Null networks ______________________________________________________________________________________________

# Create null i-sp networks to calculate null WIC/TNW for comparison to our empirical WIC/TNW (if lower, = high ind. specialization)
  ### Null assumes that all individuals have the same niche breadth as their species (WIC/TNW will be close to 1)

### Null i-sp subnetworks (using the Patefield algorithm)

spp2.rbcl<-list()
for (i in 1:length(pollinators.rbcl)){
  spp2.rbcl[[i]]<-as.character(subset(pollinators.rbcl[[i]],pollinators.rbcl[[i]]$no.ind>1)$GenusSpecies)
}
spp2.rbcl # List of species in each site with > 1 individual

subnetworks2.rbcl<-list()
for(n in 1:length(indivNet_rbcl)){
  m<-list()
  for (i in 1:length(spp2.rbcl[[n]])){
    m[[i]]<-indivNet_rbcl[[n]][grep(spp2.rbcl[[n]][i],rownames(indivNet_rbcl[[n]])),]
  }
  subnetworks2.rbcl[[n]]<-m
}
subnetworks2.rbcl # Species subnetworks for all species with > 1 individual sampled

null.networks.rbcl<-list()
for(n in 1:length(subnetworks2.rbcl)){
  NULL.subnetworks.rbcl<-list()
  for (i in 1:length(subnetworks2.rbcl[[n]])){
    NULL.subnetworks.rbcl[[i]]<-nullmodel(subnetworks2.rbcl[[n]][[i]],N=1000,method=1) # method 1 = Fixed row and column totals.
    for(j in 1:1000){
      colnames(NULL.subnetworks.rbcl[[i]][[j]])<-colnames(subnetworks.rbcl[[n]][[i]])
      rownames(NULL.subnetworks.rbcl[[i]][[j]])<-rownames(subnetworks.rbcl[[n]][[i]])
    }
  }
  null.networks.rbcl<-list()
  
  # Combine null i-sp subnetworks 
  
  for(j in 1:1000){
    m2<-list()
    null<-matrix()
    for (i in 1:length(NULL.subnetworks.rbcl)){
      m2[[i]]<-NULL.subnetworks.rbcl[[i]][[j]]
      null<-do.call(rbind,m2)
    }
    null.networks.rbcl[[j]]<-rbind(m1[[n]],null)
    null.networks.rbcl[[j]]<-null.networks.rbcl[[j]][order(rownames(null.networks.rbcl[[j]])),]
  }
  null.networks.rbcl[[n]]<-null.networks.rbcl
}

null.networks.rbcl[[1]] # List of 1,000 null i-sp networks for site 1 (will need to see how sites are assigned, alphabetical order?)
null.networks.rbcl[[2]] # List of 1,000 null i-sp networks for site 2
null.networks.rbcl[[3]] # List of 1,000 null i-sp networks for site 3
null.networks.rbcl[[4]] # List of 1,000 null i-sp networks for site 4
null.networks.rbcl[[5]] # List of 1,000 null i-sp networks for site 5
null.networks.rbcl[[6]] # List of 1,000 null i-sp networks for site 6
null.networks.rbcl[[7]] # List of 1,000 null i-sp networks for site 7
#null.networks.rbcl[[8]] # List of 1,000 null i-sp networks for site 8
#null.networks.rbcl[[9]] # List of 1,000 null i-sp networks for site 9

### Calculate null WIC/TNW

##### Mirror changes to the empirical calculations once formulas are fixed

null.wic.tnw.rbcl<-list()
for(d in 1:length(NULL.subnetworks.rbcl)){
  null.indices<-list()
  for (n in 1:length(NULL.subnetworks.rbcl[[d]])){
    results<-matrix(data=NA,nrow=length(NULL.subnetworks.rbcl[[d]][[n]]),ncol=4)
    colnames(results)<-c("WIC","TNW","BIC","WIC/TNW")
    for(i in 1:length(NULL.subnetworks.rbcl[[d]][[n]])){

   results[i,1]<-sum((rowSums(NULL.subnetworks.rbcl[[d]][[n]][[i]])/sum(rowSums(NULL.subnetworks.rbcl[[d]][[n]][[i]]))) * (diversity(NULL.subnetworks.rbcl[[d]][[n]][[i]],"shannon"))) # WIC
      
   results[i,2]<- -sum((colSums(NULL.subnetworks.rbcl[[d]][[n]][[i]])/sum(rowSums(NULL.subnetworks.rbcl[[d]][[n]][[i]]))) * log((colSums(NULL.subnetworks.rbcl[[d]][[n]][[i]])/sum(rowSums(NULL.subnetworks.rbcl[[d]][[n]][[i]]))))) # TNW
    
   results[i,3]<-results[i,2]-results[i,1] # BIC
      
   results[i,4]<-results[i,1]/results[i,2] # WIC/TNW
    }
    null.indices[[n]]<-results
  }
  null.wic.tnw.rbcl[[d]]<-null.indices
  names(null.wic.tnw.rbcl[[d]])<-sp_names.rbcl[[d]]
}

null.wic.tnw.rbcl # List with null WIC/TNW values


# Comparing empirical and null WIC/TNW __________________________________________________________________________________

# See distribution of null WIC/TNW for species at each site with histograms, mean, min
### Add t-tests to compare to empirical WIC/TNW??

# Site 1
  histnull1 <- par(mfrow=c(2,)) # Decide rows/columns based on # species at each site
  for(n in 1:length(null.wic.tnw.rbcl[[1]])){
    histnull1[[n]]<-hist(null.wic.tnw.rbcl[[1]][[n]][,4],xlab="WIC/TNW",main=sp_names[[1]][n])
  }
 
  meansnull1 <- par(mfrow=c(2,))
  for(n in 1:length(null.wic.tnw.rbcl[[1]])){
    meansnull1[[n]] <- mean(null.wic.tnw.rbcl[[1]][[n]][,4],xlab="WIC/TNW",main=sp_names[[1]][n])} # Empirical 95 % CI overlap???
  
  minsnull1 <- par(mfrow=c(2,))
  for(n in 1:length(null.wic.tnw.rbcl[[1]])){
    minsnull1[[n]] <- min(null.wic.tnw.rbcl[[1]][[n]][,4],xlab="WIC/TNW",main=sp_names[[1]][n])}
  
# Site 2 
  histnull2 <- par(mfrow=c(2,)) 
  for(n in 2:length(null.wic.tnw.rbcl[[2]])){
    histnull2[[n]]<-hist(null.wic.tnw.rbcl[[2]][[n]][,4],xlab="WIC/TNW",main=sp_names[[2]][n])
  }
  
  meansnull2 <- par(mfrow=c(2,))
  for(n in 2:length(null.wic.tnw.rbcl[[2]])){
    meansnull2[[n]] <- mean(null.wic.tnw.rbcl[[2]][[n]][,4],xlab="WIC/TNW",main=sp_names[[2]][n])} 
  
  minsnull2 <- par(mfrow=c(2,))
  for(n in 2:length(null.wic.tnw.rbcl[[2]])){
    minsnull2[[n]] <- min(null.wic.tnw.rbcl[[2]][[n]][,4],xlab="WIC/TNW",main=sp_names[[2]][n])}
  
# Site 3 
  histnull3 <- par(mfrow=c(2,)) 
  for(n in 3:length(null.wic.tnw.rbcl[[3]])){
    histnull3[[n]]<-hist(null.wic.tnw.rbcl[[3]][[n]][,4],xlab="WIC/TNW",main=sp_names[[3]][n])
  }
  
  meansnull3 <- par(mfrow=c(2,))
  for(n in 3:length(null.wic.tnw.rbcl[[3]])){
    meansnull3[[n]] <- mean(null.wic.tnw.rbcl[[3]][[n]][,4],xlab="WIC/TNW",main=sp_names[[3]][n])}
  
  minsnull3 <- par(mfrow=c(2,))
  for(n in 3:length(null.wic.tnw.rbcl[[3]])){
    minsnull3[[n]] <- min(null.wic.tnw.rbcl[[3]][[n]][,4],xlab="WIC/TNW",main=sp_names[[3]][n])}
  
# Site 4 
  histnull4 <- par(mfrow=c(2,))
  for(n in 4:length(null.wic.tnw.rbcl[[4]])){
    histnull4[[n]]<-hist(null.wic.tnw.rbcl[[4]][[n]][,4],xlab="WIC/TNW",main=sp_names[[4]][n])
  }
  
  meansnull4 <- par(mfrow=c(2,))
  for(n in 4:length(null.wic.tnw.rbcl[[4]])){
    meansnull4[[n]] <- mean(null.wic.tnw.rbcl[[4]][[n]][,4],xlab="WIC/TNW",main=sp_names[[4]][n])} 
  
  minsnull4 <- par(mfrow=c(2,))
  for(n in 4:length(null.wic.tnw.rbcl[[4]])){
    minsnull4[[n]] <- min(null.wic.tnw.rbcl[[4]][[n]][,4],xlab="WIC/TNW",main=sp_names[[4]][n])}
  
# Site 5 
  histnull5 <- par(mfrow=c(2,)) 
  for(n in 5:length(null.wic.tnw.rbcl[[5]])){
    histnull5[[n]]<-hist(null.wic.tnw.rbcl[[5]][[n]][,4],xlab="WIC/TNW",main=sp_names[[5]][n])
  }
  
  meansnull5 <- par(mfrow=c(2,))
  for(n in 5:length(null.wic.tnw.rbcl[[5]])){
    meansnull5[[n]] <- mean(null.wic.tnw.rbcl[[5]][[n]][,4],xlab="WIC/TNW",main=sp_names[[5]][n])}
  
  minsnull5 <- par(mfrow=c(2,))
  for(n in 5:length(null.wic.tnw.rbcl[[5]])){
    minsnull5[[n]] <- min(null.wic.tnw.rbcl[[5]][[n]][,4],xlab="WIC/TNW",main=sp_names[[5]][n])}
  
# Site 6 
  histnull6 <- par(mfrow=c(2,)) 
  for(n in 6:length(null.wic.tnw.rbcl[[6]])){
    histnull6[[n]]<-hist(null.wic.tnw.rbcl[[6]][[n]][,4],xlab="WIC/TNW",main=sp_names[[6]][n])
  }
  
  meansnull6 <- par(mfrow=c(2,))
  for(n in 6:length(null.wic.tnw.rbcl[[6]])){
    meansnull6[[n]] <- mean(null.wic.tnw.rbcl[[6]][[n]][,4],xlab="WIC/TNW",main=sp_names[[6]][n])}
  
  minsnull6 <- par(mfrow=c(2,))
  for(n in 6:length(null.wic.tnw.rbcl[[6]])){
    minsnull6[[n]] <- min(null.wic.tnw.rbcl[[6]][[n]][,4],xlab="WIC/TNW",main=sp_names[[6]][n])}
  
# Site 7 
  histnull7 <- par(mfrow=c(2,))
  for(n in 7:length(null.wic.tnw.rbcl[[7]])){
    histnull7[[n]]<-hist(null.wic.tnw.rbcl[[7]][[n]][,4],xlab="WIC/TNW",main=sp_names[[7]][n])
  }
  
  meansnull7 <- par(mfrow=c(2,))
  for(n in 7:length(null.wic.tnw.rbcl[[7]])){
    meansnull7[[n]] <- mean(null.wic.tnw.rbcl[[7]][[n]][,4],xlab="WIC/TNW",main=sp_names[[7]][n])}
  
  minsnull7 <- par(mfrow=c(2,))
  for(n in 7:length(null.wic.tnw.rbcl[[7]])){
    minsnull7[[n]] <- min(null.wic.tnw.rbcl[[7]][[n]][,4],xlab="WIC/TNW",main=sp_names[[7]][n])}
  
# Site 8
  histnull8 <- par(mfrow=c(2,)) 
  for(n in 8:length(null.wic.tnw.rbcl[[8]])){
    histnull8[[n]]<-hist(null.wic.tnw.rbcl[[8]][[n]][,4],xlab="WIC/TNW",main=sp_names[[8]][n])
  }
  
  meansnull8 <- par(mfrow=c(2,))
  for(n in 8:length(null.wic.tnw.rbcl[[8]])){
    meansnull8[[n]] <- mean(null.wic.tnw.rbcl[[8]][[n]][,4],xlab="WIC/TNW",main=sp_names[[8]][n])} 
  
  minsnull8 <- par(mfrow=c(2,))
  for(n in 8:length(null.wic.tnw.rbcl[[8]])){
    minsnull8[[n]] <- min(null.wic.tnw.rbcl[[8]][[n]][,4],xlab="WIC/TNW",main=sp_names[[8]][n])}
  
# Site 9 
  histnull9 <- par(mfrow=c(2,)) 
  for(n in 9:length(null.wic.tnw.rbcl[[9]])){
    histnull9[[n]]<-hist(null.wic.tnw.rbcl[[9]][[n]][,4],xlab="WIC/TNW",main=sp_names[[9]][n])
  }
  
  meansnull9 <- par(mfrow=c(2,))
  for(n in 9:length(null.wic.tnw.rbcl[[9]])){
    meansnull9[[n]] <- mean(null.wic.tnw.rbcl[[9]][[n]][,4],xlab="WIC/TNW",main=sp_names[[9]][n])} 
  minsnull9 <- par(mfrow=c(2,))
  for(n in 9:length(null.wic.tnw.rbcl[[9]])){
    minsnull9[[n]] <- min(null.wic.tnw.rbcl[[9]][[n]][,4],xlab="WIC/TNW",main=sp_names[[9]][n])}
  
# Proportional Similarity index (PSi) from Tur et al. ___________________________________________________________________

psi<-list()
for(n in 1:length(subnetworks.rbcl)){
  pij<-list()
  qj<-list()
  Qj<-list()
  for(i in 1:length(subnetworks.rbcl[[n]])){
    Qj[[i]]<-matrix(nrow=nrow(subnetworks.rbcl[[n]][[i]]),ncol=ncol(subnetworks.rbcl[[n]][[i]]))
    colnames(Qj[[i]])<-colnames(subnetworks.rbcl[[n]][[i]])
  }
  
  abs_pij_qj<-list()
  PSi<-list()
  for (i in 1:length(subnetworks.rbcl[[n]])){
    pij[[i]]<-subnetworks.rbcl[[n]][[i]]/rowSums(subnetworks.rbcl[[n]][[i]])
    qj[[i]]<-colSums(subnetworks.rbcl[[n]][[i]])/sum(rowSums(subnetworks.rbcl[[n]][[i]]))
    for(d in 1:nrow(Qj[[i]])){
      Qj[[i]][d,]<-qj[[i]]
    }
    abs_pij_qj[[i]]<-abs(pij[[i]]-Qj[[i]])
    PSi[[i]]<-1-0.5*(rowSums(abs_pij_qj[[i]]))
  }
  psi[[n]]<-PSi
}
psi # List with the PSi value for each individual
hist(unlist(psi), main="") # Histogram 


# Calculating intraspecific overlap _______________________________________________________________________________________

intrasp.overlap<-list()
for(n in 1:length(subnetworks.rbcl)){
  io<-list()
  IO<-matrix(NA,nrow=length(subnetworks.rbcl[[n]]))
  row.names(IO)<-sp_names[[n]]
  for(i in 1:length(subnetworks.rbcl[[n]])){
    io[[i]]<-degree_w((projecting_tm(subnetworks.rbcl[[n]][[i]],method="sum")),measure=c("degree","output"), type="out")
    IO[i,]<-sum(io[[i]][,3])/((ncol(subnetworks.rbcl[[n]][[i]])*(nrow(subnetworks.rbcl[[n]][[i]])-1))*nrow(subnetworks.rbcl[[n]][[i]]))
  }
  intrasp.overlap[[n]]<-IO
}
intrasp.overlap


# Pollen fidelity ___________________________________________________________________________________________________________









