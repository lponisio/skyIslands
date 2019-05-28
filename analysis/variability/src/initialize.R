library(bipartite, quietly = TRUE)
library(lme4, quietly = TRUE)
library(lmerTest, quietly = TRUE)
source('src/misc.R')
source('src/calcPca.R')
source('src/calcSpec.R')


load('../../data/nets.Rdata')
load('../../data/spec.Rdata')

save.path <- 'saved'

getNetData <- function(nets){
    sites <- sapply(strsplit(nets, "\\."),
                       function(x) x[1])
    years <-  sapply(strsplit(nets, "\\."),
                    function(x) x[2])
    SR <-  sapply(strsplit(nets, "\\."),
                        function(x) x[3])
    return(data.frame(Site=sites,
                      Year=years,
                      SR=SR))
}
