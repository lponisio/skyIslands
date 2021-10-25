mainDir <- "../skyIslands/data"
subDir <- "networks"

if (file.exists(subDir)){
} else {
    dir.create(file.path(mainDir, subDir),
               showWarnings = FALSE)
}
save.dir <- file.path(mainDir, subDir)

source('../skyIslands/dataPrep/src/illumSplit.R')
source('../skyIslands/dataPrep/src/calcFuncUniqOrig.R')
source('../skyIslands/dataPrep/src/misc.R')
load('../skyIslands/data/spec_RBCL_16s.Rdata')
