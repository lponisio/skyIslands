#!/usr/bin/env bash

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github

Rscript dataPrep/dataPrep.R

##***************************************************************
## turnover
##***************************************************************

Rscript analysis/turnover/1betalink.R "Yr" "Plant"
Rscript analysis/turnover/1betalink.R "Yr" "Parasite"

Rscript analysis/turnover/1betalink.R "YrSR" "Plant"
Rscript analysis/turnover/1betalink.R "YrSR" "Parasite"

Rscript analysis/turnover/2pca.R "Yr" "Plant" "lower.level"
Rscript analysis/turnover/2pca.R "Yr" "Plant" "higher.level"

Rscript analysis/turnover/2pca.R "Yr" "Parasite" "lower.level"
Rscript analysis/turnover/2pca.R "Yr" "Parasite" "higher.level"


Rscript analysis/turnover/2pca.R "YrSR" "Plant" "lower.level"
Rscript analysis/turnover/2pca.R "YrSR" "Plant" "higher.level"

Rscript analysis/turnover/2pca.R "YrSR" "Parasite" "lower.level"
Rscript analysis/turnover/2pca.R "YrSR" "Parasite" "higher.level"

Rscript analysis/turnover/2pca.R "YrSR" "Parasite" "lower.level"
Rscript analysis/turnover/2pca.R "YrSR" "Parasite" "higher.level"


Rscript analysis/network/1metrics.R "YrSR" "Parasite" 99
Rscript analysis/network/plotting/1metrics.R "YrSR" "Parasite"

Rscript analysis/network/1metrics.R "YrSR" "Plant" 99
Rscript analysis/network/plotting/1metrics.R "YrSR" "Plant" 
