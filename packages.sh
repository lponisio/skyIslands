#!/usr/bin/env bash

## dataprep only
Rscript -e 'install.packages("fields", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("fossil", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("RSQLite", repos="http://cran.r-project.org")'

Rscript -e 'install.packages("igraph", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("lmerTest", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("lme4", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("vegan", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("bipartite", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("tidyverse", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("betalink", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("gridExtra", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("effects", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("viridis", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("SYNCSA", repos="http://cran.r-project.org")'

