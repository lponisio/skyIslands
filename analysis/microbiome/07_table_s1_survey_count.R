rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)

setwd("skyIslands/analysis/microbiome/")

load("../../data/spec_RBCL_16s.Rdata")

## making supplemental table of the number of surveys per year, meadow 

table_s1 <- spec.net %>%
  select(Site, Year, SampleRound, Meadow) %>% 
  mutate(SampleRound = as.factor(SampleRound)) %>%
  group_by(Year, Site, Meadow) %>%
  summarize(SurveyNumber = length(unique(SampleRound)))
