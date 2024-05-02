rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)

setwd("skyIslands/analysis/microbiome/")

load("../../data/spec_RBCL_16s.Rdata")
library(dplyr)
library(gt)


spec.net.full <- spec.net

load("../../data/spec_microbes.Rdata")

these_screened_sites <- unique(spec.net$Site[spec.net$ScreenedMicrobes == 1])
these_screened_years <- c("2018", "2021")



spec.net.full$ScreenedMicrobes <- ifelse(spec.net.full$Site %in% these_screened_sites & spec.net.full %in% these_screened_years, "Yes", "No")
## making supplemental table of the number of surveys per year, meadow 


summary_data <- spec.net.full %>%
  dplyr::select(Site, Year, Meadow, SampleRound) %>% 
  group_by(Site, Meadow, Year) %>%
  mutate(SampleRound = as.factor(SampleRound)) %>%
  dplyr::summarise(Surveys = length(unique(SampleRound))) %>%
  distinct() %>%
  mutate(Sterile = ifelse(
    Site %in% these_screened_sites & Year %in% these_screened_years, "Yes", "No"
  )) %>% as.tibble()

gt_table <- summary_data %>%
  gt(rowname_col = "Year", groupname_col="Site")
gt_table  

gtsave(gt_table, "figures/tableS1_surveysummary.ltx")
gtsave(gt_table, "figures/tableS1_surveysummary.pdf")
