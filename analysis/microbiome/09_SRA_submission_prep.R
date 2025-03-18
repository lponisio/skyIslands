## Preparing SRA file submission

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path

setwd(local.path)
setwd("skyIslands_saved")

## mandatory column names from SRA
## need to pull sample name, collection date (YYYY-mm-dd), USA: Arizona or NM, lat and long (38.98 N 77.11 W)

load("../skyIslands/data/spec_microbes.Rdata")

spec <-  spec.net %>%
  select(UniqueID, Date, Long, State) %>%
  arrange(UniqueID)

spec_untrans <- read.csv("../skyIslands/data/spec_RBCL_16s.csv") %>%
  filter(UniqueID %in% spec$UniqueID) %>%
  select(UniqueID, Lat) %>%
  arrange(UniqueID)

spec$Lat <- spec_untrans$Lat

View(spec)

## fixing lat long to adhere to SRA guidelines
# Function to convert lat/lon to NCBI SRA format
convert_lat_long <- function(lat, lon) {
  lat_dir <- ifelse(lat >= 0, "N", "S")
  lon_dir <- ifelse(lon >= 0, "E", "W")
  
  sprintf("%.5f %s %.5f %s", abs(lat), lat_dir, abs(lon), lon_dir)
}

# test
convert_lat_long(spec$Lat[1], spec$Long[1])

# Apply the function to the data frame
spec$lat_long <- mapply(convert_lat_long, spec$Lat, spec$Long)

# fixing locality to match SRA guidelines
## should be either "USA: Arizona" or "USA: New Mexico"
spec$Geography <- ifelse(spec$State == "NM", "USA: New Mexico", "USA: Arizona")


## now want to filter to just samples that we actually had microbe data for
microbes <- read.table("SI_pipeline/merged/16s/maps/combined-map-2018-2021-noCtrl.txt", header=TRUE)

# drop unnecessary cols and filter to just microbes in the sequencing map
spec_clean <- spec %>%
  select(UniqueID, Date, Geography, lat_long) %>%
  filter(UniqueID %in% microbes$id)

# save out and copy and paste into ncbi table widget
write.csv(spec_clean, "SI_pipeline/SRA_submission_table_prep.csv")
