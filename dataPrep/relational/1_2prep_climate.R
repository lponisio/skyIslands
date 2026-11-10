## ================================================================
## Climate preprocessing
##
## Creates one climate record per Site × Year × SampleRound.
##
## Calculates:
##   - Cumulative precipitation (May 1 → sample round)
##   - Mean tmean, tmax, tmin
##   - Mean maximum VPD
##   - 30-year daily-normal climatology
##   - Climate anomalies
##   - Antecedent Precipitation Index (APi)
##   - Growing Degree Days (GDD)
##
## Output:
##   data/relational/original/climate.csv
## ================================================================
library(tidyverse)
library(lubridate)
source("dataPrep/relational/src/getAPi.R")
setwd("../skyIslands_saved/")

sites_shp_path <- file.path("spatial", "sites.shp")
weather_csv <- "data/relational/original/weather.csv"

# Precipitation paths
daily_precip_dir <- "data/PRISM_data/PRISM_daily_precip"

# Temperature paths
daily_temp_dir <- "data/PRISM_data/PRISM_daily_temp"
monthly_normals_file <- "data/PRISM_data/temperature/PRISM_tmin_tmax_stable_800m_200801_202212.csv"
daily_normals_file <- "data/PRISM_data/30-year_daily_normals/PRISM_ppt_tmin_tmean_tmax_vpdmax_30yr_normal_800m_daily_normals.csv"

# Max vapor pressure paths
daily_vpdmax_dir <- "data/PRISM_data/PRISM_daily_vpdmax"

## *******************************************************************
## Create a round level summary of weather (rounds extend over
## multiple days)
## *******************************************************************
weather <- read_csv(weather_csv) %>%
  mutate(
    Year = as.numeric(Year)
  )  %>%
  filter(Method == "Net")

# Summarize sample rounds
sample_windows <- weather %>%
  group_by(Site, Year, SampleRound) %>%
  summarize(
    StartDate = min(StartDate, na.rm = TRUE),
    EndDate   = max(Date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Site, Year, SampleRound) %>%
  group_by(Site, Year) %>%
  mutate(
    baseline_start = as.Date(paste0(Year, "-05-01")),
  ) %>%
  ungroup()

# Filter out sample window without matching network
sample_windows <- sample_windows %>%
  filter(
    !(Site == "RP" &
        Year == 2018 &
        SampleRound == 1)
  )

## *******************************************************************
# Calculate cumulative precipitation 
## *******************************************************************
daily_precip_files  <- list.files(daily_precip_dir, pattern = "\\.csv$",
                          full.names = TRUE)

daily_precip_data <- daily_precip_files  %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>%
            select(Name, `ppt (mm)`, Date) %>%
            rename(Site = Name, Precip = `ppt (mm)`) %>%
            mutate(Date = ymd(Date),
                   Year = year(Date)))

sample_precip <-
  daily_precip_data %>%
  inner_join(sample_windows,
             by=c("Site","Year"),
             relationship = "many-to-many")%>%
  filter(Date >= baseline_start,
         Date <= EndDate) %>%
  group_by(Site,Year,SampleRound) %>%
  summarise(
    CumPrecip=sum(Precip),
    .groups="drop")

## *******************************************************************
## Calculate temperature anomalies 
## *******************************************************************
daily_temp_files  <- list.files(daily_temp_dir, pattern = "\\.csv$", full.names = TRUE)


daily_temp_data <- daily_temp_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>% 
            rename(
              Site = Name,
              tmin = `tmin (degrees C)`,
              tmax = `tmax (degrees C)`,
              tmean = `tmean (degrees C)`,
              Date = Date
            ) %>%
            mutate(
              Date = ymd(Date),
              Year = year(Date),
              Month = month(Date)
            ))

sample_temp <-
  daily_temp_data %>%
  inner_join(sample_windows,
             by=c("Site","Year"),
             relationship = "many-to-many")%>%
  filter(Date >= baseline_start,
         Date <= EndDate) %>%
  group_by(Site,Year,SampleRound) %>%
  summarise(
    MeanTmean=mean(tmean),
    MeanTmax=mean(tmax),
    MeanTmin=mean(tmin),
    .groups="drop")

## *******************************************************************
## Read in vapor pressure data and summarize
## *******************************************************************
daily_vpdmax_files  <- list.files(daily_vpdmax_dir, pattern = "\\.csv$", full.names = TRUE)


daily_vpdmax_data <- daily_vpdmax_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>%
            select(Name, `vpdmax (hPa)`, Date) %>%
            rename(Site = Name, vpdmax = `vpdmax (hPa)`) %>%
            mutate(Date = ymd(Date),
                   Year = year(Date)))

sample_vpd <-
  daily_vpdmax_data %>%
  inner_join(sample_windows,
             by=c("Site","Year"),
             relationship = "many-to-many")%>%
  filter(Date >= baseline_start,
         Date <= EndDate) %>%
  group_by(Site,Year,SampleRound) %>%
  summarise(
    MeanVPD=mean(vpdmax),
    .groups="drop")

## *******************************************************************
## Join seasonal climate datasets
## *******************************************************************
sample_climate <-
  sample_temp %>%
  left_join(sample_precip,
            by = c("Site","Year","SampleRound")
  ) %>%
  left_join(sample_vpd,
            by = c("Site","Year","SampleRound"))

## *******************************************************************
## Calculate sample round anomalies
## *******************************************************************
daily_normals <- read_csv(daily_normals_file,
                          skip = 10)

daily_normals <- daily_normals %>%
  mutate(
    NormalDate = as.Date(
      paste0(Date, "-2000"),
      format = "%B-%d-%Y"))

sample_windows <- sample_windows %>%
  mutate(
    NormalStart = as.Date(
      paste0(
        "2000-",
        format(baseline_start, "%m-%d")
      )
    ),
    NormalEnd = as.Date(
      paste0(
        "2000-",
        format(EndDate, "%m-%d")
      )
    )
  )

climatology <-
  daily_normals %>%
  rename(Site = Name) %>%
  inner_join(
    sample_windows,
    by = "Site",
    relationship = "many-to-many"
  ) %>%
  filter(
    NormalDate >= NormalStart,
    NormalDate <= NormalEnd
  ) %>%
  group_by(
    Site,
    Year,
    SampleRound
  ) %>%
  summarise(
    MeanTmean = mean(`tmean (degrees C)`),
    MeanTmax  = mean(`tmax (degrees C)`),
    MeanTmin  = mean(`tmin (degrees C)`),
    CumPrecip = sum(`ppt (mm)`),
    MeanVPD   = mean(`vpdmax (hPa)`),
    .groups = "drop"
  )

sample_climate <-
  sample_climate %>%
  left_join(
    climatology,
    by = c("Site","Year","SampleRound"),
    suffix = c("", "_normal")
  )

sample_climate <-
  sample_climate %>%
  mutate(
    Tmean_anom = MeanTmean - MeanTmean_normal,
    Tmax_anom  = MeanTmax  - MeanTmax_normal,
    Tmin_anom  = MeanTmin  - MeanTmin_normal,
    Precip_anom = CumPrecip - CumPrecip_normal,
    VPD_anom   = MeanVPD   - MeanVPD_normal
  )

## *******************************************************************
## Calculate antecedent precipitation
## *******************************************************************
## This is a measure of cumulative precipitation to determine soil moisture
precip <- daily_precip_data %>%
  # Filter to start from May 1st of each year
  mutate(baseline_start = as.Date(paste0(Year, "-05-01")),
         # marking mid oct to account for the 2021 going into oct
         season_end = as.Date(paste0(Year, "-10-15"))) %>% 
  filter(Date >= baseline_start & Date <= season_end) %>% 
  #Calculate APi for every Site Year combo
  group_by(Site, Year) %>% 
  mutate(APi = getApi(Precip, finite = FALSE)) %>% 
  ungroup() %>% 
  rename(EndDate = Date) #change name for joining at the end

## *******************************************************************
## Calculate growing degree days 
## *******************************************************************
gdd <- daily_temp_data %>%
  
  # Filter to start from May 1st of each year
  mutate(baseline_start = as.Date(paste0(Year, "-05-01")),
         # marking mid oct to account for the 2021 going into oct
         season_end = as.Date(paste0(Year, "-10-15"))) %>%
  filter(Date >= baseline_start &  Date <= season_end) %>%
  
  # Calculate daily GDD (using base temperature of 5°C for higher elevation)
  mutate(
    base_temp = 5,  
    # Basic GDD calculation
    daily_gdd = pmax(0, tmean - base_temp)) %>%
  # Calculate cumulative GDD for each site and year
  group_by(Site, Year) %>%
  arrange(Date) %>%
  mutate(GDD = cumsum(daily_gdd)) %>%
  ungroup() %>% 
  rename(EndDate = Date)

## Combine the dataframes of gdd and api. Keep only the important columns
combined_gdd_api <-precip %>%
  select(Site, Year, EndDate, APi) %>%
  left_join(gdd,
            by = c("Site", "Year", "EndDate")) %>%
  select(Site, Year, EndDate, APi, GDD)

combined_gdd_api <- combined_gdd_api %>%
  left_join(
  sample_windows %>%
    select(Site, Year, SampleRound, EndDate),
  by = c("Site","Year","EndDate")
)

###########################################
## Combine GDD and APi with sample_climate
###########################################
climate_combined <-
  sample_climate %>%
  left_join(
    combined_gdd_api %>%
      select(Site, Year, SampleRound, APi, GDD),
    by = c("Site","Year","SampleRound"))

write.csv(climate_combined,"data/relational/original/climate.csv",row.names = FALSE)


