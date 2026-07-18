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
             by=c("Site","Year")) %>%
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
             by=c("Site","Year")) %>%
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
             by=c("Site","Year")) %>%
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


# # ---- Combine precipitation and temperature summaries ----
# climate_combined <- sample_climate %>%
#   left_join(
#     precip %>% select(Site, Year, EndDate, APi),
#     by = c("Site", "Year", "EndDate")
#   ) %>%
#   left_join(
#     gdd %>% select(Site, Year, EndDate, GDD),
#     by = c("Site", "Year", "EndDate")
#   )

# climate_combined <- tmean_combined %>%
#   left_join(precip %>% select(Site, Year, EndDate),
#     by = c("Site", "Year", "EndDate"))
# 
# climate_combined <- climate_combined %>% 
#   left_join(combined_gdd_api,
#             by = c("Site", "Year", "EndDate"))

# ###########################################
# ## Calculate precipitation and temperature variability ##
# ###########################################
# 
# ## annual standard deviation by site and year ##
# annual_precip_var <- daily_precip_data %>% 
#         # marking May 1st as time of year where temperatures are above freezing
#   mutate(baseline_start = as.Date(paste0(Year, "-05-01")),
#          # marking mid oct to account for the 2021 going into oct
#          season_end = as.Date(paste0(Year, "-10-15"))) %>%
#   filter(Date >= baseline_start &  Date <= season_end) %>%
#   group_by(Site, Year) %>% 
#   summarise(PrecipSD = sd(Precip, na.rm = TRUE))
# 
# ## Calculate annual temp and temp anomaly variability ##
# annual_temp_var <- daily_temp_data %>% 
#           # marking May 1st as time of year where temperatures are above freezing
#   mutate(baseline_start = as.Date(paste0(Year, "-05-01")),
#          # marking mid oct to account for the 2021 going into oct
#          season_end = as.Date(paste0(Year, "-10-15"))) %>%
#   filter(Date >= baseline_start &  Date <= season_end) %>%
#   group_by(Site, Year) %>% 
#   summarise(TempDifSD = sd(tmax - tmin, na.rm = TRUE), 
#             TempMeanSD = sd(Tmean_normal, na.rm = TRUE),
#             TempAnomSD = sd(Tmean_anomaly, na.rm = TRUE), .groups = "drop")
# 
# ## Merge temp and precipitation variability ##
# climate_variability <- annual_precip_var %>%
#   left_join(annual_temp_var, by = c("Site", "Year")) %>%
#   select(Site, Year, PrecipSD, TempDifSD, TempMeanSD, TempAnomSD)
# 
# climate_combined <- climate_combined %>% 
#   left_join(climate_variability,
#             by = c("Site", "Year"))
# 
# climate_combined <- climate_combined %>%
#   select(
#     Site, Year, SampleRound, StartDate, EndDate,
#     APi, WindowTmean,
#     WindowTmeanAnom,
#     RoundTmeanAnom,
#     GDD,
#     PrecipSD, TempDifSD, TempMeanSD, TempAnomSD
#   )
# 
# write.csv(climate_combined, file =
#                               "data/relational/original/climate.csv", row.names=FALSE)



# ## ---- Compute Spring Precipitation (May 1 → day before first round) ----
# spring_precip <- sample_windows %>%
#   group_by(Site, Year) %>%
#   slice_min(SampleRound) %>%         # first round per Site/Year
#   ungroup() %>%
#   select(Site, Year, baseline_start, StartDate) %>%
#   rename(window_start = baseline_start, round_start = StartDate) %>%
#   left_join(daily_precip_data, by = c("Site", "Year")) %>%
#   filter(Date >= window_start & Date < round_start) %>%
#   group_by(Site, Year) %>%
#   summarize(SpringPrecip = sum(Precip, na.rm = TRUE), .groups = "drop")

## ---- Compute Round Precipitation (during each round) ----
## round_precip <- sample_windows %>%
##   left_join(daily_precip_data, by = c("Site", "Year"), relationship = "many-to-many") %>%
##   filter(Date >= StartDate & Date <= EndDate) %>%
##   group_by(Site, Year, SampleRound, StartDate, EndDate) %>%
##   summarize(RoundPrecip = sum(Precip, na.rm = TRUE), .groups = "drop")

# round_precip <- sample_windows %>%
#   rowwise() %>%
#   mutate(
#     RoundPrecip = sum(
#       daily_precip_data$Precip[
#         daily_precip_data$Site == Site &
#           daily_precip_data$Year == Year &
#           daily_precip_data$Date >= StartDate &
#           daily_precip_data$Date <= EndDate
#       ],
#       na.rm = TRUE
#     )
#   ) %>%
#   ungroup() %>%
#   select(Site, Year, SampleRound, StartDate, EndDate, RoundPrecip)

## ---- Compute Cumulative Precipitation (spring + all previous rounds) ----
# Step 1: Combine spring + round precipitation
# precip_combined <- round_precip %>%
#   left_join(spring_precip, by = c("Site", "Year")) %>%
#   arrange(Site, Year, SampleRound) %>%
#   group_by(Site, Year) %>%
#   mutate(
#     # Cumulative sum: spring precipitation + sum of previous rounds
#     CumulativePrecip = SpringPrecip + cumsum(lag(RoundPrecip, default = 0))
#   ) %>%
#   ungroup() %>%
#   select(Site, Year, SampleRound, StartDate, EndDate,
#          SpringPrecip, RoundPrecip, CumulativePrecip)

# ## ---- Compute Spring Tmean (May 1 → day before first round) ----
# spring_tmean <- sample_windows %>%
#   group_by(Site, Year) %>%
#   slice_min(SampleRound) %>%
#   ungroup() %>%
#   select(Site, Year, baseline_start, StartDate) %>%
#   rename(window_start = baseline_start, round_start = StartDate) %>%
#   left_join(daily_temp_data, by = c("Site", "Year")) %>%
#   filter(Date >= window_start & Date < round_start) %>%
#   group_by(Site, Year) %>%
#   summarize(
#     SpringTmean = mean(tmean, na.rm = TRUE),
#     SpringTmeanAnom = mean(Tmean_anomaly, na.rm = TRUE),
#     .groups = "drop"
#   )
#  
## ---- Compute Round Tmean ----
# round_tmean <- sample_windows %>%
#   left_join(daily_temp_data, by = c("Site", "Year"), relationship = "many-to-many") %>%
#   filter(Date >= StartDate & Date <= EndDate) %>%
#   group_by(Site, Year, SampleRound, StartDate, EndDate) %>%
#   summarize(
#     RoundTmean = mean(tmean, na.rm = TRUE),
#     RoundTmeanAnom = mean(Tmean_anomaly, na.rm = TRUE),
#     .groups = "drop"
#   )

# round_tmean <- sample_windows %>%
#   rowwise() %>%
#   mutate(
#     RoundTmean = mean(
#       daily_temp_data$tmean[
#         daily_temp_data$Site == Site &
#           daily_temp_data$Year == Year &
#           daily_temp_data$Date >= StartDate &
#           daily_temp_data$Date <= EndDate
#       ],
#       na.rm = TRUE
#     ),
#     RoundTmeanAnom = mean(
#       daily_temp_data$Tmean_anomaly[
#         daily_temp_data$Site == Site &
#           daily_temp_data$Year == Year &
#           daily_temp_data$Date >= StartDate &
#           daily_temp_data$Date <= EndDate
#       ],
#       na.rm = TRUE
#     )
#   ) %>%
#   ungroup() %>%
#   select(
#     Site, Year, SampleRound, StartDate, EndDate,
#     RoundTmean, RoundTmeanAnom
#   )
# 
# ## ---- Compute Cumulative Tmean (spring + previous rounds) ----
# tmean_combined <- round_tmean %>%
#   left_join(spring_tmean, by = c("Site", "Year")) %>%
#   arrange(Site, Year, SampleRound) %>%
#   group_by(Site, Year) %>%
#   mutate(
#     CumulativeTmean = SpringTmean + cumsum(lag(RoundTmean, default = 0)),
#     CumulativeTmeanAnom = SpringTmeanAnom + cumsum(lag(RoundTmeanAnom, default = 0))
#   ) %>%
#   ungroup() %>%
#   select(Site, Year, SampleRound, StartDate, EndDate,
#          SpringTmean, RoundTmean, CumulativeTmean,
#          SpringTmeanAnom, RoundTmeanAnom, CumulativeTmeanAnom)

# # ===============================================
# # Compute cumulative and round precipitation
# # using daily PRISM data and sample rounds
# # ===============================================
# 
# rm(list = ls())
# library(tidyverse)
# library(lubridate)
# library(terra)
# library(sf)
# # ---- Paths ----
# setwd("~/")
# source("lab_paths.R")
# local.path
# dir.bombus <- file.path(local.path, "skyIslands")
# 
# sites_shp_path <- file.path(local.path, "skyIslands_saved", "spatial", "sites.shp")
# daily_precip_dir <- "data/PRISM_data/PRISM_daily_precip"
# weather_csv <- "data/relational/original/weather.csv"
# 
# # Precipitation and temperature data used in this analysis were downloaded
# # from the PRISM Climate Group’s Explorer portal:
# # https://prism.oregonstate.edu/explorer/
# #
# # - Used file "skyIslands_saved/spatial/prism_site_locations.csv"
# # - Date of download: January 08, 2026
# # - Variables: Precipitation
# # - Temporal resolution: Daily
# # - Spatial resolution: 800 meters
# # - File format: BIL (.bil)
# # - Dataset used: PRISM Climate Group, Oregon State University
# 
# # ---- Load site shapefile ----
# sites_shp <- vect(sites_shp_path)
# sites_shp <- project(sites_shp, crs(sites_shp))  # ensure CRS is consistent
# print(sites_shp)
# 
# # ---- Load weather / sample round dates ----
# setwd(file.path(local.path, "skyIslands_saved"))
# weather <- read_csv(weather_csv) %>%
#   mutate(
#     StartDate = mdy(StartDate),
#     Year = as.numeric(Year)
#   ) %>%
#   # Remove SampleRound == 0
#   filter(SampleRound != 0)
# 
# # Summarize sample rounds
# sample_windows <- weather %>%
#   group_by(Site, Year, SampleRound) %>%
#   summarize(
#     StartDate = min(StartDate, na.rm = TRUE),
#     EndDate   = max(Date, na.rm = TRUE),
#     .groups = "drop"
#   )  %>%
#   mutate(baseline_start = as.Date(paste0(Year, "-05-01"))) # Add baseline start (May 1)
# 
# # ---- Load all Daily Precipitation Data ----
# daily_files <- list.files(daily_precip_dir, pattern = "\\.csv$", full.names = TRUE)
# 
# daily_precip_list <- list()
# 
# for (f in daily_files) {
#   year_file <- str_extract(basename(f), "\\d{4}") %>% as.numeric()
#   
#   # Read CSV, skip first 10 metadata lines
#   daily_data <- read_csv(f, skip = 10) %>%
#     # Keep only columns we need: Name and ppt
#     select(Name, `ppt (mm)`, Date) %>%
#     rename(Site = Name, Precip = `ppt (mm)`) %>%
#     mutate(
#       Date = ymd(Date),
#       Year = year(Date)
#     )
#   
#   daily_precip_list[[as.character(year_file)]] <- daily_data
# }
# 
# daily_precip_data <- bind_rows(daily_precip_list) # all daily PRISM precipitation values for the sites
# 
# # ---- Tie Daily Precip to Sample Rounds ----
# # Sort and define windows for cumulative precipitation
# sample_windows <- sample_windows %>%
#   arrange(Site, Year, SampleRound) %>%
#   group_by(Site, Year) %>%
#   mutate(
#     prev_end = lag(EndDate, default = first(baseline_start)),  # baseline_start of first round in the group
#     window_start = prev_end,
#     window_end   = StartDate - days(1)  # cumulative up to day before the round starts
#   ) %>%
#   ungroup()
# 
# # Sum precipitation for each window 
# # Initialize list to store cumulative results
# cumulative_precip_list <- list()
# 
# for (i in 1:nrow(sample_windows)) {
#   rw <- sample_windows[i, ]
#   
#   site_precip <- daily_precip_data %>%
#     filter(
#       Site == rw$Site,
#       Date >= rw$window_start,
#       Date <= rw$window_end
#     ) %>%
#     summarize(CumulativePrecip = sum(Precip, na.rm = TRUE))
#   
#   cumulative_precip_list[[i]] <- tibble(
#     Site        = rw$Site,
#     Year        = rw$Year,
#     SampleRound = rw$SampleRound,
#     window_start = rw$window_start,
#     window_end   = rw$window_end,
#     CumulativePrecip = site_precip$CumulativePrecip
#   )
# }
# 
# cumulative_precip <- bind_rows(cumulative_precip_list) # vectorized table of cumulative precip leading up to each survey round
# 
# save(cumulative_precip, file = "data/PRISM_data/cumulative_precip_per_round.Rdata")
# 
# # Cumulative precipitation during the round
# round_precip <- sample_windows %>% # cumulative precipitation during the actual sample round
#   rowwise() %>%
#   mutate(
#     RoundPrecip = sum(
#       daily_precip_data$Precip[
#         daily_precip_data$Site == Site &
#           daily_precip_data$Date >= StartDate &
#           daily_precip_data$Date <= EndDate
#       ],
#       na.rm = TRUE
#     )
#   )




############################################################################################################################




# # Precipitation and temperature data used in this analysis were downloaded
# # from the PRISM Climate Group’s Explorer portal:
# # https://prism.oregonstate.edu/explorer/
# #
# # - Used file "skyIslands_saved/spatial/prism_site_locations.csv"
# # - Date of download: August 17, 2025
# # - Variables: Precipitation and Temperature
# # - Temporal resolution: Monthly
# # - Spatial resolution: 800 meters
# # - File format: BIL (.bil)
# # - Dataset used: PRISM Climate Group, Oregon State University
# 
# 
# # ---- Load Site Shapefile ----
# setwd(file.path(local.path, "skyIslands_saved"))
# 
# sites_shp <- vect("spatial/sites.shp")
# print(crs(sites_shp))
# 
# # ---- Load weather / sample round dates ----
# weather <- read_csv("data/relational/original/weather.csv")
# 
# # Convert to Date objects
# weather <- weather %>%
#   mutate(
#     StartDate = mdy(StartDate)
#   )
# 
# # Summarize to define SampleRound windows per SiteSubSite
# sample_windows <- weather %>%
#   group_by(Site, Year, SampleRound) %>%
#   summarize(
#     StartDate = min(StartDate, na.rm = TRUE),
#     EndDate   = max(Date, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # ---- Define path to PRISM data directory ----
# prism_dir <- "data/PRISM_data"
# 
# # ---- Import Precipitation Data ----
# ## ---- Summer Monsoon Precipitation (July, Aug, Sep) ----
# years <- c("2011", "2016", "2017", "2020", "2021")
# months <- c("07", "08", "09")
# base_path_summer <- file.path(prism_dir, "summer_precip")
# 
# monthly_results <- list()
# 
# for (year in years) {
#   for (month in months) {
#     
#     raster <- rast(
#       file.path(base_path_summer, year,
#                 paste0("PRISM_ppt_stable_4kmM3_", year, month, "_bil.bil"))
#     )
#     
#     precip_values <- extract(raster, sites_shp)
#     
#     site_df <- as.data.frame(sites_shp)
#     site_df$Year  <- as.numeric(year)
#     site_df$Month <- as.numeric(month)
#     site_df$Precip <- precip_values[, 2]
#     
#     monthly_results[[paste(year, month)]] <-
#       site_df %>%
#       group_by(Site, Year, Month) %>%
#       summarize(Mean_Precip = mean(Precip, na.rm = TRUE), .groups = "drop")
#   }
# }
# 
# monthly_precip_data <- bind_rows(monthly_results)
# 
# save(monthly_precip_data,
#      file = file.path("data/PRISM_data/summer_precip/monthly_precip_site_data.Rdata"))
# 
# # ---- Use sample_windows to compute cumulative summer precipitation ----
# # Define summer start
# sample_windows <- sample_windows %>%
#   mutate(
#     summer_start = as.Date(paste0(Year, "-05-01"))
#   )
# 
# # Expand each round into contributing months
# round_months <- sample_windows %>%
#   mutate(
#     summer_start = as.Date(paste0(Year, "-05-01")),
#     end_month = floor_date(StartDate - days(1), "month")
#   ) %>%
#   filter(end_month >= summer_start) %>%   # <-- critical guard
#   rowwise() %>%
#   mutate(
#     month_seq = list(seq(
#       floor_date(summer_start, "month"),
#       end_month,
#       by = "1 month"
#     ))
#   ) %>%
#   unnest(month_seq) %>%
#   mutate(
#     Month = month(month_seq)
#   ) %>%
#   ungroup()
# 
# # Join monthly PRISM precipitation
# round_precip <- round_months %>%
#   left_join(monthly_precip_data,
#             by = c("Site", "Year", "Month")) %>%
#   group_by(Site, Year, SampleRound) %>%
#   summarize(
#     CumMay1Precip = sum(Mean_Precip, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # # ---- (Old) Import Precipitation Data ----
# # 
# # ## Summer Monsoon Precipitation (July, Aug, Sep) --------------------------
# # years <- c("2011", "2016", "2017", "2020", "2021")
# # months <- c("07", "08", "09")
# # base_path_summer <- file.path(prism_dir, "summer_precip")
# # 
# # summer_rasters <- list()
# # 
# # for (year in years) {
# #   monthly_rasters <- list()
# #   for (month in months) {
# #     file_path <- file.path(base_path_summer, year, paste0("PRISM_ppt_stable_4kmM3_", year, month, "_bil.bil"))
# #     monthly_rasters[[month]] <- rast(file_path)
# #   }
# #   stacked_rasters <- rast(monthly_rasters)
# #   summer_rasters[[year]] <- sum(stacked_rasters, na.rm = TRUE)
# #   print(paste("Processed summer precipitation for year:", year))
# # }
# # 
# # ## ---- Winter Precipitation (Dec previous year + Jan current year) ----------
# # winter_years <- c("2012", "2017", "2018", "2021", "2022")
# # base_path_winter <- file.path(prism_dir, "winter_precip")
# # 
# # winter_rasters <- list()
# # 
# # for (year in winter_years) {
# #   dec_year <- as.character(as.numeric(year) - 1)
# #   dec_path <- file.path(base_path_winter, dec_year, paste0("PRISM_ppt_stable_4kmM3_", dec_year, "12_bil.bil"))
# #   jan_path <- file.path(base_path_winter, year, paste0("PRISM_ppt_stable_4kmM3_", year, "01_bil.bil"))
# # 
# #   dec_raster <- rast(dec_path)
# #   jan_raster <- rast(jan_path)
# # 
# #   winter_rasters[[year]] <- sum(c(dec_raster, jan_raster), na.rm = TRUE)
# #   print(paste("Processed winter precipitation for year:", year))
# # }
# # 
# # # ---- Reproject sites to match raster CRS ------------------------------
# # sites_shp <- project(sites_shp, crs(summer_rasters[[1]]))
# # 
# # # ---- Extract Precipitation Data -------------
# # 
# # ## ---- Extract Summer Monsoon Precipitation values at sites -----------------
# # monsoon_results <- list()
# # 
# # for (year in years) {
# #   monsoon_raster <- summer_rasters[[year]]
# #   precip_values <- extract(monsoon_raster, sites_shp)
# # 
# #   site_df <- as.data.frame(sites_shp)
# #   site_df$Monsoon_Precipitation <- precip_values[, 2]  # Extracted values
# #   site_df$Year <- as.numeric(year)
# # 
# #   site_summary <- site_df %>%
# #     group_by(Site, Year) %>%
# #     summarize(Mean_Monsoon_Precip = mean(Monsoon_Precipitation, na.rm = TRUE), .groups = "drop")
# # 
# #   monsoon_results[[year]] <- site_summary
# # }
# # 
# # monsoon_precip_data <- bind_rows(monsoon_results)
# # 
# # ## ---- Extract Winter Precipitation values at sites --------------------------
# # winter_results <- list()
# # 
# # for (year in winter_years) {
# #   winter_raster <- winter_rasters[[year]]
# #   precip_values <- extract(winter_raster, sites_shp)
# # 
# #   site_df <- as.data.frame(sites_shp)
# #   site_df$Winter_Precipitation <- precip_values[, 2]
# #   site_df$Winter_Year <- as.numeric(year)
# # 
# #   site_summary <- site_df %>%
# #     group_by(Site, Winter_Year) %>%
# #     summarize(Mean_Winter_Precip = mean(Winter_Precipitation, na.rm = TRUE), .groups = "drop")
# # 
# #   winter_results[[year]] <- site_summary
# # }
# # 
# # winter_precip_data <- bind_rows(winter_results)
# # 
# # # ---- Final Prep ----------
# # 
# # # Rename Winter_Year to Year for consistency
# # winter_precip_data <- winter_precip_data %>%
# #   rename(Year = Winter_Year)
# # 
# # # Shift monsoon data by one year (since it's antecedent year monsoon precip)
# # monsoon_precip_data <- monsoon_precip_data %>%
# #   mutate(Year = Year + 1)
# # 
# # # ---- Save processed precipitation data for downstream analyses -------------
# # save(monsoon_precip_data, winter_precip_data, file = file.path('data/PRISM_data/precipitation_site_data.Rdata'))
# 
# # # ---- Join precipitation data to spec_net --------
# # 
# # ## ---- Load specimen table ----
# # setwd(dir.bombus)
# # load("data/spec_net.Rdata")
# # 
# # ## ---- Load the processed climate data ----
# # setwd(file.path(local.path, "skyIslands_saved"))
# # load("data/PRISM_data/precipitation_site_data.Rdata")
# # 
# # ## ---- Join climate datasets to spec.net ---------
# # spec.net <- spec.net %>%
# #   left_join(monsoon_precip_data, by = c("Site", "Year")) %>%
# #   left_join(winter_precip_data,  by = c("Site", "Year"))
