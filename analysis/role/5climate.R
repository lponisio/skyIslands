######################################################################
## -- This script calculates annual climate variability -- ##
#####################################################################
setwd("C:/")
source("lab_paths.R")
local.path
dir.bombus <- file.path(local.path, "skyIslands")
setwd(dir.bombus)

setwd("../skyIslands_saved/")


sites_shp_path <- file.path("spatial", "sites.shp")

# Precipitation paths
daily_precip_dir <- "data/PRISM_data/PRISM_daily_precip"

# Temperature paths
daily_temp_dir <- "data/PRISM_data/PRISM_daily_temp"
monthly_normals_file <- "data/PRISM_data/temperature/PRISM_tmin_tmax_stable_800m_200801_202212.csv"

daily_files <- list.files(daily_precip_dir, pattern = "\\.csv$",
                          full.names = TRUE)

daily_precip_data <- daily_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>%
            select(Name, `ppt (mm)`, Date) %>%
            rename(Site = Name, Precip = `ppt (mm)`) %>%
            mutate(Date = ymd(Date),
                   Year = year(Date)))

daily_files <- list.files(daily_temp_dir, pattern = "\\.csv$", full.names = TRUE)


daily_temp_data <- daily_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>%   # adjust skip if needed
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

## -- annual standard deviation by site and year -- ##
annualPrecip.var <- daily_precip_data %>% 
  group_by(Site, Year) %>% 
  summarise(precipSD = sd(Precip, na.rm = TRUE))

## -- Load long-term monthly normals -- ##
monthly_normals <- read_csv(monthly_normals_file, skip = 10) %>%
  rename(
    Site = Name,
    Tmin = `tmin (degrees C)`,
    Tmax = `tmax (degrees C)`
  ) %>%
  mutate(
    Date = as.Date(paste0(Date, "-01")),
    Month = month(Date),
    Tmean = (Tmin + Tmax)/2
  ) %>%
  group_by(Site, Month) %>%
  summarize(Tmean_normal = mean(Tmean, na.rm = TRUE), .groups = "drop")

## -- Merge normals to daily data for anomalies -- ##
daily_temp_data <- daily_temp_data %>%
  left_join(monthly_normals, by = c("Site", "Month")) %>%
  mutate(Tmean_anomaly = tmean - Tmean_normal)

## -- Calculate annual temp and temp anomoly variability -- ##
annualTemp.var <- daily_temp_data %>% 
  group_by(Site, Year) %>% 
  summarise(tempSD = sd(Tmean_normal, na.rm = TRUE),
            tempAnomSD = sd(Tmean_anomaly, na.rm = TRUE), .groups = "drop")

## -- Merge temp and precip data -- ##
climateVar <- annualPrecip.var %>%
  left_join(annualTemp.var, by = c("Site", "Year")) %>%
  select(Site, Year, precipSD, tempSD, tempAnomSD)
getwd()

write.csv(climateVar, file = "saved/traits/climateVariability.csv")
