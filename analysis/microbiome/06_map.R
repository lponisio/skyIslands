rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("skyIslands_saved/spatial")

## map making for SI project

## packages
library(sf)
library(terra)
library(ggplot2)
devtools::install_github("16EAGLE/basemaps")
library(basemaps)

## load in data

## set projection
crdref <- "+proj=latlong +datum=NAD83"
## site points spatvector
site_points <- terra::vect("sites.shp", crs=crdref)



#plot(site_points)
basemap_rast <- basemap_raster(site_points, map_service = "esri", map_type = "natgeo_world_map")

site_points <- project(site_points, crs(basemap_rast))

## site points sf
site_points_sf <- sf::st_as_sf(site_points, crs=crs(basemap_rast))

basemap_gg <- basemap_ggplot(site_points, map_service = "esri", map_type = "natgeo_world_map")

basemap_gg + geom_sf(data=site_points_sf, aes(geometry=geometry)) + labs(x="",y="")
## base layer

# or as ggplot2 layer:
library(ggplot2)
ggplot() + 
  basemap_gglayer(site_points) +
  scale_fill_identity() + 
  geom_point(aes(site_points$Lat, site_points$Long))



