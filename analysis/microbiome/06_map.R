rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("skyIslands_saved/spatial")

## map making for SI project

## packages
library(sf)
library(ggrepel)
# install.packages("devtools")
devtools::install_github("yutannihilation/ggsflabel")
library(ggsflabel)
library(terra)
library(ggplot2)
devtools::install_github("16EAGLE/basemaps")
library(basemaps)

## load in data

crs.std <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

## site points spatvector
sites_sf <- terra::vect("sites.shp", crs=crs.std) %>%
  st_as_sf()

ggplot() +
  geom_sf(data = sites_sf) +
  ggtitle("Map of Plot Locations")


## set defaults for the basemap
set_defaults(map_service = "mapbox",
             map_type = "satellite",
             map_token = "pk.eyJ1IjoicmhheWVzNyIsImEiOiJjbHZteHh5b3QwN3k5MnJucHFicWs0NHBuIn0.QJxQToKiZQesEtHIF-x9Zg")


## Change the CRS to match the basemaps
site_points <- st_transform(sites_sf, crs = st_crs(3857))

## Create a new bounding box to avoid points in the corners
bbox_new <- st_bbox(site_points) # current bounding box

xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values



bbox_new[1] <- bbox_new[1] - (0.1 * xrange) # xmin - left
bbox_new[3] <- bbox_new[3] + (0.1 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.1 * yrange) # ymin - bottom
bbox_new[4] <- bbox_new[4] + (0.1 * yrange) # ymax - top
bbox_new <- bbox_new %>%  # take the bounding box make it a spatial object
  st_as_sfc()

map <-ggplot() +
  basemap_gglayer(bbox_new) + # Use new bbox to download the basemap
  geom_sf(data = site_points, 
          color = "black",
          fill = "green",
          pch=25,
          size=3,
          stroke=1.1) +
  coord_sf(xlim = st_coordinates(bbox_new)[c(1,2),1], # min & max of x values
           ylim = st_coordinates(bbox_new)[c(2,3),2], expand = FALSE) +
  scale_fill_identity()+ 
  xlab("Longitude") + ylab("Latitude") +
  annotation_north_arrow(location = "tl", style = north_arrow_fancy_orienteering)+
  annotation_scale(location = "br") +
  theme_minimal()

map 



