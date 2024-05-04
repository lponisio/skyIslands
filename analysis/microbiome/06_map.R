rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("skyIslands_saved/spatial")

load("../../skyIslands/data/spec_microbes.Rdata")

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
library(ggspatial)
## load in data

crs.std <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

microbes_screened_sites <- c("CH", "HM","JC", "MM", "PL", "RP", "SC", "SM")

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
site_points$ScreenedMicrobes <- ifelse(site_points$Site %in% microbes_screened_sites, 1, 0)

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

## setting up plot labels
labels_to_plot <- site_points %>%
  group_by(Site) %>%
  select(Site, ScreenedMicrobes) %>%
  distinct() 

## dropping sites where we didn't sample
site_points <- site_points[!(site_points$Site %in% c("JM", "CC", "SH")),]


these_rows <- c()

for (i in site_points$Site){
  if (!(i %in% these_rows)){
    these_rows <- base::append(match(i, site_points$Site), these_rows)
  } else {}
  these_rows <- unique(these_rows)
}



map <-ggplot() +
  basemap_gglayer(bbox_new) + 
  geom_sf(data = subset(site_points[these_rows,], ScreenedMicrobes == "1"), 
          color = "black",
          fill = "green",
          pch=25,
          size=3,
          stroke=1.1) +
  geom_sf(data = subset(site_points[these_rows,], ScreenedMicrobes == "0"), 
          color = "black",
          fill = "orange",
          pch=21,
          size=2,
          stroke=1.1) +
  geom_sf_label_repel(data=subset(site_points[these_rows,]),
                      aes(label=Site, geometry=geometry),
                      #point.padding = 10,
                      min.segment.length = 0,
                      box.padding = 0.5,
                      nudge_x=1,
                      seed=21) +
  coord_sf(xlim = st_coordinates(bbox_new)[c(1,2),1], # min & max of x values
           ylim = st_coordinates(bbox_new)[c(2,3),2], expand = FALSE) +
  scale_fill_identity() + 
  xlab("Longitude") + ylab("Latitude") +
  annotation_north_arrow(location = "tl", style = north_arrow_fancy_orienteering)+
  annotation_scale(location = "br") +
  theme_minimal()

map
setwd("../../skyIslands/analysis/microbiome/figures/")
ggsave(map, file="map.pdf", height=6, width=4)



