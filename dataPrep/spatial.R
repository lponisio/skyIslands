library(sf)
library(ggplot2)
library(terra)
library(ggspatial)
library(basemaps)
library(tidyverse)
library(gdalUtilities)
library(rgdal)

library(spatstat)

source("lab_paths.R")
setwd(file.path(local.path))

#sp.dir <- file.path(local.path, "skyIslands_saved/spatial")

geo <- read.csv("skyIslands_saved/data/relational/original/geography.csv")

crs.std <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"




## ***********************************************************************
## merge nm and az files and reproject
## ***********************************************************************
nm <- raster(file.path(sp.dir,
                       "nm/nm_landcover/nm_landcover.tif"))
nm2 <- projectRaster(nm, crs=crs.std)

az <- raster(file.path(sp.dir,
                       "az/az_landcover/az_landcover.tif"))
az2 <- projectRaster(az, crs=crs.std)

writeRaster(nm2, filename=file.path(sp.dir, "nm.tif"),
                  format="GTiff", overwrite=TRUE)

writeRaster(az2, filename=file.path(sp.dir, "az.tif"),
                  format="GTiff", overwrite=TRUE)

nm.az <- c("az.tif", "nm.tif")

e.nm <- bbox(nm2)
e.az <- bbox(az2)

e <- extent(min(c(e.az[1,1], e.nm[1,1])),
                max(c(e.az[1,2], e.nm[1,2])),
              min(c(e.az[2,1], e.nm[2,1])),
                  max(c(e.az[2,2], e.nm[2,2])))

template <- raster(e)
projection(template) <- crs.std
writeRaster(template, file="nm_az.tif", format="GTiff",
            overwrite=TRUE)

mosaic_rasters(gdalfile=nm.az,
               dst_dataset="nm_az.tif",of="GTiff")
gdalinfo("nm_az.tif")

nm.az.r <- raster(file.path(sp.dir,
                       "nm_az.tif"))

## ***********************************************************************
## sample sites
## ***********************************************************************
## Create dataframe with the points and save as a shapefile
#sites <- SpatialPointsDataFrame(coords = cbind(geo$Long, geo$Lat),
#                               data = geo,
#                              proj4string = CRS(crs.std))

#writeOGR(sites, file.path(sp.dir, "sites.shp"), "sites",
#         driver="ESRI Shapefile")

## We only one one point per mountain so filter to subsite 0. 
## That way Mt with two meadows only appear as one point. 
geo <- geo %>% 
  filter(!is.na(MtRange) & SubSite == 1 & MtRange != "Jemez" & Site != "UK") 
## Create shapefile
sites_sf<- st_as_sf(geo,
         coords = c("Long", "Lat"),
         crs = crs.std)
ggplot() +
  geom_sf(data = sites_sf) +
  ggtitle("Map of Plot Locations")

## Get basemap

## set defaults for the basemap
set_defaults(map_service = "esri", map_type = "usa_topo_maps")
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

map <- ggplot() +
  basemap_gglayer(bbox_new) + # Use new bbox to download the basemap
  geom_sf(data = site_points, color = "black") +
  coord_sf(xlim = st_coordinates(bbox_new)[c(1,2),1], # min & max of x values
           ylim = st_coordinates(bbox_new)[c(2,3),2], expand = FALSE) +
  geom_sf_text(data = site_points, aes(label = MtRange), size = 2.5,
    color = "black", nudge_y = 19500) +
  scale_fill_identity()+ 
  xlab("Longitude") + ylab("Latitude") +
  annotation_north_arrow(location = "tl", style = north_arrow_fancy_orienteering)+
  annotation_scale(location = "br")+
  theme_minimal()

ggsave("skyIslands_sites.jpg", path = "skyIslands_saved/spatial", height=4, width=5)

## ***********************************************************************
## all the sw
## ***********************************************************************

sw <- raster(file.path(sp.dir, "swgap_landcover/swgap_landcover.tif"))
sw2 <- projectRaster(sw, crs=crs.std)

sw3 <- crop(sw2, y=e)

writeRaster(sw3, filename=file.path(sp.dir, "sw.tif"),
            format="GTiff", overwrite=TRUE)

## ***********************************************************************
## meadows
## ***********************************************************************

## meadows <- calc(nm.az.r, fun=function(x){x[x %in% c(70, 71, 86, 91)] <- NA;
##     return(x)} )

meadows <- raster(file.path(sp.dir,
                            "meadows.tif"))
meadows2 <- projectRaster(meadows, crs=crs.std)

meadows2 <- projectRaster(meadows2, crs=CRS("+proj=utm +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m
 +no_defs"))

sites <- spTransform(sites,
                     CRS("+proj=utm +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

## https://nowosad.github.io/post/lsm-bp1/
meadow.area <- sample_lsm(meadows2,
           y = sites,
           size =5000,
           what = "lsm_c_ca",
           shape = "circle")

meadow.area <- as.data.frame(meadow.area)

meadow.area$Site <- sites@data$Site

spplot(nm.az.r, scales = list(draw = TRUE),
       sp.layout = list("sp.points",
                        sites, pch=20, cex=2, col='black'))


spplot(meadows2, scales = list(draw = TRUE),
       sp.layout = list("sp.points",
                        sites, pch=1, cex=1, col='black'))




