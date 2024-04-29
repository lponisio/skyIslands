library(sf)
library(ggplot2)
library(raster)

library(gdalUtilities)
library(rgdal)

library(spatstat)
#library(maptools)
source("lab_paths.R")
#setwd("~/Dropbox/skyIslands_saved/data/spatial")

sp.dir <- local.path

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

sites <- SpatialPointsDataFrame(coords = cbind(geo$Long, geo$Lat),
                               data = geo,
                               proj4string = CRS(crs.std))

writeOGR(sites, file.path(sp.dir, "sites.shp"), "sites",
         driver="ESRI Shapefile")

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




