library(sf)
library(ggplot2)

sp.dir <- "~/Documents/skyIslandsSpatial"

## nm.n <- st_read(file.path(sp.dir,
##                           "NMN/V5_Ecological_Response_Unit_NMN.shp"))
## nm.s <- st_read(file.path(sp.dir,
##                           "NMS/V5_Ecological_Response_Unit_NMS.shp"))
## az.n <- st_read(file.path(sp.dir,
##                           "AZN/V5_Ecological_Response_Unit_AZN.shp"))
## az.s <- st_read(file.path(sp.dir,
##                           "AZS/V5_Ecological_Response_Unit_AZS.shp"))

si.habitats <- c("Spruce-Fir Forest",
                 "Alpine and Tundra",
                 "Ponderosa Pine Forest")

## nm.n.si <- nm.n[nm.n$R3ERU %in% si.habitats,]

## ggplot() +
##   geom_sf(data = nm.n.si, size = 3) +
##   ggtitle("NM N") +
##   coord_sf()


veg <- sf::st_read(dsn =
                       file.path(sp.dir, "MidScaleVeg_LifeForm.gdb/MidScaleVeg_LifeForm.gdb"))
