# Here, I imported the raw spatial data layers for Gorongosa
# then cropped them to the study area for the fecal parasite samples
# then ran some manipulations (calculating distance to feature, values in focal area)
# exported these in original units in "raw" folder
# then normalized within the study area (set mean to 0, SD to 1) and exported these in "normalized" folder

setwd("~/Dropbox/projects/GORONGOSA2/Camera Trap Grid/R")

library(rgdal)
library(tidyverse)
library(broom)
library(sp)
library(ggplot2)
library(maptools)
library(plyr)
library(viridis)
library(raster)
library(glmulti)
library(RStoolbox)

# import raw data layers 

## roads
roadsmajor.lyr <- readOGR("~/Dropbox/projects/GORONGOSA2/Camera Trap Grid/R/GIS/Spatial data/GIS layers", "Roads_major") %>% spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84")) 
plot(roadsmajor.lyr)

## Lake Urema
urema.lyr <- readOGR("~/Dropbox/projects/GORONGOSA2/Camera Trap Grid/R/GIS/Spatial data/GIS layers", "lake_urema_latlong") %>% spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84")) 

## rivers
rivers.lyr <- readOGR("~/Dropbox/projects/GORONGOSA2/Camera Trap Grid/R/GIS/Spatial data/GIS layers", "gnp_main_rivers_latlong") %>% spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

## Conservative Pans (no channels and floodplain)
panscon.lyr <- readOGR("~/Dropbox/projects/GORONGOSA2/Camera Trap Grid/R/GIS/Spatial data/Pan_conservative", "PanOutline_130429_2") %>% spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

## Only large pans >1km2 (no channels and floodplain)
panslarge.lyr <- readOGR("~/Dropbox/projects/GORONGOSA2/Camera Trap Grid/R/GIS/Spatial data/Pan_conservative", "Pan_largerthan1km2") %>% spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

## only take pans outside of floodplain
# large
landscapes <- readOGR("~/Dropbox/projects/GORONGOSA2/Camera Trap Grid/R/GIS/Data from Marc", "gnp_landscapes_park&buffer_latlong") %>% 
  spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))
floodplain <- landscapes[landscapes@data$NAME == "Rift Valley Riverine & Floodplain",]
panslarge.offflood <- gDifference(panslarge.lyr, floodplain)
# all
pans.offflood <- gDifference(panscon.lyr, floodplain)

# park boundary
boundary.lyr <- readOGR("~/Dropbox/projects/GORONGOSA2/Camera Trap Grid/R/GIS/Spatial data/GIS layers", "gnp_boundary_west_straight_latlong") %>%  spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))


# We want to crop the raster layers to the study area (where the samples are located) so we don't have huge data files. Doing termites first, then use as template for others so the extent AND resolution matches

# manually specify extent based on exploration of fecal data

# termites
termites <- raster("GIS/Spatial data/Termites/JulAugEdgeClip_MosNull_161026Reclass.tif") %>%
  projectRaster(crs = CRS("+proj=utm +south +zone=36 +ellps=WGS84")) %>%
  raster::crop(extent(c(630700, 660000, 7897000, 7920000)))
writeRaster(termites, 'Matt_LOD/Spatial data/termites.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW")) # export raster as .tif

# fire
fire <- raster("GIS/Spatial data/Fire2/tif1.tif") %>%
  projectRaster(crs = CRS("+proj=utm +south +zone=36 +ellps=WGS84")) %>%
  raster::resample(termites, method = 'bilinear')
writeRaster(fire, 'Matt_LOD/Spatial data/fire.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW")) # export raster as .tif

# Hansen tree layer
tree.hansen <- raster("GIS/Spatial data/Tree cover/Hansen_10S_030E_treecover2010_v3.tif") %>%
  projectRaster(crs = CRS("+proj=utm +south +zone=36 +ellps=WGS84")) %>%
  raster::resample(termites, method = 'bilinear')
writeRaster(tree.hansen, 'landscape-of-disgust/data/spatial/raw/tree.hansen.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW")) # export raster as .tif

# rasterize pan layers

panscon.raster <- rasterize(panscon.lyr, termites, field = 1, background = 0) 
panslarge.raster <- rasterize(panslarge.lyr, termites, field = 1, background = 0) 
pans.offflood.raster <- rasterize(pans.offflood, termites, field = 1, background = 0) 
panslarge.offflood.raster <- rasterize(panslarge.offflood, termites, field = 1, background = 0) 

writeRaster(panscon.raster, 'landscape-of-disgust/data/spatial/raw/panscon.raster.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(panslarge.raster, 'landscape-of-disgust/data/spatial/raw/panslarge.raster.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(pans.offflood.raster, 'landscape-of-disgust/data/spatial/raw/pans.offflood.raster.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(panslarge.offflood.raster, 'landscape-of-disgust/data/spatial/raw/panslarge.offflood.raster.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

# calculate distance from vector features

# roads
roads.pts = as(roadsmajor.lyr, "SpatialPointsDataFrame")
road.dist <- distanceFromPoints(termites, roads.pts)
writeRaster(road.dist, 'landscape-of-disgust/data/spatial/raw/road.dist.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

# rivers
rivers.pts = as(rivers.lyr, "SpatialPointsDataFrame") 
rivers.dist <- distanceFromPoints(termites, rivers.pts)
writeRaster(rivers.dist, 'landscape-of-disgust/data/spatial/raw/rivers.dist.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

# lake
urema.pts = as(urema.lyr, "SpatialLinesDataFrame") %>% as("SpatialPointsDataFrame") # first have to convert to lines, then to points (can't go right from polygon to points)
urema.dist <- distanceFromPoints(termites, urema.pts)
writeRaster(urema.dist, 'landscape-of-disgust/data/spatial/raw/urema.dist.tif', format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))




# calculate focal rasters

## termite mounds - didn't run all of these just because it was taking ages

termites.100m <- focal(termites, w = matrix(1/1201, nc = 101, nr = 101))
writeRaster(termites.crop.100m, 'landscape-of-disgust/data/spatial/raw/termites.100m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

termites.250m <- focal(termites, w = matrix(1/251001, nc = 251, nr = 251))
writeRaster(termites.250m, 'landscape-of-disgust/data/spatial/raw/termites.250m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

#termites.500m <- focal(termites, w = matrix(1/63001, nc = 501, nr = 501))
#writeRaster(termites.crop.500m, 'landscape-of-disgust/data/spatial/raw/termites.500m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

#termites.1km <- focal(termites, w = matrix(1/1002001, nc = 1001, nr = 1001))
#writeRaster(termites.crop.1km, 'landscape-of-disgust/data/spatial/raw/termites.1000m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

## pans - didn't run all of these just because it was taking ages

#panscon.100m <- focal(panscon.raster, w = matrix(1/1201, nc = 101, nr = 101))
#writeRaster(panscon.100m, 'landscape-of-disgust/data/spatial/raw/panscon.100m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#panscon.250m <- focal(panscon.raster, w = matrix(1/251001, nc = 251, nr = 251))
#writeRaster(panscon.250m, 'landscape-of-disgust/data/spatial/raw/panscon.250m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#panscon.500m <- focal(panscon.raster, w = matrix(1/63001, nc = 501, nr = 501))
#writeRaster(panscon.500m, 'landscape-of-disgust/data/spatial/raw/panscon.500m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#panscon.1km <- focal(panscon.raster, w = matrix(1/1002001, nc = 1001, nr = 1001))
#writeRaster(panscon.1km, 'landscape-of-disgust/data/spatial/raw/panscon.1km.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#panslarge.100m <- focal(panslarge.raster, w = matrix(1/1201, nc = 101, nr = 101))
#writeRaster(panslarge.100m, 'landscape-of-disgust/data/spatial/raw/panslarge.100m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#panslarge.250m <- focal(panslarge.raster, w = matrix(1/251001, nc = 251, nr = 251))
#writeRaster(panslarge.250m, 'landscape-of-disgust/data/spatial/raw/panslarge.250m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#panslarge.500m <- focal(panslarge.raster, w = matrix(1/63001, nc = 501, nr = 501))
#writeRaster(panslarge.500m, 'landscape-of-disgust/data/spatial/raw/panslarge.500m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#panslarge.1km <- focal(panslarge.raster, w = matrix(1/1002001, nc = 1001, nr = 1001))
#writeRaster(panslarge.1km, 'landscape-of-disgust/data/spatial/raw/panslarge.1km.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

pans.offflood.100m <- focal(pans.offflood.raster, w = matrix(1/1201, nc = 101, nr = 101))
writeRaster(pans.offflood.100m, 'landscape-of-disgust/data/spatial/raw/pans.offflood.100m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

pans.offflood.250m <- focal(pans.offflood.raster, w = matrix(1/251001, nc = 251, nr = 251))
writeRaster(pans.offflood.250m, 'landscape-of-disgust/data/spatial/raw/pans.offflood.250m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

#pans.offflood.500m <- focal(pans.offflood.raster, w = matrix(1/63001, nc = 501, nr = 501))
#writeRaster(pans.offflood.500m, 'landscape-of-disgust/data/spatial/raw/pans.offflood.500m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#pans.offflood.1km <- focal(pans.offflood.raster, w = matrix(1/1002001, nc = 1001, nr = 1001))
#writeRaster(pans.offflood.1km, 'landscape-of-disgust/data/spatial/raw/pans.offflood.1km.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

panslarge.offflood.100m <- focal(panslarge.offflood.raster, w = matrix(1/1201, nc = 101, nr = 101))
writeRaster(panslarge.offflood.100m, 'landscape-of-disgust/data/spatial/raw/panslarge.offflood.100m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

panslarge.offflood.250m <- focal(panslarge.offflood.raster, w = matrix(1/251001, nc = 251, nr = 251))
writeRaster(panslarge.offflood.250m, 'landscape-of-disgust/data/spatial/raw/panslarge.offflood.250m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

#panslarge.offflood.500m <- focal(panslarge.offflood.raster, w = matrix(1/63001, nc = 501, nr = 501))
#writeRaster(panslarge.offflood.500m, 'landscape-of-disgust/data/spatial/raw/panslarge.offflood.500m.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
#
#panslarge.offflood.1km <- focal(panslarge.offflood.raster, w = matrix(1/1002001, nc = 1001, nr = 1001))
#writeRaster(panslarge.offflood.1km, 'landscape-of-disgust/data/spatial/raw/panslarge.offflood.1km.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))


# normalize within study area (rectangular)

# first need to bring in all of the rasters again
fire <- raster("landscape-of-disgust/data/spatial/raw/fire.tif")
tree.hansen <- raster("landscape-of-disgust/data/spatial/raw/tree.hansen.tif")
rivers.dist <- raster("landscape-of-disgust/data/spatial/raw/rivers.dist.tif")
road.dist <- raster("landscape-of-disgust/data/spatial/raw/road.dist.tif")
urema.dist <- raster("landscape-of-disgust/data/spatial/raw/urema.dist.tif")
termites.100m <- raster("landscape-of-disgust/data/spatial/raw/termites.100m.tif")
termites.250m <- raster("landscape-of-disgust/data/spatial/raw/termites.250m.tif")
pans.offflood.100m <- raster("landscape-of-disgust/data/spatial/raw/pans.offflood.100m.tif")
pans.offflood.250m <- raster("landscape-of-disgust/data/spatial/raw/pans.offflood.250m.tif")
panslarge.offflood.100m <- raster("landscape-of-disgust/data/spatial/raw/panslarge.offflood.100m.tif")
panslarge.offflood.250m <- raster("landscape-of-disgust/data/spatial/raw/panslarge.offflood.250m.tif")

# normalize
fire.norm <- normImage(fire, norm = TRUE)
tree.hansen.norm <- normImage(tree.hansen, norm = TRUE)
rivers.dist.norm <- normImage(rivers.dist, norm = TRUE)
road.dist.norm <- normImage(road.dist, norm = TRUE)
urema.dist.norm <- normImage(urema.dist, norm = TRUE)
termites.100m.norm <- normImage(termites.100m, norm = TRUE)
termites.250m.norm <- normImage(termites.250m, norm = TRUE)
pans.offflood.100m.norm <- normImage(pans.offflood.100m, norm = TRUE)
pans.offflood.250m.norm <- normImage(pans.offflood.250m, norm = TRUE)
panslarge.offflood.100m.norm <- normImage(panslarge.offflood.100m, norm = TRUE)
panslarge.offflood.250m.norm <- normImage(panslarge.offflood.250m, norm = TRUE)

# export
writeRaster(fire.norm, 'landscape-of-disgust/data/spatial/normalized/fire.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(tree.hansen.norm, 'landscape-of-disgust/data/spatial/normalized/tree.hansen.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(rivers.dist.norm, 'landscape-of-disgust/data/spatial/normalized/rivers.dist.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(road.dist.norm, 'landscape-of-disgust/data/spatial/normalized/road.dist.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(urema.dist.norm, 'landscape-of-disgust/data/spatial/normalized/urema.dist.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(termites.100m.norm, 'landscape-of-disgust/data/spatial/normalized/termites.100m.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(termites.250m.norm, 'landscape-of-disgust/data/spatial/normalized/termites.250m.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(pans.offflood.100m.norm, 'landscape-of-disgust/data/spatial/normalized/pans.offflood.100m.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(pans.offflood.250m.norm, 'landscape-of-disgust/data/spatial/normalized/pans.offflood.250m.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(panslarge.offflood.100m.norm, 'landscape-of-disgust/data/spatial/normalized/panslarge.offflood.100m.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
writeRaster(panslarge.offflood.250m.norm, 'landscape-of-disgust/data/spatial/normalized/panslarge.offflood.250m.norm.tif', format = "GTiff", overwrite = TRUE, options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))