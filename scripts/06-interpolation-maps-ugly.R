library(gstat)

# create a bounding box polygon
study.area <- as(raster::extent(all_locs), "SpatialPolygons")
proj4string(study.area) <- "+proj=utm +south +zone=36 +ellps=WGS84"
plot(study.area)
plot(all_locs, add =T)

# convert to raster
# fire doesn't work well as a template because it has a hole at the river bend; switched to tree
#fire <- raster("GIS/Spatial data/Fire2/tif1.tif")
#plot(fire)
#study.area.raster <- intersect(fire, study.area)
#study.area.raster <- projectRaster(study.area.raster, crs = CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

# first bring in a raster representing tree cover for the entire park, to use as a template
tree <- raster("GIS/Spatial data/Tree cover/BinMosFilterFinal.tif")
plot(tree)

# crop tree to the study area
study.area.raster <- intersect(tree, study.area)
study.area.raster <- projectRaster(study.area.raster, crs = CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

# make coarser?
study.area.raster.coarse <- aggregate(study.area.raster, fact=100)


### was failing due to missing values in object - there was a single NA record in there (18TRSC_11) and I'm not sure why; should have been 0. Can investigate later but for now just replacing wiht 0
all_locs@data[is.na(all_locs@data)] <- 0

# inverse distance weighted interpolation
gs <- gstat(formula=EPG~1, locations=all_locs)
idw <- interpolate(study.area.raster.coarse, gs, na.rm=TRUE)
plot(idw)
plot(log(idw))

