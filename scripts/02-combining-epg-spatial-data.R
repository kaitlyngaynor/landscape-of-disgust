library(here) # rather than set working directory; see https://www.tidyverse.org/articles/2017/12/workflow-vs-script/

source(here::here('scripts', 'r_workarounds.r')) # bring in R workarounds (with fuck_you_kml function)

library(rgdal)
library(tidyverse)
library(broom)
library(sp)
library(ggplot2)
library(maptools)
library(plyr)
library(viridis)
library(raster)

# bring in fecal sample location KML files

locs17 <- fuck_you_kml(here::here('data', 'para_locs_2017.kml'))
locs17$Latitude <- as.numeric(as.character(locs17$Latitude))
locs17$Longitude <- as.numeric(as.character(locs17$Longitude))
locs17$ID <- sapply(locs17$ID, function(x) paste('17_', x, sep=''))

locs18 <- fuck_you_kml(here::here('data', 'para_locs_2018.kml'))
locs18$Latitude <- as.numeric(as.character(locs18$Latitude))
locs18$Longitude <- as.numeric(as.character(locs18$Longitude))
locs18$ID <- sapply(locs18$ID, function(x) paste('18_', x, sep=''))

locs19 <- fuck_you_kml(here::here('data', 'para_locs_2019.kml'))
locs19$Latitude <- as.numeric(as.character(locs19$Latitude))
locs19$Longitude <- as.numeric(as.character(locs19$Longitude))
locs19$ID <- sapply(locs19$ID, function(x) paste('19_', x, sep=''))

all_locs <- rbind(locs17, locs18, locs19)

# bring in fecal egg count info
e17 <- read.csv(here::here('data', 'fecal_egg_counts_2017.csv'), header=TRUE)
e18 <- read.csv(here::here('data', 'fecal_egg_counts_2018.csv'), header=TRUE)
e19 <- read.csv(here::here('data', 'fecal_egg_counts_2019.csv'), header=TRUE)
e17$Year <- 2017
e18$Year <- 2018
e19$Year <- 2019

e18 <- e18[,c(1:8,10)] # need to get rid of extraneous 'X' notes column that is present in these two years
e19 <- e19[,c(1:8,10)]

# change names for 2019 data so they match (the columns ARE the same, they are just called different things)
names(e19) <- names(e18)

e17$SampleID <- sapply(e17$SampleID, function(x) paste('17_', x, sep=''))
e18$SampleID <- sapply(e18$SampleID, function(x) paste('18_', x, sep=''))
e19$SampleID <- sapply(e19$SampleID, function(x) paste('19_', x, sep=''))

para <- rbind(e17, e18, e19)

# some more data cleaning (Matt code)
uLers <- grep('uL', para$notes)
para$notes <- as.character(para$notes)
para$notes[uLers[1]] <- '23mL'
para$notes[uLers[2]] <- 'Collar #32011, 21mL'
para$notes[uLers[3]] <- '10mL'
para <- para[-which(is.na(para$total)),]
para <- para[-which(para$SampleID == '18_TRSC11'),]

# calculate egg count
mf <- function(x, mL=28) 1/((x/mL)*0.3)

parasites <- matrix(ncol=3, nrow=0)
samples <- unique(as.character(para$SampleID))
for(i in 1:length(samples)){
  dat <- para[grep(paste(samples[i],'$',sep=''), para$SampleID),]
  if(all(dat$notes == "")==FALSE){
    prob <- as.character(dat$notes[which(dat$notes != '')])
    found <- grep('mL', prob)
    if(length(found) == 0){
      epg <- mean(dat$total * mf(dat$fecal_weight), na.rm=TRUE)
      spp <- substr(dat$SampleID[1],4,7)
      parasites <- rbind(parasites, data.frame(Sample=dat$SampleID[1], Species=spp, EPG=epg, Date=dat$Date_Processed))
      
    }else{
      dat$vol <- NA
      dat$vol[found] <- as.numeric(substr(prob[found], 1,2))
      rows <- 1:nrow(dat)
      dat$vol[rows[-found]] <- 28
      mf_v <- function(x, mL) 1/((x/mL)*0.3)
      epg <- mean(dat$total * mf_v(dat$fecal_weight, dat$vol))
      spp <- substr(dat$SampleID[1],4,7)
      parasites <- rbind(parasites, data.frame(Sample=dat$SampleID[1], Species=spp, EPG=epg, Date=dat$Date_Processed))
    }
  }else{
    epg <- mean(dat$total * mf(dat$fecal_weight))
    spp <- substr(dat$SampleID[1],4,7)
    parasites <- rbind(parasites, data.frame(Sample=dat$SampleID[1], Species=spp, EPG=epg, Date=dat$Date_Processed))
  }
}

all_locs$EPG <- parasites$EPG[match(all_locs$ID, parasites$Sample)]

# convert samples into spatial points data frame
all_locs <- all_locs[, c('Longitude', 'Latitude')] %>% 
  SpatialPointsDataFrame(all_locs[, c('ID', 'EPG')], proj4string=CRS("+proj=longlat +ellps=WGS84")) %>%
  spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

# plot them
plot(all_locs)

# first bring in rasters
fire.norm <- raster(here::here("data", "spatial", "normalized", "fire.norm.tif"))
tree.hansen.norm <- raster(here::here("data", "spatial", "normalized", "tree.hansen.norm.tif"))
rivers.dist.norm <- raster(here::here("data", "spatial", "normalized", "rivers.dist.norm.tif"))
road.dist.norm <- raster(here::here("data", "spatial", "normalized", "road.dist.norm.tif"))
urema.dist.norm <- raster(here::here("data", "spatial", "normalized", "urema.dist.norm.tif"))
termites.100m.norm <- raster(here::here("data", "spatial", "normalized", "termites.100m.norm.tif"))
termites.250m.norm <- raster(here::here("data", "spatial", "normalized", "termites.250m.norm.tif"))
pans.offflood.100m.norm <- raster(here::here("data", "spatial", "normalized", "pans.offflood.100m.norm.tif"))
pans.offflood.250m.norm <- raster(here::here("data", "spatial", "normalized", "pans.offflood.250m.norm.tif"))
panslarge.offflood.100m.norm <- raster(here::here("data", "spatial", "normalized", "panslarge.offflood.100m.norm.tif"))
panslarge.offflood.250m.norm <- raster(here::here("data", "spatial", "normalized", "panslarge.offflood.250m.norm.tif"))

# create raster stack at 2 meter x 2 meter resolution, masked to camera grid, and normalized (each layer with mean of 0 and SD of 1)
raster.stack.norm <- raster::stack(fire.norm, tree.hansen.norm, rivers.dist.norm, 
                                   road.dist.norm, urema.dist.norm, termites.100m.norm, 
                                   termites.250m.norm, pans.offflood.100m.norm, pans.offflood.250m.norm,
                                   panslarge.offflood.100m.norm, panslarge.offflood.250m.norm)

# change names for ease later on (or ESRI driver will do it for you)
names(raster.stack.norm) <- c("fire", "tree", "river", "road", "urema", "term100", "term250", 
                              "pan100", "pan250", "panlg100", "panlg250")

# extract values at egg locations
all_locs <- raster::extract(raster.stack.norm, all_locs, sp = TRUE) # extract raster values
head(all_locs@data) # view metadata
write.csv(all_locs@data, file = here::here('data', 'epg_with_metadata_norm.csv'), row.names = FALSE)
          
# write as shape file 
writeOGR(all_locs, here::here('data', 'epg_with_metadata_norm.shp'), driver = 'ESRI Shapefile', layer = 'epg_with_metadata_norm')
