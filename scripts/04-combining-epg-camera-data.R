### map fecal egg counts to camera locations and compare with camera records

library(here)
library(rgdal)
library(tidyverse)
library(broom)
library(sp)
library(ggplot2)
library(maptools)
library(plyr)
library(viridis)

# read in camera hexes
cams <- readOGR(here::here('data', 'spatial'), 'CameraGridHexes')  %>% 
  spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

# bring in egg SPDF
all_locs <- readOGR(here::here('data', 'epg_with_metadata_norm.shp'))  %>% 
  spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

# see the overlap of cameras and samples
cam_labs <- as.data.frame(do.call(rbind, lapply(cams@polygons, function(x) rbind(x@Polygons[[1]]@labpt))))
cam_labs$lab <- cams@data$StudySite
plot(cams)
points(all_locs, col = "red")
text(cam_labs[,1], cam_labs[,2], labels=as.character(cam_labs[,3]))

# let's pare things down to the samples that within hexs
all_locs <- all_locs[cams,]

# extract column (ID code with hex row and column), merge with ID and EPG
all_locs <- cbind(all_locs[,c(1,2)], over(all_locs, cams)[,6])
names(all_locs) <- c("ID", "EPG", "StudySite")

# let's see the distribution by camera and remove those with less than 3 samples
table(all_locs$StudySite)
less_than <- names(which(table(all_locs$StudySite) < 3))
all_locs <- all_locs[-which(as.character(all_locs$StudySite) %in% less_than),]

# remove all hexes that don't have samples in them
cams <- cams[all_locs,]

# plot everything again, looks good
plot(cams)
points(all_locs, col = "red")
text(cam_labs[,1], cam_labs[,2], labels=as.character(cam_labs[,3]))

# determine species for each sample
all_locs$Species <- substr(all_locs$ID, 4, 7)

# calculate overall mean, median EPG and sample count for each camera hex
cam_epg <- as.data.frame(all_locs) %>% group_by(StudySite) %>% dplyr::summarise(mean_EPG = mean(EPG), median_EPG = median(EPG), count = n()) %>% as.data.frame()

# calculate mean, median EPG and sample count for floodplain species only, for each camera hex
cam_epg_flood <- as.data.frame(all_locs) %>% 
  filter(Species == "KOEL" | Species == "AEME" | Species == "REAR" | Species == "OUOU" | Species == "PHAF") %>% 
  group_by(StudySite) %>% 
  dplyr::summarise(mean_EPG_flood = mean(EPG), median_EPG_flood = median(EPG), count_flood = n()) %>% 
  as.data.frame()

# calculate mean, median EPG and sample count for each species in each camera hex
cam_epg_species <- as.data.frame(all_locs) %>% group_by(StudySite, Species) %>% dplyr::summarise(mean_EPG = mean(EPG), median_EPG = median(EPG), count = n()) %>% as.data.frame()
cam_epg_species <- cam_epg_species[cam_epg_species$count > 2,] # remove all species-hexes with less than 3 samples
cam_epg_species_med <- spread(cam_epg_species[, c(1,2,4)], key = Species, value = median_EPG)
cam_epg_species_mean <- spread(cam_epg_species[, c(1,2,3)], key = Species, value = mean_EPG)
names(cam_epg_species_med) <- c("StudySite", "AEME_mean", "KOEL_mean", "OUOU_mean", "PHAF_mean", "REAR_mean", "TRAN_mean", "TRSC_mean", "TRST_mean")
names(cam_epg_species_mean) <- c("StudySite", "AEME_med", "KOEL_med", "OUOU_med", "PHAF_med", "REAR_med", "TRAN_med", "TRSC_med", "TRST_med")

cam_epg <- join_all(list(cam_epg, cam_epg_flood, cam_epg_species_med, cam_epg_species_mean), by = "StudySite", type = "left")

# now bring in camera trap data
rai <- read.csv(here::here("data", "RAI_LOD_by_year_wide.csv"))

#  merge individual data points with RAI for the camera
all_locs_rai <- left_join(as.data.frame(all_locs), rai)
write.csv(all_locs_rai, here::here("data", "LOD_RAI_and_EPG.csv"))

# merge camera data with EPG (will drop cameras for which there is no EPG data)
rai <- left_join(cam_epg, rai)


# make spatial
row.names(cams) <- row.names(rai)
cams <- SpatialPolygonsDataFrame(cams, rai)

# get SPDF into format that can be plotted in ggplot2
cams.df <- tidy(cams)
cams$polyID <- sapply(slot(cams, "polygons"), function(x) slot(x, "ID"))
cams.df <- merge(cams.df, cams, by.x = "id", by.y="polyID")

# plot count
ggplot(cams.df, aes(long,lat,group=group,fill=count))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Number of Samples')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))
ggsave(here::here('figures', 'hex-sample-count.pdf'))

# plot medians
ggplot(cams.df, aes(long,lat,group=group,fill=median_EPG))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15)) +
  ggtitle('Parasite load')
ggsave(here::here('figures', 'hex-median-epg.pdf'))

# plot mean
ggplot(cams.df, aes(long,lat,group=group,fill=mean_EPG))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Mean EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15)) +
  ggtitle('Parasite load')
ggsave(here::here('figures', 'hex-mean-epg.pdf'))


# plot medians - floodplain species only
ggplot(cams.df, aes(long,lat,group=group,fill=median_EPG_flood))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15)) +
  ggtitle('Parasite load (floodplain species)')
ggsave(here::here('figures', 'hex-median-epg-flood.pdf'))

# plot mean - floodplain species only
ggplot(cams.df, aes(long,lat,group=group,fill=mean_EPG_flood)) + 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Mean EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15)) +
  ggtitle('Parasite load (floodplain species)')
ggsave(here::here('figures', 'hex-mean-epg-flood.pdf'))

## Activity of all major species (those included here: AEME, KOEL, OUOU, PHAF, REAR, TRAN, TRSC)
ggplot(cams.df, aes(long,lat,group=group,fill=MajorSpecies_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Animal activity')
ggsave(here::here('figures', 'hex-rai-allspecies.pdf'))


## Activity of all floodplain species (those included here: AEME, KOEL, OUOU, PHAF, REAR)
ggplot(cams.df, aes(long,lat,group=group,fill=FloodSpecies_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Animal activity (floodplain species)')
ggsave(here::here('figures', 'hex-rai-floodspecies.pdf'))

# relationship between activity and parasites for "Major Species"
ggplot(cams.df, aes(MajorSpecies_all_years, median_EPG)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('All species') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-allspecies-vs-epg-median.pdf'))

ggplot(cams.df, aes(MajorSpecies_all_years, mean_EPG)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('All species') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-allspecies-vs-epg-mean.pdf'))

# relationship between activity and parasites for "Floodplain Species"
ggplot(cams.df, aes(FloodSpecies_all_years, median_EPG_flood)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Floodplain species') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-floodspecies-vs-epg-median.pdf'))

ggplot(cams.df, aes(FloodSpecies_all_years, mean_EPG_flood)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Floodplain species') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-floodspecies-vs-epg-mean.pdf'))

### Species by species

### AEME
# eggs
ggplot(cams.df, aes(long,lat,group=group,fill=AEME_med))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Impala parasites')
ggsave(here::here('figures', 'hex-aeme-epg-median.pdf'))

# activity
ggplot(cams.df, aes(long,lat,group=group,fill=Impala_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Impala activity')
ggsave(here::here('figures', 'hex-aeme-rai.pdf'))

# relationship between activity and parasites
ggplot(cams.df, aes(Impala_all_years, AEME_med)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Impala') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-vs-epg-med-aeme.pdf'))

### KOEL
# eggs
ggplot(cams.df, aes(long,lat,group=group,fill=KOEL_med))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Waterbuck parasites')
ggsave(here::here('figures', 'hex-koel-epg-median.pdf'))

# activity
ggplot(cams.df, aes(long,lat,group=group,fill=Waterbuck_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Waterbuck activity')
ggsave(here::here('figures', 'hex-koel-rai.pdf'))

# relationship between activity and parasites
ggplot(cams.df, aes(Waterbuck_all_years, KOEL_med)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Waterbuck') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-vs-epg-med-koel.pdf'))

### OUOU
# eggs
ggplot(cams.df, aes(long,lat,group=group,fill=OUOU_med))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Oribi parasites')
ggsave(here::here('figures', 'hex-ouou-epg-median.pdf'))

# activity
ggplot(cams.df, aes(long,lat,group=group,fill=Oribi_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Oribi activity')
ggsave(here::here('figures', 'hex-ouou-rai.pdf'))

# relationship between activity and parasites
ggplot(cams.df, aes(Oribi_all_years, OUOU_med)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Oribi') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-vs-epg-med-ouou.pdf'))

### PHAF
# eggs
ggplot(cams.df, aes(long,lat,group=group,fill=PHAF_med))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Warthog parasites')
ggsave(here::here('figures', 'hex-phaf-epg-median.pdf'))

# activity
ggplot(cams.df, aes(long,lat,group=group,fill=Warthog_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Warthog activity')
ggsave(here::here('figures', 'hex-phaf-rai.pdf'))

# relationship between activity and parasites
ggplot(cams.df, aes(Warthog_all_years, PHAF_med)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Warthog') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-vs-epg-med-phaf.pdf'))

### REAR
# eggs
ggplot(cams.df, aes(long,lat,group=group,fill=REAR_med))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Reedbuck parasites')
ggsave(here::here('figures', 'hex-rear-epg-median.pdf'))

# activity
ggplot(cams.df, aes(long,lat,group=group,fill=Reedbuck_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Reedbuck Activity')
ggsave(here::here('figures', 'hex-rear-rai.pdf'))

# relationship between activity and parasites
ggplot(cams.df, aes(Reedbuck_all_years, REAR_med)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Reedbuck') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-vs-epg-med-rear.pdf'))

### TRAN
# eggs
ggplot(cams.df, aes(long,lat,group=group,fill=TRAN_med))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Nyala parasites')
ggsave(here::here('figures', 'hex-tran-epg-median.pdf'))

# activity
ggplot(cams.df, aes(long,lat,group=group,fill=Nyala_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Nyala activity')
ggsave(here::here('figures', 'hex-tran-rai.pdf'))

# relationship between activity and parasites
ggplot(cams.df, aes(Nyala_all_years, TRAN_med)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Nyala') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-vs-epg-med-tran.pdf'))

### TRSC
# eggs
ggplot(cams.df, aes(long,lat,group=group,fill=TRSC_med))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='Median EPG')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Bushbuck parasites')
ggsave(here::here('figures', 'hex-trsc-epg-median.pdf'))

# activity
ggplot(cams.df, aes(long,lat,group=group,fill=Bushbuck_all_years))+ 
  geom_polygon()+
  geom_path(color="white")+
  coord_equal()+
  scale_fill_viridis(option='magma', name='RAI')+
  theme_void()+
  theme(
    legend.position=c(0.85,0.15))+
  ggtitle('Bushbuck activity')
ggsave(here::here('figures', 'hex-trsc-rai.pdf'))

# relationship between activity and parasites
ggplot(cams.df, aes(Bushbuck_all_years, TRSC_med)) +
  geom_point() +
  geom_smooth(method='lm') +
  ggtitle('Bushbuck') +
  xlab('Relative Activity Index')
ggsave(here::here('figures', 'rai-vs-epg-med-trsc.pdf'))
