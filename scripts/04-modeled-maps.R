library(here)
library(dplyr)
library(glmulti)
library(raster)
library(rgdal)
library(corrplot)

epg.metadata <- read.csv(here::here('data', 'epg_with_metadata_norm.csv'))

# determine species for each sample
epg.metadata$Species <- substr(epg.metadata$ID, 4, 7)

# bring in RAI at each data point (for those that fall within grid)
rai <- read.csv(here::here('data', 'LOD_RAI_and_EPG.csv'))
epg.metadata <- left_join(epg.metadata, rai)

# test all environmental covariates
fit <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                 term100 + term250 + panlg100 + panlg250, 
               data = epg.metadata, # for all species together
               method = "h", 
               crit="aic", # use AIC as criterion
               level=1) # only main effects are to be used, no interactions

# take a look at best model
summary(fit@objects[[1]])

# tree cover is important (fewer trees = more parasites)
# distance to lake is important (closer to lake = fewer parasites)
# fire frequency is important, but less so (more fire = fewer parasites)
# termite mound density is very slightly important (more termite mounds = fewer parasites)

# save to txt file
model.output <- capture.output(summary(fit@objects[[1]]))
cat("lod-model", model.output, file= here::here('results', 'lod-model.txt'), sep="\n", append=F)

# now do one species at a time

# aeme
fit.aeme <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                 term100 + term250 + panlg100 + panlg250, 
               data = subset(epg.metadata, Species == "AEME"), 
               method = "h", 
               crit="aic", # use AIC as criterion
               level=1) # only main effects are to be used, no interactions
model.output.aeme <- capture.output(summary(fit.aeme@objects[[1]]))
cat("Impala model", model.output.aeme, file= here::here('results', 'lod-model-by-species.txt'), sep="\n", append=F)
summary(fit.aeme@objects[[1]])
# termite mound density is the only predictor of impala parasite load (weak; slightly lower in areas with termite mounds)

# koel
fit.koel <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                      term100 + term250 + panlg100 + panlg250, 
                    data = subset(epg.metadata, Species == "KOEL"), 
                    method = "h", 
                    crit="aic", # use AIC as criterion
                    level=1) # only main effects are to be used, no interactions
model.output.koel <- capture.output(summary(fit.koel@objects[[1]]))
cat("Waterbuck model", model.output.koel, file= here::here('results', 'lod-model-by-species.txt'), sep="\n", append=T)
summary(fit.koel@objects[[1]])
# pan density (w/in 250 meters) is only predictor for waterbuck; parasite load is higher with more pans (very very weak effect)

# ouou
fit.ouou <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                      term100 + term250 + panlg100 + panlg250, 
                    data = subset(epg.metadata, Species == "OUOU"), 
                    method = "h", 
                    crit="aic", # use AIC as criterion
                    level=1) # only main effects are to be used, no interactions
model.output.ouou <- capture.output(summary(fit.ouou@objects[[1]]))
cat("Oribi model", model.output.ouou, file= here::here('results', 'lod-model-by-species.txt'), sep="\n", append=T)
summary(fit.ouou@objects[[1]])
# no important spatial predictors of oribi parasite load

# phaf
fit.phaf <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                      term100 + term250 + panlg100 + panlg250, 
                    data = subset(epg.metadata, Species == "PHAF"), 
                    method = "h", 
                    crit="aic", # use AIC as criterion
                    level=1) # only main effects are to be used, no interactions
model.output.phaf <- capture.output(summary(fit.phaf@objects[[1]]))
cat("Warthog model", model.output.phaf, file= here::here('results', 'lod-model-by-species.txt'), sep="\n", append=T)
summary(fit.phaf@objects[[1]])
# warthog parasite load is higher further from Lake Urema

# rear
fit.rear <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                      term100 + term250 + panlg100 + panlg250, 
                    data = subset(epg.metadata, Species == "REAR"), 
                    method = "h", 
                    crit="aic", # use AIC as criterion
                    level=1) # only main effects are to be used, no interactions
model.output.rear <- capture.output(summary(fit.rear@objects[[1]]))
cat("Reedbuck model", model.output.rear, file= here::here('results', 'lod-model-by-species.txt'), sep="\n", append=T)
summary(fit.rear@objects[[1]])
# all of the pan layers are associated with reedbuck but now that I think about it we really should only use one of these in each model!! 

# tran
fit.tran <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                      term100 + term250 + panlg100 + panlg250, 
                    data = subset(epg.metadata, Species == "TRAN"), 
                    method = "h", 
                    crit="aic", # use AIC as criterion
                    level=1) # only main effects are to be used, no interactions
model.output.tran <- capture.output(summary(fit.tran@objects[[1]]))
cat("Nyala model", model.output.tran, file= here::here('results', 'lod-model-by-species.txt'), sep="\n", append=T)
summary(fit.tran@objects[[1]])
# lots of covariates pop up here, but weak effects; nyala parasite load is higher close to the lake, close to rivers, in areas with less tree cover, and near pans

# trsc
fit.trsc <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                      term100 + term250 + panlg100 + panlg250, 
                    data = subset(epg.metadata, Species == "TRSC"), 
                    method = "h", 
                    crit="aic", # use AIC as criterion
                    level=1) # only main effects are to be used, no interactions
model.output.trsc <- capture.output(summary(fit.trsc@objects[[1]]))
cat("Bushbuck model", model.output.trsc, file= here::here('results', 'lod-model-by-species.txt'), sep="\n", append=T)
summary(fit.trsc@objects[[1]])
# only distance to river - bushbuck parasite load is higher close to rivers



# spatially map model results (only did this for the all-species-combined model)

# import raster files
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

# stack raster
# create raster stack at 2 meter x 2 meter resolution, masked to camera grid, and normalized (each layer with mean of 0 and SD of 1)
raster.stack.norm <- raster::stack(fire.norm, tree.hansen.norm, rivers.dist.norm, 
                                   road.dist.norm, urema.dist.norm, termites.100m.norm, 
                                   termites.250m.norm, pans.offflood.100m.norm, pans.offflood.250m.norm,
                                   panslarge.offflood.100m.norm, panslarge.offflood.250m.norm)

# change names to match
names(raster.stack.norm) <- c("fire", "tree", "river", "road", "urema", "term100", "term250", 
                              "pan100", "pan250", "panlg100", "panlg250")

# predict values throughout study area (using best model)
predict.best <- raster::predict(raster.stack.norm, fit@objects[[1]], type = "response")

# plot output
plot(predict.best)

# unfortunately, it ends up looking like a tree cover map

# bring in the sample data
all_locs <- readOGR(here::here('data', 'epg_with_metadata_norm.shp'))  %>% 
  spTransform(CRS("+proj=utm +south +zone=36 +ellps=WGS84"))

points(all_locs, add =T)
# hmm, I'm not sure why some of these samples fall outside of the raster extent!! it doesn't make sense, given how I did things, but ah well
# can fix this but it will take a few days to regenerate the rasters for the study area, not a priority

# save
pdf(here::here('figures', 'lod-model.pdf'))
plot(predict.best)
points(all_locs, add =T)
dev.off()

# what if we use animal activity as a predictor of EPG, rather than consider it an outcome?
# of course, some of the points fall outside of the camera grid; I did not yet account for that

# throw in animal activity also
fit2 <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                  term100 + term250 + panlg100 + panlg250 +
                  MajorSpecies_all_years, 
                data = epg.metadata,
                method = "h", 
                crit="aic", # use AIC as criterion
                level=1) # only main effects are to be used, no interactions

# take a look at best model
summary(fit2@objects[[1]])

# but of course, wildlife activity is also likely highly correlated to some of these factors - let's see
correlation.matrix <- cor(epg.metadata[,c(2:13, 88)],use="complete.obs")

# plot heatmap of correlations
pdf(here::here('figures', 'correlation-plot.pdf'))
corrplot(correlation.matrix, type="upper", method = "color", addCoef.col = "black", diag = F, tl.col="black", tl.srt=45) 
dev.off()

write.csv(correlation.matrix, here::here('results', 'metadata-correlation-matrix.csv'))

# we see here that animal activity is highly correlated to distance to Lake Urema (more animal activity near the lake)

# cannot spatially map output of this model, since we don't have a raster layer for animal RAI