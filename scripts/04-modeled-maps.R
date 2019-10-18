library(here)
library(dplyr)
library(glmulti)
library(raster)
library(corrplot)

epg.metadata <- read.csv(here::here('data', 'epg_with_metadata_norm.csv'))

# test all environmental covariates
fit <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                 term100 + term250 + panlg100 + panlg250, 
               data = epg.metadata, # change for each species
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


# spatially map model results

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

# save
pdf(here::here('figures', 'lod-model.pdf'))
plot(predict.best)
dev.off()

# what if we use animal activity as a predictor of EPG, rather than consider it an outcome?
# of course, some of the points fall outside of the camera grid; I did not yet account for that

# bring in overall RAI
rai <- read.csv(here::here('data', 'LOD_RAI_and_EPG.csv'))
epg.metadata <- left_join(epg.metadata, rai)

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