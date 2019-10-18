library(here)
library(glmulti)

epg.metadata <- read.csv(here::here('data', 'epg_with_metadata_norm.csv'))

# bring in overall RAI
rai <- read.csv("Matt_LOD/LOD_RAI_and_EPG.csv")

metadata <- left_join(epg.metadata, rai)

# test all environmental covariates
fit <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                 term100 + term250 + panlg100 + panlg250, 
               data = epg.metadata, # change for each species
               method = "h", 
               crit="aic", # use AIC as criterion
               level=1) # only main effects are to be used, no interactions

# throw in animal activity also
fit2 <- glmulti(EPG ~ tree + fire + river + urema + pan100 + pan250 + 
                  term100 + term250 + panlg100 + panlg250 +
                  MajorSpecies_all_years, 
                data = epg.metadata, # change for each species
                method = "h", 
                crit="aic", # use AIC as criterion
                level=1) # only main effects are to be used, no interactions

# spatially map model results

## Extract raster values at fecal sample locations
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