# Data-centric SDM Step 2) Import and process explanatory data
## Sara Hansen
## Modified July 25, 2024

library(tidyverse)
library(sf)
library(raster) # only for select functions
library(corrplot)
library(terra)
library(geodata)
library(exactextractr)
library(tidyterra)
library(fasterize)

################################# READ IN DATA #################################

#Import study areas created in Step 1
mich <- sf::read_sf("output/mich.shp")

sag <- sf::read_sf("output/sag.shp")

#Import the downloaded data and crop to each study area

# WorldClim bioclimatic variables 6, 10, 11, and 18
wclist <- c("output/TRUE/wc2.1_30s/wc2.1_30s_bio_6.tif",
            "output/TRUE/wc2.1_30s/wc2.1_30s_bio_10.tif",
            "output/TRUE/wc2.1_30s/wc2.1_30s_bio_11.tif",
            "output/TRUE/wc2.1_30s/wc2.1_30s_bio_18.tif")
wcM <- terra::rast(wclist) %>%
  terra::crop(mich, mask = T)

wcS <- terra::rast(wclist) %>%
  terra::crop(sag, mask = T)

# Solar radiation (mean across entire year)
solarlist <- list.files(path = "output/TRUE/wc2.1_30s", pattern = "srad",
                        full.names = TRUE) #should be 12
solarM <- terra::rast(solarlist) %>%
  terra::crop(mich, mask = T) %>%
  terra::app(fun = mean)
solarS <- terra::rast(solarlist) %>%
  terra::crop(sag, mask = T) %>%
  terra::app(fun = mean)

# Wind speed (mean across entire year)
wind_speedlist <- list.files(path = "output/TRUE/wc2.1_30s", pattern = "wind",
                             full.names = TRUE)
wind_speedM <- terra::rast(wind_speedlist) %>%
  terra::crop(mich, mask = T) %>%
  terra::app(fun = mean)
wind_speedS <- terra::rast(wind_speedlist) %>%
  terra::crop(sag, mask = T) %>%
  terra::app(fun = mean)

# Water vapor pressure (mean across entire year)
vaporlist <- list.files(path = "output/TRUE/wc2.1_30s", pattern = "vapr",
                        full.names = TRUE)
vaporM <- terra::rast(vaporlist) %>%
  terra::crop(mich, mask = T) %>%
  terra::app(fun = mean)
vaporS <- terra::rast(vaporlist) %>%
  terra::crop(sag, mask = T) %>%
  terra::app(fun = mean)

# Population density (2020)
popM <- terra::rast("output/TRUE/pop/gpw_v4_population_density_rev11_2020_30s.tif") %>%
  terra::crop(mich, mask = T)
popS <- rast("output/TRUE/pop/gpw_v4_population_density_rev11_2020_30s.tif") %>%
  terra::crop(sag, mask = T)

# Elevation
elevM <- terra::rast("output/TRUE/wc2.1_30s/wc2.1_30s_elev.tif") %>%
  terra::crop(mich, mask = T) 
elevS <- terra::rast("output/TRUE/wc2.1_30s/wc2.1_30s_elev.tif") %>%
  terra::crop(sag, mask = T)

# Nitrogen at soil
nitrogenM <- terra::rast("output/TRUE/soil_world/nitrogen_0-5cm_mean_30s.tif") %>%
  terra::crop(mich, mask = T)
nitrogenS <- terra::rast("output/TRUE/soil_world/nitrogen_0-5cm_mean_30s.tif") %>%
  terra::crop(sag, mask = T)

# Phosphorus fertilizer application (resample to match resolution and projection)
phosphorusM <- raster::raster("output/phosphorus0/Global_Phosphorus_Dist_Map_1223/data/pforms_den.nc") %>%
  raster::projectRaster(to = wcM[[1]], method = "bilinear") %>% #choose another layer with appropriate attributes
  terra::rast() %>%
  terra::crop(mich, mask = T)

phosphorusS <- raster::raster("output/phosphorus0/Global_Phosphorus_Dist_Map_1223/data/pforms_den.nc",
                      var = "tot") %>%
  raster::projectRaster(to = wcS[[1]], method = "bilinear") %>%
  terra::rast() %>%
  terra::crop(sag, mask = T)

# Land cover
# categorical variables cannot be simply project to a new resolution
# but we can use a resampling approach with the exactextractr package
# Reference: https://cran.r-project.org/web/packages/exactextractr/vignettes/vig2_categorical.html

# First, empty raster objects for a template
geodata_path("output")
wc0 <- geodata::worldclim_global(var = "bio", res = 0.5, version = "2.1")
# wc0 is just a template here 

emptyRaster <- terra::rast(ext(wc0), resolution = res(wc0))
crs(emptyRaster) <- crs(wc0)

# empty rasters for each study area
michRaster <- terra::rasterize(mich, emptyRaster) %>%
  terra::crop(mich, mask = T)

sagRaster <- terra::rasterize(sag, emptyRaster) %>%
  terra::crop(sag, mask = T)

michRaster
sagRaster

# Convert those rasters to simple features
michSF <- rasterToPolygons(raster(michRaster), n = 4, na.rm = FALSE) %>%
  st_as_sf()

sagSF <- rasterToPolygons(raster(sagRaster), n = 4, na.rm = FALSE) %>%
  st_as_sf()

michSF
sagSF

# Note the simple features are larger than the actual mich extent
# because of the conversion to sf, which doesn't preserve cropping here
# Layers relying on the simple features will be cropped, though

# Import each landcover raster layer and resample it by finding the majority class
# aka mode in each polygon of the template raster
landcoverM <- raster::raster("output/michLandcover/nlcd_NLCD_Land_Cover_2019.tif")
landcoverM_resampled <- exact_extract(landcoverM, michSF, "mode")

landcoverM_dataframe <- michSF %>%
  mutate(LandClass = landcoverM_resampled)
# And rasterize to the original raster so we have the kind of layer we need

landcoverM_raster <- fasterize(sf = landcoverM_dataframe, raster = raster(michRaster),
                               field = "LandClass") %>%
  terra::rast() %>%
  terra::crop(mich, mask = T)

# Visualize and check out details
plot(landcoverM_raster)
landcoverM_raster
# Reference:
wcM
plot(wcM[[1]])

# Same protocol for Saginaw Bay
landcoverS <- raster::raster("output/sagLandcover/nlcd_NLCD_Land_Cover_2019.tif")
landcoverS_resampled <- exact_extract(landcoverS, sagSF, "mode")

landcoverS_dataframe <- sagSF %>%
  mutate(LandClass = landcoverS_resampled)
# And rasterize to the original raster so we have the kind of layer we need

landcoverS_raster <- fasterize(sf = landcoverS_dataframe, raster = raster(sagRaster),
                               field = "LandClass") %>%
  terra::rast() %>%
  terra::crop(sag, mask = T)

# Visualize and check out details
plot(landcoverS_raster)
landcoverS_raster
# Reference:
wcS
plot(wcS[[1]])

# Distance from nearest boat launch
boat_launch_distM <- terra::rast("output/boat_launch_distM0.tif") %>%
  terra::crop(mich, mask = T)
boat_launch_distS <- terra::rast("output/boat_launch_distS0.tif") %>%
  terra::crop(sag, mask = T)

################################################################################

###################### APPLY A QUICK CHECK ON ALL LAYERS #######################
plot(c(boat_launch_distM, elevM, landcoverM_raster, nitrogenM, phosphorusM,
           popM, solarM, vaporM, wcM, wind_speedM))

plot(c(boat_launch_distS, elevS, landcoverS_raster, nitrogenS, phosphorusS,
           popS, solarS, vaporS, wcS, wind_speedS))

################################################################################

rm(mich)
rm(sag)

############# CHECK FOR CORRELATION BETWEEN EXPLANATORY VARIABLES ##############

# First will be using a correlation threshold of 0.7 (stackMich and stackSag)
# Second will be using a correlation threshold of 0.5 (stackMichB and stackSagB)

# Not including land cover because it is explanatory, but remember to add it to final stacks

# Michigan
stackMich0 <- c(boat_launch_distM0,
                elevM,
                nitrogenM,
                phosphorusM,
                popM,
                solarM,
                vaporM,
                wcM$wc2.1_30s_bio_6,
                wcM$wc2.1_30s_bio_10,
                wcM$wc2.1_30s_bio_11,
                wcM$wc2.1_30s_bio_18,
                wind_speedM) %>%
  setNames(c("boat", "elev", "nit", "phos",
             "pop", "solar", "vapor", "wc6",
             "wc10", "wc11", "wc18", "wind"))

stackMich0_df <- terra::extract(stackMich0, michSF, xy = FALSE) %>% dplyr::select(-ID)
# Many values will be NA because the michSF object covers a larger area
# than the actual Michigan study area
# This is not an issue though as we can just drop the na rows below
# (which is necessary for a correlation matrix anyway)

michCor <- stats::cor(stackMich0_df %>% tidyr::drop_na())

# corrplot::corrplot(michCor, method = "circle", addCoef.col = "black", type = "upper")

# This will return variables to remove at cutoff of 0.7
caret::findCorrelation(michCor, cutoff = 0.7, names = TRUE, exact = TRUE)

stackMich <- terra::subset(stackMich0, c("wc10", "vapor", "wc11", "wc6"),
                           negate = TRUE) %>% append(landcoverM_raster)

# plot(stackMich)

# This will return variables to remove at cutoff of 0.5
caret::findCorrelation(michCor, cutoff = 0.5, names = TRUE, exact = TRUE)

stackMichB <- terra::subset(stackMich0, c("wc10", "vapor", "wc11", "wc6",
                                         "solar", "wc18", "wind"),
                            negate = TRUE) %>% append(landcoverM_raster)

# plot(stackMichB)


# Saginaw Bay
stackSag0 <- c(boat_launch_distS,
               elevS,
               nitrogenS,
               phosphorusS,
               popS,
               solarS,
               vaporS,
               wcS$wc2.1_30s_bio_6,
               wcS$wc2.1_30s_bio_10,
               wcS$wc2.1_30s_bio_11,
               wcS$wc2.1_30s_bio_18,
               wind_speedS) %>%
  setNames(c("boat", "elev", "nit", "phos",
             "pop", "solar", "vapor", "wc6",
             "wc10", "wc11", "wc18", "wind"))

stackSag0_df <- terra::extract(stackSag0, sagSF, xy = FALSE) %>% dplyr::select(-ID)

stackSag0_df %>% tidyr::drop_na() %>% count()

sagCor <- stats::cor(stackSag0_df %>% tidyr::drop_na())

# corrplot::corrplot(sagCor, method = "circle", addCoef.col = "black", type = "upper")

# This will return variables to remove at cutoff of 0.7
caret::findCorrelation(sagCor, cutoff = 0.7, names = TRUE, exact = TRUE)

stackSag <- terra::subset(stackSag0, c("phos", "wc6", "wc18", "wc11",
                                       "solar", "vapor", "wind"),
                           negate = TRUE) %>% append(landcoverS_raster)
# plot(stackSag)

# This will return variables to remove at cutoff of 0.5
caret::findCorrelation(sagCor, cutoff = 0.5, names = TRUE, exact = TRUE)
# returns same variables as the 0.7 cutoff, so no need to make another stack

################################################################################

# Write out stacks
terra::writeRaster(stackMich, "output/stackMich.tif")
terra::writeRaster(stackSag, "output/stackSag.tif")

terra::writeRaster(stackMichB, "output/stackMichB.tif")

rm(list = ls())
