# Data-centric SDM Step 9) Delineating the large-scale study region
## Sara Hansen
## Modified July 25, 2024

library(conflicted)
library(tidyverse)
library(terra)
library(raster) # used only for a few select functions
library(sf)
library(tidyterra)
library(concaveman)
library(USA.state.boundaries)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(enmSdmX)
library(geodata)
library(FedData)
library(jsonlite)
library(ncdf4)
library(ggspatial)
library(ggpubr)

library(tools)
library(mgcv)
library(earth)
library(dismo)
# To use the dismo::maxent function, you need to have the MaxEnt software installed
# https://biodiversityinformatics.amnh.org/open_source/maxent/
# You should also make sure Java is installed and compatible with your version of R
# https://www.java.com/en/download/
library(rJava)
library(gbm)
library(randomForest)
library(party)
library(maxnet)
library(neuralnet)
library(xgboost)
library(caret)
library(Metrics)

# Note that nothing except two "metrics" files, "alt" supplementary table, and "alt" figures
# should be written out. Some objects have the same names as existing objects 
# from earlier steps. Do not write out anything additional without changing names.

################################################################################

# The following code outlines the process for choosing a large scale study area,
# either Michigan or North America. This is a prerequisite to the entire study,
# outlined as Decision Point 3a in the manuscript.
# However, because we go through the entire modeling process here,
# it is easiest to run this step after all other steps, so as not to 
# complicate the linear modeling process too much. 
# So we treat this code as a "final" step for clarity.

################################################################################

# Step 1: Data Download
# Already completed at global scale for all 11 variables to be included
# Note that land cover and boat launch are excluded because they are US/Michigan-specific

# Michigan extent/scale has been defined:
mich <- sf::read_sf("output/mich.shp")

# However, we do need to define the bigger invaded extent
# Which we will call "North America" because it's the entire core invaded area
# of the continent

# Study area = North America
# The spatial scale of explanatory variables will be 30 arc seconds, or about 1 km
# so we will limit our uncertainty to that amount
# We will also exclude the points used in the Saginaw Bay study area

h <- 1000 # threshold (in meters) for excluding records

# We will also exclude a few outliers from the core invaded region

noramPoints <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                          sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::filter(occurrenceStatus == "PRESENT" & 
                  stateProvince != "Washington" & stateProvince != "Florida" & 
                  # ^ outliers: not confirmed wild populations, way outside the normal range
                  !is.na(decimalLongitude) & !is.na(decimalLatitude) &
                  (is.na(coordinateUncertaintyInMeters)
                   | coordinateUncertaintyInMeters <= h) &
                  decimalLongitude > -90) %>%
  distinct(decimalLongitude, decimalLatitude) %>%
  dplyr::select(x = decimalLongitude, y = decimalLatitude)

noramPoints2 <- terra::vect(noramPoints, geom = c("x", "y"))
crs(noramPoints2) <- "+proj=longlat"

# Define a concave hull around points to establish the sampled area
noramHull <- concaveman::concaveman(noramPoints2 %>% 
                                      sf::st_as_sf()) %>%
  terra::vect()

#plot(noramHull)
#plot(noramPoints2, add = TRUE)

# Then establish the 1-km buffer as the study area
noram0 <- terra::buffer(noramHull, width = 1000)
crs(noram0) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

#plot(noram0)

# Finally, create an inner buffer on the Laurentian Great Lakes
# so that the study area does not extend more than 1km into the open water
# and use it to clip the noram study area

# Lakes layer
lakes <- rnaturalearth::ne_download(scale = 110, 
                                    type = "lakes", 
                                    category = "physical",
                                    returnclass = "sf")
crs(lakes) # same as other layers

lakesNegativeBuffer <- lakes %>% terra::vect() %>% terra::buffer(width = -1000)
crs(lakesNegativeBuffer) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

# Final North America study area:
noram <- noram0 %>% terra::erase(lakesNegativeBuffer) %>% st_as_sf()
#plot(noram)

################################################################################

rm(list = ls(pattern = "lakes|noram0|Hull|Points"))

################################################################################

# Step 2: Explanatory Data Import and Processing
#Import the downloaded data and crop to each study area

# WorldClim bioclimatic variables 6, 10, 11, and 18
wclist <- c("output/TRUE/wc2.1_30s/wc2.1_30s_bio_6.tif",
            "output/TRUE/wc2.1_30s/wc2.1_30s_bio_10.tif",
            "output/TRUE/wc2.1_30s/wc2.1_30s_bio_11.tif",
            "output/TRUE/wc2.1_30s/wc2.1_30s_bio_18.tif")
wcM <- terra::rast(wclist) %>%
  terra::crop(mich, mask = T)

wcN <- terra::rast(wclist) %>%
  terra::crop(noram, mask = T)

# Solar radiation (mean across entire year)
solarlist <- list.files(path = "output/TRUE/wc2.1_30s", pattern = "srad",
                        full.names = TRUE) #should be 12
solarM <- terra::rast(solarlist) %>%
  terra::crop(mich, mask = T) %>%
  terra::app(fun = mean)
solarN <- terra::rast(solarlist) %>%
  terra::crop(noram, mask = T) %>%
  terra::app(fun = mean)

# Wind speed (mean across entire year)
wind_speedlist <- list.files(path = "output/TRUE/wc2.1_30s", pattern = "wind",
                             full.names = TRUE)
wind_speedM <- terra::rast(wind_speedlist) %>%
  terra::crop(mich, mask = T) %>%
  terra::app(fun = mean)
wind_speedN <- terra::rast(wind_speedlist) %>%
  terra::crop(noram, mask = T) %>%
  terra::app(fun = mean)

# Water vapor pressure (mean across entire year)
vaporlist <- list.files(path = "output/TRUE/wc2.1_30s", pattern = "vapr",
                        full.names = TRUE)
vaporM <- terra::rast(vaporlist) %>%
  terra::crop(mich, mask = T) %>%
  terra::app(fun = mean)
vaporN <- terra::rast(vaporlist) %>%
  terra::crop(noram, mask = T) %>%
  terra::app(fun = mean)

# Population density (2020)
popM <- terra::rast("output/TRUE/pop/gpw_v4_population_density_rev11_2020_30s.tif") %>%
  terra::crop(mich, mask = T)
popN <- rast("output/TRUE/pop/gpw_v4_population_density_rev11_2020_30s.tif") %>%
  terra::crop(noram, mask = T)

# Elevation
elevM <- terra::rast("output/TRUE/wc2.1_30s/wc2.1_30s_elev.tif") %>%
  terra::crop(mich, mask = T) 
elevN <- terra::rast("output/TRUE/wc2.1_30s/wc2.1_30s_elev.tif") %>%
  terra::crop(noram, mask = T)

# Nitrogen at soil
nitrogenM <- terra::rast("output/TRUE/soil_world/nitrogen_0-5cm_mean_30s.tif") %>%
  terra::crop(mich, mask = T)
nitrogenN <- terra::rast("output/TRUE/soil_world/nitrogen_0-5cm_mean_30s.tif") %>%
  terra::crop(noram, mask = T)

# Phosphorus fertilizer application (resample to match resolution and projection)
phosphorusM <- raster::raster("output/phosphorus0/Global_Phosphorus_Dist_Map_1223/data/pforms_den.nc",
                              var = "tot") %>%
  raster::projectRaster(to = wcM[[1]], method = "bilinear") %>% #choose another layer with appropriate attributes
  terra::rast() %>%
  terra::crop(mich, mask = T)

phosphorusN <- raster::raster("output/phosphorus0/Global_Phosphorus_Dist_Map_1223/data/pforms_den.nc",
                              var = "tot") %>%
  raster::projectRaster(to = wcN[[1]], method = "bilinear") %>%
  terra::rast() %>%
  terra::crop(noram, mask = T)
 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Now define the raster stacks for both of these, at both thresholds
stackMich0 <- c(elevM,
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
  setNames(c("elev", "nit", "phos",
             "pop", "solar", "vapor", "wc6",
             "wc10", "wc11", "wc18", "wind"))

stackMich0_df <- terra::extract(stackMich0, mich, xy = FALSE) %>% dplyr::select(-ID)

michCor <- stats::cor(stackMich0_df %>% tidyr::drop_na())

# corrplot::corrplot(michCor, method = "circle", addCoef.col = "black", type = "upper")

# This will return variables to remove at cutoff of 0.7
caret::findCorrelation(michCor, cutoff = 0.7, names = TRUE, exact = TRUE)

stackMich <- terra::subset(stackMich0, c("wc10", "vapor", "wc11", "wc6"),
                           negate = TRUE)

# Cutoff of 0.5
caret::findCorrelation(michCor, cutoff = 0.5, names = TRUE, exact = TRUE)

stackMichB <- terra::subset(stackMich0, c("wc10", "vapor", "wc11", "wc6",
                                         "solar", "wc18", "wind"),
                            negate = TRUE)

# And the same for North America
stackNoram0 <- c(elevN,
                 nitrogenN,
                 phosphorusN,
                 popN,
                 solarN,
                 vaporN,
                 wcN$wc2.1_30s_bio_6,
                 wcN$wc2.1_30s_bio_10,
                 wcN$wc2.1_30s_bio_11,
                 wcN$wc2.1_30s_bio_18,
                 wind_speedN) %>%
  setNames(c("elev", "nit", "phos",
             "pop", "solar", "vapor", "wc6",
             "wc10", "wc11", "wc18", "wind"))

stackNoram0_df <- terra::extract(stackNoram0, noram, xy = FALSE) %>% dplyr::select(-ID)

noramCor <- stats::cor(stackNoram0_df %>% tidyr::drop_na())

# corrplot::corrplot(michCor, method = "circle", addCoef.col = "black", type = "upper")

caret::findCorrelation(noramCor, cutoff = 0.7, names = TRUE, exact = TRUE)

stackNoram <- terra::subset(stackNoram0, c("wc10", "vapor", "wc11", "wc6"),
                            negate = TRUE)


caret::findCorrelation(noramCor, cutoff = 0.5, names = TRUE, exact = TRUE)

stackNoramB <- terra::subset(stackNoram0, c("wc10", "vapor", "wc11", "wc6",
                                           "solar", "wc18"), 
                            negate = TRUE)

################################################################################

rm(list = ls(pattern = "elev|nitrogen|phosphorus|pop|solar|vapor|wc|wind"))
rm(list = ls(pattern = "Cor|0"))

################################################################################

# Step 3: Response Data Download, Import, and Processing
# We do not need to get the Michigan response data again (but remember to remove some columns)
# But we do need to get the North America response data

efbNoramP0 <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                          sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::filter(occurrenceStatus == "PRESENT" & 
                  stateProvince != "Washington" & stateProvince != "Florida" & 
                  # ^ outliers: not confirmed wild populations, way outside the normal range
                  !is.na(decimalLongitude) & !is.na(decimalLatitude) &
                  (is.na(coordinateUncertaintyInMeters)
                   | coordinateUncertaintyInMeters <= h) &
                  decimalLongitude > -90) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  mutate(occurrence = as.integer(1),
         id = row_number()) %>%
  dplyr::select(x = decimalLongitude,
                y = decimalLatitude, occurrence, id)

# Use gridSample to thin occurrences
# and remove any occurrences outside the defined study region
# To find ones outside of the study region, we can extract values of variables
# and any points with only NA are outside the region
# Use the noram stack with more variables
names(stackNoram)
nonNApoints <- bind_cols(efbNoramP0, terra::extract(stackNoram, efbNoramP0[, 1:2])) %>%
  dplyr::select(-ID) %>%
  dplyr::filter(!is.na(elev) | !is.na(nit) | !is.na(phos) | !is.na(pop) |
                  !is.na(solar) | !is.na(wc18) | !is.na(wind))
# NA points are excluded because of how we clipped, ie they are on islands or in open water

# Example plots:
ggplot() +
  geom_sf(data = noram) +
  geom_point(data = nonNApoints, aes(x = x, y = y))

ggplot() +
  geom_sf(data = noram) +
  geom_point(data = efbNoramP0 %>% anti_join(nonNApoints, by = "id"),
             aes(x = x, y = y))
# These make < 4% of the points

set.seed(345)
efbNoramP_gs <- dismo::gridSample(nonNApoints %>% 
                                   dplyr::select(x, y),
                                 r = raster::raster(stackNoram), n = 1) %>% # n = 1334
  mutate(occurrence = as.integer(1),
         id = row_number())

efbNoramP <- efbNoramP_gs 

################################################################################

rm(efbNoramP_gs, efbNoramP0)

################################################################################

# Step 4: Get Background Points

# Ratio is 1:1 presence to absences, just as with original models

# Get raster for noram study area
noramRaster <- terra::subset(stackNoram, 1) %>% raster::raster() # can be any of the layers

drawBackground <- function(i, extentOb) {
  
  tempNum <- as.numeric(count(i))
  
  set.seed(789)
  
  background <- dismo::randomPoints(mask = extentOb, n = tempNum,
                                    p = i %>% dplyr::select(x, y),
                                    excludep = TRUE) %>%
    as.data.frame() %>%
    mutate(occurrence = as.integer(0),
           id = paste0("background_", row_number()))
  
  rbind(i, background) %>%
    sample_frac(1)
}

efbNoramPB <- drawBackground(efbNoramP, noramRaster)

# See that no two points are in the same cell (n is always 1):
#terra::cells(rast(noramRaster), vect(efbNoramPB, geom = c("x", "y"))) %>%
  #as.data.frame() %>% count(cell) %>% summarize(max(n))

# Extract explanatory variables at point locations
extractEx <- function(i, stack) {
  bind_cols(i, terra::extract(stack, i[, 1:2])) %>%
    dplyr::select(-ID)
  }

noramPB0 <- extractEx(efbNoramPB, stackNoram)

# And scale:
scaleVars <- function(i) {
  
  r <- ncol(i)
  
  i %>%
    dplyr::mutate(across(5:all_of(r),
                         ~ c(scale(.)))) #c() makes output a data frame, not a series of matrices
}

noramPB <- scaleVars(noramPB0)

# Do this again for the 0.5 findCorrelation threshold "B" dataset
noramPB_B0 <- extractEx(efbNoramPB, stackNoramB)
noramPB_B <- scaleVars(noramPB_B0)

# And the presence datasets
noramP0 <- extractEx(efbNoramP, stackNoram)
noramP <- scaleVars(noramP0)

noramP_B0 <- extractEx(efbNoramP, stackNoramB)
noramP_B <- scaleVars(noramP_B0)

# Because North America was not ultimately chosen for further decision points,
# we do not need the data parsed as "all," specimen, and observation
# Optional code is below for parsing and creating the relevant datasets

# Specimens
#efbNoramP0_spec <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                         #sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  #dplyr::filter(occurrenceStatus == "PRESENT" & 
                  #stateProvince != "Washington" & stateProvince != "Florida" & 
                  #!is.na(decimalLongitude) & !is.na(decimalLatitude) &
                  #(is.na(coordinateUncertaintyInMeters)
                   #| coordinateUncertaintyInMeters <= h) &
                  #decimalLongitude > -90) %>%
  #dplyr::filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
  #distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  #mutate(occurrence = as.integer(1),
         #id = row_number()) %>%
  #dplyr::select(x = decimalLongitude,
                #y = decimalLatitude, occurrence, id) %>%
  #inner_join(nonNApoints %>% dplyr::select(x, y), by = c("x", "y"))

#efbNoramP_spec <- dismo::gridSample(efbNoramP0_spec %>% 
                                      #dplyr::select(x, y),
                                  #r = raster::raster(stackNoram), n = 1) %>% # n = 200
  #mutate(occurrence = as.integer(1),
         #id = row_number())

#efbNoramPB_spec <- drawBackground(efbNoramP_spec, noramRaster)
#noramPB0_spec <- extractEx(efbNoramPB_spec, stackNoram)
#noramPB_spec <- scaleVars(noramPB0_spec)

#noramPB_B0_spec <- extractEx(efbNoramPB_spec, stackNoramB)
#noramPB_B_spec <- scaleVars(noramPB_B0_spec)

#noramP0_spec <- extractEx(efbNoramP_spec, stackNoram)
#noramP_spec <- scaleVars(noramP0_spec)

#noramP_B0_spec <- extractEx(efbNoramP_spec, stackNoramB)
#noramP_B_spec <- scaleVars(noramP_B0_spec)

# Observation
#efbNoramP0_obs <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                              #sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  #dplyr::filter(occurrenceStatus == "PRESENT" & 
                  #stateProvince != "Washington" & stateProvince != "Florida" & 
                  #!is.na(decimalLongitude) & !is.na(decimalLatitude) &
                  #(is.na(coordinateUncertaintyInMeters)
                   #| coordinateUncertaintyInMeters <= h) &
                  #decimalLongitude > -90) %>%
  #dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION") %>%
  #distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  #mutate(occurrence = as.integer(1),
         #id = row_number()) %>%
  #dplyr::select(x = decimalLongitude,
                #y = decimalLatitude, occurrence, id) %>%
  #inner_join(nonNApoints %>% dplyr::select(x, y), by = c("x", "y"))

#efbNoramP_obs <- dismo::gridSample(efbNoramP0_obs %>% 
                                      #dplyr::select(x, y),
                                    #r = raster::raster(stackNoram), n = 1) %>% # n = 200
  #mutate(occurrence = as.integer(1),
         #id = row_number())

#efbNoramPB_obs <- drawBackground(efbNoramP_obs, noramRaster)
#noramPB0_obs <- extractEx(efbNoramPB_obs, stackNoram)
#noramPB_obs <- scaleVars(noramPB0_obs)

#noramPB_B0_obs <- extractEx(efbNoramPB_obs, stackNoramB)
#noramPB_B_obs <- scaleVars(noramPB_B0_obs)

#noramP0_obs <- extractEx(efbNoramP_obs, stackNoram)
#noramP_obs <- scaleVars(noramP0_obs)

#noramP_B0_obs <- extractEx(efbNoramP_obs, stackNoramB)
#noramP_B_obs <- scaleVars(noramP_B0_obs)

# For Michigan, we can read in existing data and just exclude extra variables (see below)
  
################################################################################

rm(list = ls(pattern = "efb|0|Raster"))

################################################################################

# Step 5: Fit and Evaluate Models

# Create functions (exactly the same as in Step 5, just no "land")
# (converMetrics has modifications to assign correct grouping
# Do both the "row remove" and "column remove" methods, even though the 0.5
# threshold datasets are really small and will only have a couple variables left
# Just to be consistent with earlier methods

# Exact same functions as Step 5:
# Function to build presence-absence (or presence-background) models
buildPAmodels <- function(i, gam = TRUE, brtRate = 0.01, bagFraction = 0.75) {

  ### Data prep 
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  
  rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  # Datasets with NA columns removed
  set.seed(987)
  colrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
    dplyr::select(where(~ !any(is.na(.)))) %>% 
    dplyr::select(-c(x, y, id))
  
  colrm_train1 <- colrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  colrm_train2 <- colrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  colrm_train3 <- colrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  colrm_train4 <- colrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  datList <- list(rowrm_train1 = rowrm_train1,
                  rowrm_train2 = rowrm_train2,
                  rowrm_train3 = rowrm_train3,
                  rowrm_train4 = rowrm_train4,
                  colrm_train1 = colrm_train1,
                  colrm_train2 = colrm_train2,
                  colrm_train3 = colrm_train3,
                  colrm_train4 = colrm_train4)
  
  lapply(datList, FUN = function(dat) {
    
      maxPred <- dat %>% dplyr::select(-occurrence)
      
      annDat <- dat %>% as.data.frame()
      
      xgDat <- annDat %>% as.matrix() 
      
      gamForm <- as.formula(paste0("occurrence ~", paste(
        paste0("s(", colnames(dat)[-1], ", k = 8)"), collapse = " + ")))
    
    annForm <- as.formula(paste0("occurrence ~ ", paste(colnames(annDat)[-1], collapse = " + ")))
    
    r <- ncol(dat)
    
    q <- ncol(xgDat)
    
    ### Models
    print(paste0("USING ", nrow(dat), " ROWS AND ", (ncol(dat)-1), " PREDICTORS"))
    
    set.seed(654)
    # Generalized Linear Model (GLM)
    glm = stats::glm(occurrence ~ ., family = binomial(link = "logit"), data = dat)
    
    set.seed(127)
    # Generalized Additive Model (GAM)
    gam = 
      if(gam == TRUE) {
        mgcv::gam(gamForm, family = binomial(link = "logit"), data = dat, method = "REML") }
    
    else if(gam == FALSE) { NULL }
    
    set.seed(996)
    # Multivariate Adaptive Regression Splines (MARS)
    mars = earth::earth(occurrence ~ ., data = dat)
    
    set.seed(125)
    # Boosted Regression Trees (BRT) aka Gradient Boosting Machine
    brt = dismo::gbm.step(data = as.data.frame(dat),
                          gbm.x = 2:r, gbm.y = "occurrence",
                          family = "bernoulli", #binomial
                          tree.complexity = 5, learning.rate = brtRate,
                          bag.fraction = bagFraction, max.trees = 10000,
                          plot.main = FALSE, verbose = FALSE)
    
    set.seed(424)
    # Random Forest (RF)
    rf = randomForest::randomForest(as.factor(occurrence) ~ ., data = dat,
                                    importance = TRUE, ntree = 5000, replace = TRUE)
    
    set.seed(193)
    # Conditional Forest (Cforest)
    cforest = party::cforest(as.factor(occurrence) ~ ., data = dat,
                             controls = cforest_unbiased(ntree = 5000, mtry = 5))
    
    set.seed(368)
    # MaxEnt
    maxent = dismo::maxent(x = maxPred, p = dat$occurrence)
    
    set.seed(705)
    # MaxNet
    maxnet = maxnet::maxnet(p = dat$occurrence, data = maxPred,
                            f = maxnet.formula(dat$occurrence, maxPred))
    
    set.seed(678)
    # Artificial Neural Network (ANN)
    ann =  neuralnet(annForm, data = annDat, hidden = c(5, 2),
                     linear.output = FALSE, lifesign = "minimal", startweights = NULL)
    
    set.seed(512)
    # eXtreme Gradient Boosting Machine (XGBoost)
    xgboost = xgboost::xgboost(data = xgDat[, 2:q], label = xgDat[, 1],
                               nrounds = 5000, max_depth = 10, eta = 0.8,
                               objective = "binary:logistic", verbose = 0)
    
    list(glm = glm, gam = gam, mars = mars, brt = brt, rf = rf, cforest = cforest, 
         maxent = maxent, maxnet = maxnet, ann = ann, xgboost = xgboost) # close the list
  } # close the mini-function
  ) # close the lapply
} # close the overall function

# Function to evaluate presence-absence models
evaluatePAmodels <- function(i, modList, threshold) {

  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  
  rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  rowrm_test1 <- rowrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group)
  rowrm_test2 <- rowrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group)
  rowrm_test3 <- rowrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group)
  rowrm_test4 <- rowrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group)
  
  # Datasets with NA columns removed
  set.seed(987)
  colrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
    dplyr::select(where(~ !any(is.na(.)))) %>% 
    dplyr::select(-c(x, y, id))
  
  colrm_train1 <- colrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  colrm_train2 <- colrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  colrm_train3 <- colrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  colrm_train4 <- colrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  colrm_test1 <- colrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group)
  colrm_test2 <- colrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group)
  colrm_test3 <- colrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group)
  colrm_test4 <- colrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group)
  
  lapply(1:length(modList),
         FUN = function(mod) {
           
           # Setting which data to use to test each model
           if(str_detect(names(modList[mod]), "row") &
              str_detect(names(modList[mod]), "1")) { testDat <- rowrm_test1 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "2")) { testDat <- rowrm_test2 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "3")) { testDat <- rowrm_test3 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "4")) { testDat <- rowrm_test4 }
           
           else if(str_detect(names(modList[mod]), "col") &
                   str_detect(names(modList[mod]), "1")) { testDat <- colrm_test1 }
           
           else if(str_detect(names(modList[mod]), "col") &
                   str_detect(names(modList[mod]), "2")) { testDat <- colrm_test2 }
           
           else if(str_detect(names(modList[mod]), "col") &
                   str_detect(names(modList[mod]), "3")) { testDat <- colrm_test3 }
           
           else if(str_detect(names(modList[mod]), "col") &
                   str_detect(names(modList[mod]), "4")) { testDat <- colrm_test4 }
           
             maxTestDat <- testDat
             
             annTestDat <- testDat
             
             xgTestDat <- testDat %>% as.matrix()
           
           set.seed(456)
           
           # Testing models based on their type
           if(str_detect(names(modList[mod]), "glm|gam|mars|brt")) {
             
             if(is.null(class(modList[mod]))) {predictDF <- NULL}
             else{
               predictDF <- stats::predict(object = modList[[mod]],
                                           newdata = testDat[, -1],
                                           type = "response") %>%
                 as.data.frame() %>%
                 dplyr::select(prediction = 1) } }
           
           else if(str_detect(names(modList[mod]), "rf")) {
             predictDF <- stats::predict(object = modList[[mod]],
                                         newdata = testDat[, -1],
                                         type = "prob") %>%
               as.data.frame() %>%
               dplyr::select(prediction = 2) }
           
           else if(str_detect(names(modList[mod]), "cf")) {
             predictDF <- stats::predict(object = modList[[mod]],
                                         newdata = testDat[, -1],
                                         type = "prob") %>%
               tibble::enframe(name = NULL,
                               value = "value") %>%
               tidyr::separate(col = "value",
                               sep = ",",
                               into = c("zero", "one")) %>%
               dplyr::mutate(one = str_remove_all(one, "\\)"),
                             one = as.numeric(one)) %>%
               dplyr::select(prediction = one) }
           
           else if(str_detect(names(modList[mod]), "maxent")) {
             predictDF <- dismo::predict(object = modList[[mod]],
                                         x = maxTestDat[, -1]) %>%
               as.data.frame() %>%
               dplyr::select(prediction = 1) }
           
           else if(str_detect(names(modList[mod]), "maxnet")) {
             predictDF <- stats::predict(object = modList[[mod]],
                                         newdata = maxTestDat[, -1],
                                         type = "cloglog") %>%
               as.data.frame() %>%
               dplyr::select(prediction = 1) }
           
           else if(str_detect(names(modList[mod]), "ann")) {
             predictDF <- stats::predict(object = modList[mod], # single bracket subset
                                         newdata = annTestDat[, -1],
                                         type = "response") %>%
               as.data.frame() %>%
               dplyr::select(prediction = 1) }
           
           else if(str_detect(names(modList[mod]), "xgboost")) {
             predictDF <- stats::predict(object = modList[[mod]],
                                         newdata = xgTestDat[, -1],
                                         type = "response") %>%
               as.data.frame() %>%
               dplyr::select(prediction = 1) }
           
           # Prediction vector
           prediction  <- as.vector(predictDF$prediction)
           
           # Defining the table to calculate metrics from        
           if(str_detect(names(modList[mod]), "glm|gam|mars|brt|rf|cf")) {
             
             m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                             truth = as.numeric(testDat$occurrence)) }
           
           else if(str_detect(names(modList[mod]), "max")) {
             
             m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                             truth = as.numeric(maxTestDat$occurrence))  }
           
           else if(str_detect(names(modList[mod]), "ann")) {
             
             m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                             truth = as.numeric(annTestDat$occurrence))  }
           
           else if(str_detect(names(modList[mod]), "xgboost")) {
             
             m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                             truth = as.numeric(xgTestDat[, 1]))  }
           
           list(prediction = prediction,
                
                metrics = data.frame(c(
                  TN = m %>% dplyr::filter(prediction == 0 & truth == 0) %>% count(),
                  FP = m %>% dplyr::filter(prediction == 1 & truth == 0) %>% count(),
                  FN = m %>% dplyr::filter(prediction == 0 & truth == 1) %>% count(),
                  TP = m %>% dplyr::filter(prediction == 1 & truth == 1) %>% count())) %>%
                  rename_with(~str_remove(., '.n')) %>%
                  mutate(sensitivity = TP/(TP+FN),
                         specificity = TN/(TN+FP),
                         threshold = threshold,
                         mae = Metrics::mae(actual = m$truth,
                                            predicted = prediction),
                         # note that AUC ROC is not reported in the manuscript
                         aucROC = mltools::auc_roc(preds = prediction,
                                                   actuals = m$truth)))
         } # close the internal function
  ) # close the lapply
} # close the overall function

# Function to build presence-only models
buildPmodels <- function(i) {
  
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, occurrence, id))
  
  rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  # Datasets with NA columns removed
  set.seed(987)
  colrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
    dplyr::select(where(~ !any(is.na(.)))) %>% 
    dplyr::select(-c(x, y, occurrence, id))
  
  colrm_train1 <- colrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  colrm_train2 <- colrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  colrm_train3 <- colrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  colrm_train4 <- colrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  datList <- list(rowrm_train1 = rowrm_train1,
                  rowrm_train2 = rowrm_train2,
                  rowrm_train3 = rowrm_train3,
                  rowrm_train4 = rowrm_train4,
                  colrm_train1 = colrm_train1,
                  colrm_train2 = colrm_train2,
                  colrm_train3 = colrm_train3,
                  colrm_train4 = colrm_train4)
  
  lapply(datList, FUN = function(dat) {
    
    r <- ncol(dat)
    
    print(paste0("USING ", nrow(dat), " ROWS AND ", (ncol(dat)-1), " PREDICTORS"))
    
    set.seed(656)
    # BIOCLIM
    bioclim = dismo::bioclim(dat)
    
    set.seed(123)
    # DOMAIN
    domain = dismo::domain(dat)
    
    list(bioclim = bioclim, domain = domain) # close the list
  } # close the mini-function
  ) # close the lapply
} # close the overall function

# Function to evaluate presence-only models 
evaluatePmodels <- function(i, modList, threshold) {
  
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  
  rowrm_test1 <- rowrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group)
  rowrm_test2 <- rowrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group)
  rowrm_test3 <- rowrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group)
  rowrm_test4 <- rowrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group)
  
  # Datasets with NA columns removed
  set.seed(987)
  colrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
    dplyr::select(where(~ !any(is.na(.)))) %>% 
    dplyr::select(-c(x, y, id))
  
  colrm_test1 <- colrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group)
  colrm_test2 <- colrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group)
  colrm_test3 <- colrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group)
  colrm_test4 <- colrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group)
  
  lapply(1:length(modList),
         FUN = function(mod) {
           
           # Setting which data to use to test each model
           if(str_detect(names(modList[mod]), "row") &
              str_detect(names(modList[mod]), "1")) { testDat <- rowrm_test1 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "2")) { testDat <- rowrm_test2 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "3")) { testDat <- rowrm_test3 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "4")) { testDat <- rowrm_test4 }
           
           else if(str_detect(names(modList[mod]), "col") &
                   str_detect(names(modList[mod]), "1")) { testDat <- colrm_test1 }
           
           else if(str_detect(names(modList[mod]), "col") &
                   str_detect(names(modList[mod]), "2")) { testDat <- colrm_test2 }
           
           else if(str_detect(names(modList[mod]), "col") &
                   str_detect(names(modList[mod]), "3")) { testDat <- colrm_test3 }
           
           else if(str_detect(names(modList[mod]), "col") &
                   str_detect(names(modList[mod]), "4")) { testDat <- colrm_test4 }
           
           set.seed(456)
           
           # Testing models
           prediction <- dismo::predict(object = modList[[mod]],
                                        x = testDat[, -1])
           
           # Defining the table to calculate metrics from 
           m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                           truth = as.numeric(testDat$occurrence))
           
           list(prediction = prediction,
                
                metrics = data.frame(c(
                  TN = m %>% dplyr::filter(prediction == 0 & truth == 0) %>% count(),
                  FP = m %>% dplyr::filter(prediction == 1 & truth == 0) %>% count(),
                  FN = m %>% dplyr::filter(prediction == 0 & truth == 1) %>% count(),
                  TP = m %>% dplyr::filter(prediction == 1 & truth == 1) %>% count())) %>%
                  rename_with(~str_remove(., '.n')) %>%
                  mutate(sensitivity = TP/(TP+FN),
                         specificity = TN/(TN+FP),
                         threshold = threshold,
                         mae = Metrics::mae(actual = m$truth,
                                            predicted = prediction)))
         } # close the internal function
  ) # close the lapply
} # close the overall function

#Function to convert evaluation metrics to dataframe
convertMetrics <- function(eval, modList) {
  
  evalName <- deparse(substitute(eval))
  
  dplyr::bind_rows(lapply(eval, `[[`, 2)) %>%
    dplyr::mutate(overall = (sensitivity + specificity + (1-mae)),
                  model = names(modList),
                  modelType = case_when(str_detect(model, "glm") ~ "GLM",
                                        str_detect(model, "gam") ~ "GAM",
                                        str_detect(model, "mars") ~ "MARS",
                                        str_detect(model, "brt") ~ "BRT",
                                        str_detect(model, "rf") ~ "RF",
                                        str_detect(model, "cforest") ~ "Cforest",
                                        str_detect(model, "maxent") ~ "MaxEnt",
                                        str_detect(model, "maxnet") ~ "MaxNet",
                                        str_detect(model, "ann") ~ "ANN",
                                        str_detect(model, "xgboost") ~ "XGBoost",
                                        str_detect(model, "bioclim") ~ "BIOCLIM",
                                        str_detect(model, "domain") ~ "DOMAIN"),
                  removalType = case_when(str_detect(model, "rowrm") ~ "row remove",
                                          str_detect(model, "colrm") ~ "col remove"),
                  scale = case_when(str_detect(evalName, "mich") ~ "mich",
                                    str_detect(evalName, "noram") ~ "noram"),
                  dataSource = "all",
                  dataType = case_when(str_detect(evalName, "PB") ~ "presence-background",
                                       str_detect(evalName, "P") &
                                         !str_detect(evalName, "PB|PA") ~ "presence-only",
                                       str_detect(evalName, "PA") ~ "presence-absence",
                                       str_detect(evalName, "Ab") ~ "abundance"),
                  modelClass = case_when(str_detect(modelType, "BIOCLIM|DOMAIN") ~ "Envelope",
                                         str_detect(modelType, "GLM|GAM|MARS|MaxEnt|MaxNet") ~ "Regression",
                                         str_detect(modelType, "BRT|RF|Cforest|ANN|XGBoost") ~ "Machine learning"))
}

# Now build the set of models for each study area
# Michigan, 0.7 findCorrelation threshold (remove extra variables)
michPB0 <- read.delim("output/michPB.txt", header = TRUE,
                     sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")
colnames(michPB0)
names(stackMich)
michPB <- michPB0 %>%
  dplyr::select(x, y, occurrence, id, elev, nit, phos, pop, solar, wc18, wind)

# Build models
michPBmodels0 <- buildPAmodels(i = michPB)

# Unnest the list for the evaluate functions
michPBmodels <- unlist(michPBmodels0, recursive = FALSE)

michPBmodels[["colrm_train1.ann"]] <- NULL # this one built but not enough in group
michPBmodels[["colrm_train2.ann"]] <- NULL
michPBmodels[["colrm_train3.ann"]] <- NULL
michPBmodels[["colrm_train4.ann"]] <- NULL

# Evaluate at threshold of 0.5
michPBeval <- evaluatePAmodels(i = michPB, modList = michPBmodels, threshold = 0.5)

# Convert metrics to a data frame for saving and analyzing 
michPBmetrics <- convertMetrics(eval = michPBeval, modList = michPBmodels)

# Write out metrics
write_delim(michPBmetrics, "output/michPBmetrics_alt.txt", delim = "\t", quote = "none")

#saveRDS(michPBmodels, "output/michPBmodels_alt.rds")
#saveRDS(michPBeval, "output/michPBeval_alt.rds")

rm(list = ls(pattern = "michPB"))

# Michigan, 0.5 findCorrelation threshold (remove extra variables)
michPB_B0 <- read.delim("output/michPB_B.txt", header = TRUE,
                      sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")
colnames(michPB_B0)
names(stackMichB)
michPB_B <- michPB_B0 %>%
  dplyr::select(x, y, occurrence, id, elev, nit, phos, pop)

# Build models
michPB_Bmodels0 <- buildPAmodels(i = michPB_B)

# Unnest the list for the evaluate functions
michPB_Bmodels <- unlist(michPB_Bmodels0, recursive = FALSE)

michPB_Bmodels[["colrm_train1.ann"]] <- NULL

# Evaluate at threshold of 0.5
michPB_Beval <- evaluatePAmodels(i = michPB_B, modList = michPB_Bmodels, threshold = 0.5)

# Convert metrics to a data frame for saving and analyzing 
michPB_Bmetrics <- convertMetrics(eval = michPB_Beval, modList = michPB_Bmodels)

# Write out metrics
write_delim(michPB_Bmetrics, "output/michPB_Bmetrics_alt.txt", delim = "\t", quote = "none")

#saveRDS(michPB_Bmodels, "output/michPB_Bmodels_alt.rds")
#saveRDS(michPB_Beval, "output/michPB_Beval_alt.rds")

rm(list = ls(pattern = "michPB_B"))

# The two sets of Michigan presence-only models
# Read in data
michP <- read.delim("output/michP.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM") %>%
  dplyr::select(x, y, occurrence, id, elev, nit, phos, pop, solar, wc18, wind)

# Build models
michPmodels0 <- buildPmodels(michP)
michPmodels <- unlist(michPmodels0, recursive = FALSE)

# Evaluate
michPeval <- evaluatePmodels(michP, michPmodels, 0.5)
michPmetrics <- convertMetrics(michPeval, michPmodels)

# Write out metrics
write_delim(michPmetrics, "output/michPmetrics_alt.txt", delim = "\t", quote = "none")

#saveRDS(michPmodels, "output/michPmodels_alt.rds")
#saveRDS(michPeval, "output/michPeval_alt.rds")

rm(list = ls(pattern = "michP"))

# Read in data
michP_B <- read.delim("output/michP_B.txt", header = TRUE,
                      sep = "\t", quote = "", fileEncoding = "UTF-8-BOM") %>%
  dplyr::select(x, y, occurrence, id, elev, nit, phos, pop)

# Build models
michP_Bmodels0 <- buildPmodels(michP_B)
michP_Bmodels <- unlist(michP_Bmodels0, recursive = FALSE)

# Evaluate models
michP_Beval <- evaluatePmodels(michP_B, michP_Bmodels, 0.5)
michP_Bmetrics <- convertMetrics(michP_Beval, michP_Bmodels)

# Write out metrics
write_delim(michP_Bmetrics, "output/michP_Bmetrics_alt.txt", delim = "\t", quote = "none")

#saveRDS(michP_Bmodels, "output/michP_Bmodels_alt.rds")
#saveRDS(michP_Beval, "output/michP_Beval_alt.rds")

rm(list = ls(pattern = "michP_B"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# North America

# Build models
#noramPBmodels0 <- buildPAmodels(i = noramPB)
# All columns are missing at least some data, so we actually need 
# to redefine the modeling functions just without column removal

buildPAmodels_nocolremove <- function(i, gam = TRUE, brtRate = 0.01, bagFraction = 0.75) {
  
  ### Data prep 
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  
  rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

  datList <- list(rowrm_train1 = rowrm_train1,
                  rowrm_train2 = rowrm_train2,
                  rowrm_train3 = rowrm_train3,
                  rowrm_train4 = rowrm_train4)
  
  lapply(datList, FUN = function(dat) {
    
    maxPred <- dat %>% dplyr::select(-occurrence)
    
    annDat <- dat %>% as.data.frame()
    
    xgDat <- annDat %>% as.matrix() 
    
    gamForm <- as.formula(paste0("occurrence ~", paste(
      paste0("s(", colnames(dat)[-1], ", k = 8)"), collapse = " + ")))
    
    annForm <- as.formula(paste0("occurrence ~ ", paste(colnames(annDat)[-1], collapse = " + ")))
    
    r <- ncol(dat)
    
    q <- ncol(xgDat)
    
    ### Models
    print(paste0("USING ", nrow(dat), " ROWS AND ", (ncol(dat)-1), " PREDICTORS"))
    
    set.seed(654)
    # Generalized Linear Model (GLM)
    glm = stats::glm(occurrence ~ ., family = binomial(link = "logit"), data = dat)
    
    set.seed(127)
    # Generalized Additive Model (GAM)
    gam = 
      if(gam == TRUE) {
        mgcv::gam(gamForm, family = binomial(link = "logit"), data = dat, method = "REML") }
    
    else if(gam == FALSE) { NULL }
    
    set.seed(996)
    # Multivariate Adaptive Regression Splines (MARS)
    mars = earth::earth(occurrence ~ ., data = dat)
    
    set.seed(125)
    # Boosted Regression Trees (BRT) aka Gradient Boosting Machine
    brt = dismo::gbm.step(data = as.data.frame(dat),
                          gbm.x = 2:r, gbm.y = "occurrence",
                          family = "bernoulli", #binomial
                          tree.complexity = 5, learning.rate = brtRate,
                          bag.fraction = bagFraction, max.trees = 10000,
                          plot.main = FALSE, verbose = FALSE)
    
    set.seed(424)
    # Random Forest (RF)
    rf = randomForest::randomForest(as.factor(occurrence) ~ ., data = dat,
                                    importance = TRUE, ntree = 5000, replace = TRUE)
    
    set.seed(193)
    # Conditional Forest (Cforest)
    cforest = party::cforest(as.factor(occurrence) ~ ., data = dat,
                             controls = cforest_unbiased(ntree = 5000, mtry = 5))
    
    set.seed(368)
    # MaxEnt
    maxent = dismo::maxent(x = maxPred, p = dat$occurrence)
    
    set.seed(705)
    # MaxNet
    maxnet = maxnet::maxnet(p = dat$occurrence, data = maxPred,
                            f = maxnet.formula(dat$occurrence, maxPred))
    
    set.seed(678)
    # Artificial Neural Network (ANN)
    ann =  neuralnet(annForm, data = annDat, hidden = c(5, 2),
                     linear.output = FALSE, lifesign = "minimal", startweights = NULL)
    
    set.seed(512)
    # eXtreme Gradient Boosting Machine (XGBoost)
    xgboost = xgboost::xgboost(data = xgDat[, 2:q], label = xgDat[, 1],
                               nrounds = 5000, max_depth = 10, eta = 0.8,
                               objective = "binary:logistic", verbose = 0)
    
    list(glm = glm, gam = gam, mars = mars, brt = brt, rf = rf, cforest = cforest, 
         maxent = maxent, maxnet = maxnet, ann = ann, xgboost = xgboost) # close the list
  } # close the mini-function
  ) # close the lapply
} # close the overall function
evaluatePAmodels_nocolremove <- function(i, modList, threshold) {
  
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  
  rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  rowrm_test1 <- rowrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group)
  rowrm_test2 <- rowrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group)
  rowrm_test3 <- rowrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group)
  rowrm_test4 <- rowrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group)
  
  lapply(1:length(modList),
         FUN = function(mod) {
           
           # Setting which data to use to test each model
           if(str_detect(names(modList[mod]), "row") &
              str_detect(names(modList[mod]), "1")) { testDat <- rowrm_test1 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "2")) { testDat <- rowrm_test2 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "3")) { testDat <- rowrm_test3 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "4")) { testDat <- rowrm_test4 }
           
           maxTestDat <- testDat
           
           annTestDat <- testDat
           
           xgTestDat <- testDat %>% as.matrix()
           
           set.seed(456)
           
           # Testing models based on their type
           if(str_detect(names(modList[mod]), "glm|gam|mars|brt")) {
             
             if(is.null(class(modList[mod]))) {predictDF <- NULL}
             else{
               predictDF <- stats::predict(object = modList[[mod]],
                                           newdata = testDat[, -1],
                                           type = "response") %>%
                 as.data.frame() %>%
                 dplyr::select(prediction = 1) } }
           
           else if(str_detect(names(modList[mod]), "rf")) {
             predictDF <- stats::predict(object = modList[[mod]],
                                         newdata = testDat[, -1],
                                         type = "prob") %>%
               as.data.frame() %>%
               dplyr::select(prediction = 2) }
           
           else if(str_detect(names(modList[mod]), "cf")) {
             predictDF <- stats::predict(object = modList[[mod]],
                                         newdata = testDat[, -1],
                                         type = "prob") %>%
               tibble::enframe(name = NULL,
                               value = "value") %>%
               tidyr::separate(col = "value",
                               sep = ",",
                               into = c("zero", "one")) %>%
               dplyr::mutate(one = str_remove_all(one, "\\)"),
                             one = as.numeric(one)) %>%
               dplyr::select(prediction = one) }
           
           else if(str_detect(names(modList[mod]), "maxent")) {
             predictDF <- dismo::predict(object = modList[[mod]],
                                         x = maxTestDat[, -1]) %>%
               as.data.frame() %>%
               dplyr::select(prediction = 1) }
           
           else if(str_detect(names(modList[mod]), "maxnet")) {
             predictDF <- stats::predict(object = modList[[mod]],
                                         newdata = maxTestDat[, -1],
                                         type = "cloglog") %>%
               as.data.frame() %>%
               dplyr::select(prediction = 1) }
           
           else if(str_detect(names(modList[mod]), "ann")) {
             predictDF <- stats::predict(object = modList[mod], # single bracket subset
                                         newdata = annTestDat[, -1],
                                         type = "response") %>%
               as.data.frame() %>%
               dplyr::select(prediction = 1) }
           
           else if(str_detect(names(modList[mod]), "xgboost")) {
             predictDF <- stats::predict(object = modList[[mod]],
                                         newdata = xgTestDat[, -1],
                                         type = "response") %>%
               as.data.frame() %>%
               dplyr::select(prediction = 1) }
           
           # Prediction vector
           prediction  <- as.vector(predictDF$prediction)
           
           # Defining the table to calculate metrics from        
           if(str_detect(names(modList[mod]), "glm|gam|mars|brt|rf|cf")) {
             
             m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                             truth = as.numeric(testDat$occurrence)) }
           
           else if(str_detect(names(modList[mod]), "max")) {
             
             m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                             truth = as.numeric(maxTestDat$occurrence))  }
           
           else if(str_detect(names(modList[mod]), "ann")) {
             
             m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                             truth = as.numeric(annTestDat$occurrence))  }
           
           else if(str_detect(names(modList[mod]), "xgboost")) {
             
             m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                             truth = as.numeric(xgTestDat[, 1]))  }
           
           list(prediction = prediction,
                
                metrics = data.frame(c(
                  TN = m %>% dplyr::filter(prediction == 0 & truth == 0) %>% count(),
                  FP = m %>% dplyr::filter(prediction == 1 & truth == 0) %>% count(),
                  FN = m %>% dplyr::filter(prediction == 0 & truth == 1) %>% count(),
                  TP = m %>% dplyr::filter(prediction == 1 & truth == 1) %>% count())) %>%
                  rename_with(~str_remove(., '.n')) %>%
                  mutate(sensitivity = TP/(TP+FN),
                         specificity = TN/(TN+FP),
                         threshold = threshold,
                         mae = Metrics::mae(actual = m$truth,
                                            predicted = prediction),
                         # note that AUC ROC is not reported in the manuscript
                         aucROC = mltools::auc_roc(preds = prediction,
                                                   actuals = m$truth)))
         } # close the internal function
  ) # close the lapply
} # close the overall function
buildPmodels_nocolremove <- function(i) {
  
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, occurrence, id))
  
  rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

  datList <- list(rowrm_train1 = rowrm_train1,
                  rowrm_train2 = rowrm_train2,
                  rowrm_train3 = rowrm_train3,
                  rowrm_train4 = rowrm_train4)
  
  lapply(datList, FUN = function(dat) {
    
    r <- ncol(dat)
    
    print(paste0("USING ", nrow(dat), " ROWS AND ", (ncol(dat)-1), " PREDICTORS"))
    
    set.seed(656)
    # BIOCLIM
    bioclim = dismo::bioclim(dat)
    
    set.seed(123)
    # DOMAIN
    domain = dismo::domain(dat)
    
    list(bioclim = bioclim, domain = domain) # close the list
  } # close the mini-function
  ) # close the lapply
} # close the overall function
evaluatePmodels_nocolremove <- function(i, modList, threshold) {
  
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  
  rowrm_test1 <- rowrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group)
  rowrm_test2 <- rowrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group)
  rowrm_test3 <- rowrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group)
  rowrm_test4 <- rowrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group)
  
  lapply(1:length(modList),
         FUN = function(mod) {
           
           # Setting which data to use to test each model
           if(str_detect(names(modList[mod]), "row") &
              str_detect(names(modList[mod]), "1")) { testDat <- rowrm_test1 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "2")) { testDat <- rowrm_test2 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "3")) { testDat <- rowrm_test3 }
           
           else if(str_detect(names(modList[mod]), "row") &
                   str_detect(names(modList[mod]), "4")) { testDat <- rowrm_test4 }
           
           set.seed(456)
           
           # Testing models
           prediction <- dismo::predict(object = modList[[mod]],
                                        x = testDat[, -1])
           
           # Defining the table to calculate metrics from 
           m <- data.frame(prediction = ifelse(prediction >= threshold, 1, 0),
                           truth = as.numeric(testDat$occurrence))
           
           list(prediction = prediction,
                
                metrics = data.frame(c(
                  TN = m %>% dplyr::filter(prediction == 0 & truth == 0) %>% count(),
                  FP = m %>% dplyr::filter(prediction == 1 & truth == 0) %>% count(),
                  FN = m %>% dplyr::filter(prediction == 0 & truth == 1) %>% count(),
                  TP = m %>% dplyr::filter(prediction == 1 & truth == 1) %>% count())) %>%
                  rename_with(~str_remove(., '.n')) %>%
                  mutate(sensitivity = TP/(TP+FN),
                         specificity = TN/(TN+FP),
                         threshold = threshold,
                         mae = Metrics::mae(actual = m$truth,
                                            predicted = prediction)))
         } # close the internal function
  ) # close the lapply
} # close the overall function

# 0.7 findCorrelation threshold
# Build models
noramPBmodels0 <- buildPAmodels_nocolremove(i = noramPB)

# Unnest the list for the evaluate functions
noramPBmodels <- unlist(noramPBmodels0, recursive = FALSE)

# If any models fail to converge or only partially build, exclude them
noramPBmodels[["rowrm_train1.ann"]] <- NULL
noramPBmodels[["rowrm_train2.ann"]] <- NULL
noramPBmodels[["rowrm_train3.ann"]] <- NULL
noramPBmodels[["rowrm_train4.ann"]] <- NULL

# Evaluate at threshold of 0.5
noramPBeval <- evaluatePAmodels_nocolremove(i = noramPB, modList = noramPBmodels, threshold = 0.5)

# Convert metrics to a data frame for saving and analyzing 
noramPBmetrics <- convertMetrics(eval = noramPBeval, modList = noramPBmodels)

# Write out metrics
write_delim(noramPBmetrics, "output/noramPBmetrics_alt.txt", delim = "\t", quote = "none")

# Write out models because they take so long to run
#saveRDS(noramPBmodels, "output/noramPBmodels.rds")
#saveRDS(noramPBeval, "output/noramPBeval.rds")

# North America, 0.5 findCorrelation threshold
# Build models
noramPB_Bmodels0 <- buildPAmodels_nocolremove(i = noramPB_B)

# Unnest the list for the evaluate functions
noramPB_Bmodels <- unlist(noramPB_Bmodels0, recursive = FALSE)

# Evaluate at threshold of 0.5
noramPB_Beval <- evaluatePAmodels(i = noramPB_B, modList = noramPB_Bmodels, threshold = 0.5)

# Convert metrics to a data frame for saving and analyzing 
noramPB_Bmetrics <- convertMetrics(eval = noramPB_Beval, modList = noramPB_Bmodels)

# Write out metrics
write_delim(noramPB_Bmetrics, "output/noramPB_Bmetrics_alt.txt", delim = "\t", quote = "none")

#saveRDS(noramPB_Bmodels, "output/noramPB_Bmodels.rds")
#saveRDS(noramPB_Beval, "output/noramPB_Beval.rds")

rm(list = ls(pattern = "noramPB_B"))

# 0.7 findCorrelation threshold, presence-only models
# Build models
noramPmodels0 <- buildPmodels_nocolremove(i = noramP)

# Unnest the list for the evaluate functions
noramPmodels <- unlist(noramPmodels0, recursive = FALSE)

# Evaluate at threshold of 0.5
noramPeval <- evaluatePmodels_nocolremove(i = noramP, modList = noramPmodels, threshold = 0.5)

# Convert metrics to a data frame for saving and analyzing 
noramPmetrics <- convertMetrics(eval = noramPeval, modList = noramPmodels)

# Write out metrics
write_delim(noramPmetrics, "output/noramPmetrics_alt.txt", delim = "\t", quote = "none")

#saveRDS(noramPmodels, "output/noramPmodels.rds")
#saveRDS(noramPeval, "output/noramPeval.rds")

# 0.5 findCorrelation threshold, presence-only models
# Build models
noramP_Bmodels0 <- buildPmodels_nocolremove(i = noramP_B)

# Unnest the list for the evaluate functions
noramP_Bmodels <- unlist(noramP_Bmodels0, recursive = FALSE)

# Evaluate at threshold of 0.5
noramP_Beval <- evaluatePmodels_nocolremove(i = noramP_B, modList = noramP_Bmodels, threshold = 0.5)

# Convert metrics to a data frame for saving and analyzing 
noramP_Bmetrics <- convertMetrics(eval = noramP_Beval, modList = noramP_Bmodels)

# Write out metrics
write_delim(noramP_Bmetrics, "output/noramP_Bmetrics_alt.txt", delim = "\t", quote = "none")

#saveRDS(noramP_Bmodels, "output/noramP_Bmodels.rds")
#saveRDS(noramP_Beval, "output/noramP_Beval.rds")

################################################################################

# Step 6: Analysis and Visualization
metricsList <- list.files(path = "output", pattern = "metrics_alt.txt", full.names = TRUE)
names(metricsList) <- tools::file_path_sans_ext(basename(metricsList))

metrics <- lapply(metricsList, read.delim, header = TRUE, sep = "\t",
                  fileEncoding = "UTF-8", row.names = NULL) %>%
  bind_rows(.id = "id") %>%
  rename(evalThreshold = threshold) %>%
  mutate(varRemovalThreshold = case_when(str_detect(id, "_B") ~ 0.5,
                                         TRUE ~ 0.7))
rm(metricsList)

# Controlling for the best method of removing missing values and best findCorrelation threshold
# in each set (just as we did in Step 6 for Questions 3, 4, and 5)

# Missing value removal method (removalType)
betterMAEgroup <- metrics %>%
  group_by(modelClass, modelType, scale, dataType, dataSource, removalType, varRemovalThreshold) %>%
  mutate(count = n()) %>%
  dplyr::filter(count >= 3) %>%
  dplyr::select(-count) %>%
  mutate(removalType = case_when(removalType == "row remove" ~ "rowrm",
                                 removalType == "col remove" ~ "colrm")) %>%
  summarize(mean = mean(mae)) %>%
  pivot_wider(names_from = removalType, values_from = mean) %>%
  mutate(betterMAE = case_when(rowrm < colrm | (is.na(colrm) & !is.na(rowrm)) ~ "row remove",
                               rowrm > colrm | (is.na(rowrm) & !is.na(colrm)) ~ "col remove")) %>%
  dplyr::select(-c(rowrm, colrm))

metrics2 <- metrics %>%
  right_join(betterMAEgroup, by = c("modelClass", "modelType", "scale", "dataType", "dataSource",
                                    "varRemovalThreshold", "removalType" = "betterMAE"))

rm(betterMAEgroup)

#findCorrelation threshold (varRemovalThreshold)
betterMAEgroup2 <-  metrics2 %>%
  group_by(modelClass, modelType, scale, dataType, dataSource, varRemovalThreshold) %>%
  mutate(count = n()) %>%
  #exclude any where one of the groups has fewer than 3
  dplyr::filter(count >= 3) %>%
  dplyr::select(-count) %>%
  summarize(mean = mean(mae)) %>%
  pivot_wider(names_from = varRemovalThreshold, values_from = mean) %>%
  mutate(betterMAE = case_when(`0.5` < `0.7` | (is.na(`0.7`) & !is.na(`0.5`)) ~ 0.5,
                               `0.5` > `0.7` | (is.na(`0.5`) & !is.na(`0.7`)) ~ 0.7)) %>%
  dplyr::select(-c(`0.5`, `0.7`))

metrics3 <- metrics2 %>%
  right_join(betterMAEgroup2, by = c("modelClass", "modelType", "scale", "dataType", "dataSource",
                                     "varRemovalThreshold" = "betterMAE"))

metrics3 %>%
  count(removalType, varRemovalThreshold)
# removing by row was always better, which makes sense given the small explanatory variable sets
# and the 0.7 threshold was usually better, which makes sense for the same reason

# This is an alternative to Question 3 (scale)
question3 <- bind_rows(
  metrics3 %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(mae ~ scale, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, scale) %>%
                summarize(meanMAE = mean(mae)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanMAE"),
              by = "modelType"),
  metrics3 %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(sensitivity ~ scale, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, scale) %>%
                summarize(meanSens = mean(sensitivity)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanSens"),
              by = "modelType"),
  metrics3 %>%
    dplyr::filter(modelClass != "Envelope") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(specificity ~ scale, paired = FALSE) %>%
    left_join(metrics3 %>%
                group_by(modelType, scale) %>%
                summarize(meanSpec = mean(specificity)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanSpec"),
              by = "modelType")) %>%
  mutate(sig = case_when(p <= 0.05 ~ "s",
                         p > 0.05 ~ "ns")) %>%
  dplyr::rename(metric = .y.) %>%
  mutate(metric = case_when(metric == "mae" ~ "Mean absolute error (MAE)",
                            metric == "sensitivity" ~ "Sensitivity (true positive rate)",
                            metric == "specificity" ~ "Specificity (true negative rate)"))

# Supplementary table for Question 3 (Appendix S2 in manuscript)
suppTable3 <- question3 %>%
  mutate(across(.cols = c(mich, noram),
                ~ round(., 4))) %>%
  mutate(sig = str_to_sentence(sig)) %>%
  dplyr::select("Model type" = modelType, "Metric" = metric, 
                "Group 1 n" = n1, "Group 2 n" = n2,
                "Wilcox statistic" = statistic, "P" = p, "Significance" = sig,
                "Mean value for Michigan" = mich, 
                "Mean value for North America" = noram)

#write_delim(suppTable3, "tables/suppTable3_alt.txt", delim = "\t", quote = "none")

# Visualize differences in each metric between scales
question3Long <- question3 %>%
  pivot_longer(cols = c(mich, noram),
               names_to = "scale",
               values_to = "meanValue") %>%
  mutate(meanValue = case_when(metric == "Mean absolute error (MAE)" ~ 1 - meanValue,
                               TRUE ~ meanValue),
         metric = case_when(metric == "Mean absolute error (MAE)" ~ "1 - mean absolute error (MAE)",
                            TRUE ~ metric))

# Plot (Appendix S1 in manuscript)
plot3 <- ggplot(question3Long, aes(x = meanValue, 
                                   y = reorder(modelType, desc(modelType)))) +
  geom_line() +
  geom_point(aes(color = scale, group = metric),
             size = 4) +
  geom_text(data = question3Long %>%
              mutate(sig = case_when(sig == "s" ~ "*",
                                     sig == "ns" ~ "")) %>%
              arrange(metric, modelType, meanValue) %>%
              distinct(modelType, metric, sig, .keep_all = TRUE),
            aes(label = sig), position = position_nudge(x = -0.085, y = -0.1),
            size = 6) +
  facet_wrap(vars(metric)) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
                     expand = c(0.04, 0.04)) +
  labs(x = "Mean value", y = "Model type",
       color = "Scale") +
  scale_color_manual(values = c("#37406B", "#96BBDA"),
                     labels = c("Michigan", "North America")) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color = "#000000", fill = NA,
                                    linewidth = 1),
        panel.grid.major.y = element_line(color = "#CDCBCF"),
        panel.spacing = unit(0.49, "cm"),
        strip.background = element_rect(color = "#000000", fill = NA,
                                        linewidth = 1),
        legend.position = "bottom",
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = ggplot2::margin(t = 0.25, r = 0, b = 0, l = 0, unit = "cm")),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0.25, b = 0, l = 0, unit = "cm")),
        axis.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.25, "in"),
        legend.key.width = unit(0.25, "in"),
        plot.margin = ggplot2::margin(t = 9, r = 9, b = 9, l = 9, unit = "pt"))

plot3
#ggsave("figures/plot3_alt.pdf", width = 7.25, height = 5, units = "in")
#ggsave("figures/plot3_alt.png", width = 7.25, height = 5, units = "in", dpi = 600)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Now we want to compare the mdoel outcomes, ie predicted probabilities
# across the Michigan study area

# Let's say 0.7 varRemovalThreshold, row remove, XGBoost is the best model for both
# which it seems to be
# Now we want to built full dataset models with those parameters to predict
# suitability for the Michigan scale (basically what I did at the end of Step 6)

# The question is: Does the performance differ between models built on only
# Michigan data versus the entire core invaded region of North America?
# It's a question of scale, so really an alternative Question 3 (Decision Point 3a in manuscript)

michPB_train <- read.delim("michPB.txt", header = TRUE,
                           sep = "\t", quote = "", fileEncoding = "UTF-8-BOM") %>%
  dplyr::select(occurrence, elev, nit, phos, pop, solar, wc18, wind) %>%
  drop_na() %>%
  as.matrix()

noramPB_train <- noramPB %>%
  drop_na() %>%
  dplyr::select(-c(x, y, id)) %>%
  as.matrix()

q1 <- ncol(michPB_train)

q2 <- ncol(noramPB_train)

buildTopModels <- function() {
  
  # Build models (keeping parameters the same as before but using entire dataset)
  # just like we did in Step 6 for Michigan and Saginaw Bay models
  
  # Michigan
  set.seed(512)
  michPB_XGBoost = xgboost::xgboost(data = michPB_train[, 2:q1],
                                    label = michPB_train[, 1],
                                    nrounds = 5000,
                                    max_depth = 10,
                                    eta = 0.8,
                                    objective = "binary:logistic",
                                    verbose = 0)
  
  michPB_varImportance = xgboost::xgb.importance(colnames(michPB_XGBoost),
                                                 model = michPB_XGBoost)
  
  michPB_XGBoost_n = nrow(michPB_train)
  
  
  # North America
  set.seed(512)
  noramPB_XGBoost = xgboost::xgboost(data = noramPB_train[, 2:q2],
                                     label = noramPB_train[, 1],
                                     nrounds = 5000,
                                     max_depth = 10,
                                     eta = 0.8,
                                     objective = "binary:logistic",
                                     verbose = 0)
  
  noramPB_varImportance = xgboost::xgb.importance(colnames(noramPB_XGBoost),
                                                  model = noramPB_XGBoost)
  
  noramPB_XGBoost_n = nrow(noramPB_train)
  
  
  list(michPB_XGBoost = michPB_XGBoost, michPB_varImportance = michPB_varImportance,
       michPB_XGBoost_n = michPB_XGBoost_n, noramPB_XGBoost = noramPB_XGBoost,
       noramPB_varImportance = noramPB_varImportance, noramPB_XGBoost_n = noramPB_XGBoost_n)
}

topModels <- buildTopModels()

# See how variables contributed
topModelsVarImportance <- topModels[c(2, 5)] %>%
  bind_rows(.id = "id") %>%
  mutate(Feature = case_when(Feature == "elev" ~ "Elevation",
                             Feature == "nit" ~ "Nitrogen",
                             Feature == "phos" ~ "Phosphorus application",
                             Feature == "pop" ~ "Population density",
                             Feature == "solar" ~ "Solar radiation",
                             Feature == "wc18" ~ "WorldClim BIO18",
                             Feature == "wind" ~ "Wind speed"))

topModelsVarImportance  %>%
  count(id, Feature, Gain) %>%
  arrange(desc(Gain))

ggplot(topModelsVarImportance, aes(x = Feature, y = Gain)) +
  geom_col(aes(fill = id), position = "dodge") +
  geom_label(aes(label = round(Gain*100, 0), group = id), 
             position = position_dodge(0.9))
# Note that we won't report these variable contributions, because they are 
# not the same variables used for the rest of the Decision Points
# Takeaway: Elevation and phosphorus contributed a little more for Michigan
# and climate variables solar, wc18, and wind contributed a little more for North America
# Nitrogen and populationd density were very close

# Now predict onto the full Michigan dataset with both models
fullMichScaled <- read.delim("fullMichScaled.txt", header = TRUE,
                             sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::select(x, y, id, elev, nit, phos, pop, solar, wc18, wind) %>%
  drop_na()

predValues <- function(i, model) {
  
  i2 <- i %>% dplyr::select(-c(x, y, id)) %>% as.matrix()
  
  set.seed(630)
  stats::predict(object = model,
                 newdata = i2,
                 type = "response") %>%
    bind_cols(i %>% drop_na() %>% dplyr::select(c(x, y))) %>%
    relocate(1, .after = y) %>%
    rename("pred" = 3)
  
}

michPB_predValues <- predValues(fullMichScaled, topModels[["michPB_XGBoost"]])
noramPB_predValues <- predValues(fullMichScaled, topModels[["noramPB_XGBoost"]])

summary(michPB_predValues$pred)
summary(noramPB_predValues$pred)

# Are predicted values significantly different?
# Paired Wilcox test
michPB_predValues %>%
  rename(michPred = pred) %>%
  inner_join(noramPB_predValues %>% rename(noramPred = pred),
             by = c("x", "y")) %>%
  pivot_longer(cols = c(michPred, noramPred),
               names_to = "scale", values_to = "pred") %>%
  rstatix::wilcox_test(pred ~ scale, paired = TRUE)

mean(noramPB_predValues$pred - michPB_predValues$pred)
# 0.05537 -> Michigan model predicts 5.5% lower on average

median(noramPB_predValues$pred - michPB_predValues$pred)
# median difference 0.000649 -> which is < 0.1% median difference

# Are predicted values different enough to change response?
# Wilcox effect size
michPB_predValues %>%
  rename(michPred = pred) %>%
  inner_join(noramPB_predValues %>% rename(noramPred = pred),
             by = c("x", "y")) %>%
  pivot_longer(cols = c(michPred, noramPred),
               names_to = "scale", values_to = "pred") %>%
  rstatix::wilcox_effsize(pred ~ scale, paired = TRUE)
# Effect size = 0.214 which is "small"

# Conclusion: While significant, the differences between predicted values are
# likely too small to really influence an application of these models

# Figures for manuscript
data("state_boundaries_wgs84")
baseMich <- subset(state_boundaries_wgs84, NAME == "Michigan" & TYPE == "Land",
                   select = Shape) %>% 
  terra::vect()
crs(baseMich) <- "+proj=longlat"

michBoundingBox <- terra::ext(c(-85, -82.2, 41.5, 46.5)) %>% terra::vect()
crs(michBoundingBox) <- crs(baseMich)

# Mapping functions (adjusted a little from Step 6)
predMap <- function(i, baseMap, boundingBox, title) {
  
  predXY <- i
  
  a = -85.0
  b = -82.5
  d = 42.0
  e = 46.0 # the Michigan map boundaries
  
  ggplot() +
    geom_sf(data = baseMap %>% terra::crop(boundingBox), fill = "#FFFFFF") +
    geom_sf(data = boundingBox, fill = NA, linewidth = 2/.pt) +
    geom_spatraster(data = predXY) +
    scale_x_continuous(breaks = c(a, b), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(d, e), expand = c(0, 0)) +
    scale_fill_gradient(low = "#DAE7F2", high = "#293051", na.value = NA) +
    labs(fill = "Probability of suitable habitat", 
         title = title) +
    theme_classic() +
    theme(axis.text = element_text(size = 8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 8.5),
          legend.position = "bottom",
          legend.title = element_text(size = 8.5),
          legend.text = element_text(size = 8))
}

comparisonMap <- function(i, baseMap, boundingBox, title) {
  diffXY <- i
  
  a = -85
  b = -82.5
  d = 42
  e = 46
  
  ggplot() +
    geom_sf(data = baseMap %>% terra::crop(boundingBox), fill = "#FFFFFF") +
    geom_sf(data = boundingBox, fill = NA, linewidth = 2/.pt) +
    geom_spatraster(data = diffXY) +
    scale_x_continuous(breaks = c(a, b), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(d, e), expand = c(0, 0)) +
    scale_fill_gradient2(low = "#f6999c", mid = "#FFFFFF",
                         high = "#28560b", na.value = NA) +
    labs(fill = "Difference in predicted probability of suitable habitat", 
         title = title) +
    theme_classic() +
    theme(axis.text = element_text(size = 8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 8.5),
          legend.position = "bottom",
          legend.title = element_text(size = 8.5),
          legend.text = element_text(size = 8))
}

# Prep data
michPB_predValues_rast <- michPB_predValues %>%
  terra::rast(type = "xyz", crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")

noramPB_predValues_rast <- noramPB_predValues %>%
  terra::rast(type = "xyz", crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")

# Prediction maps (Appendix S3 in manuscript)
michPBmap <- predMap(michPB_predValues_rast, 
                     baseMich, michBoundingBox,
                     "Prediction from Michigan points")

noramPBmap <- predMap(noramPB_predValues_rast,
                      baseMich, michBoundingBox,
                      "Prediction from North America points")

altPredMaps <- ggpubr::ggarrange(michPBmap, noramPBmap,
                                  nrow = 1, ncol = 2, common.legend = TRUE,
                                  legend = "bottom")
altPredMaps
#ggsave("figures/altPredMaps.pdf", width = 7.25, height = 8, units = "in")
#ggsave("figures/altPredMaps.png", width = 7.25, height = 8, units = "in", dpi = 600)

# Comparison map (Appendix S4 in manuscript)
comp_rast <- michPB_predValues %>% rename(michPred = pred) %>%
  inner_join(noramPB_predValues %>% rename(noramPred = pred),
             by = c("x", "y")) %>%
  mutate(diff = noramPred - michPred) %>%
  dplyr::select(-c(michPred, noramPred)) %>%
  terra::rast(type = "xyz", crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")

altCompMap <- comparisonMap(comp_rast, baseMich, michBoundingBox, 
                            "North America - Michigan")

altCompMap
#ggsave("figures/altCompMap.pdf", width = 7.25, height = 8, units = "in")
#ggsave("figures/altCompMap.png", width = 7.25, height = 8, units = "in", dpi = 600)

# noram models predict higher values in some pockets, particularly inland
# and mich models predict higher values in other pockets, particularly southeast area

# Depending on the goals of the model and where the area of interest is,
# you might make different decisions based on these two models
# But the more conservative model just depends on which area you're looking at

################################################################################

rm(list = ls())
