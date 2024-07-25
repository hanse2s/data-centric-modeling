# Data-centric SDM Step 4) Draw background points - final data processing step before creating models
## Sara Hansen
## Modified July 25, 2024

library(conflicted)
library(tidyverse)
library(terra)
library(dismo)
library(rio)

################################# READ IN DATA #################################

# Read in the European frog-bit datasets created in Step 3

efbMichP <- read.delim("output/efbMichP.txt", header = TRUE,
                       sep = "\t", quote = "", fileEncoding = "UTF-8")
efbMichSpecP <- read.delim("output/efbMichSpecP.txt", header = TRUE,
                           sep = "\t", quote = "", fileEncoding = "UTF-8")
efbMichObsP <- read.delim("output/efbMichObsP.txt", header = TRUE,
                          sep = "\t", quote = "", fileEncoding = "UTF-8")
efbSagPA <- read.delim("output/efbSagPA.txt", header = TRUE,
                       sep = "\t", quote = "", fileEncoding = "UTF-8")
efbSagP <- read.delim("output/efbSagP.txt", header = TRUE,
                      sep = "\t", quote = "", fileEncoding = "UTF-8")

# And raster stacks, at both correlation thresholds
# findCorrelation threshold 0.7
stackMich <- terra::rast("output/stackMich.tif")
stackSag <- terra::rast("output/stackSag.tif")

#findCorrelation threshold 0.5
stackMichB <- terra::rast("output/stackMichB.tif")

mich <- terra::subset(stackMich, 1) %>% raster::raster() # can be any of the layers
sag <- terra::subset(stackSag, 1) %>% raster::raster() # can be any of the layers

################################################################################

############################ DRAW BACKGROUND POINTS ############################

# We will draw background points for the four presence datasets
# The proportion of background to presence points should match as much as possible
# Use a simple 1:1 to reduce bias in either direction
# The only dataset that will violate this standard is the Saginaw Bay presence-absence
# data, which has a close to 1:1 ratio in the entire population but 
# fewer absences in the random sample taken in Step 3

# Draw random points using dismo::randomSample, which excludes presences
# and samples not more than 1 point per cell

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

# Run it on each presence dataset
efbMichPB <- drawBackground(efbMichP, mich)
efbMichSpecPB <- drawBackground(efbMichSpecP, mich)
efbMichObsPB <- drawBackground(efbMichObsP, mich)
efbSagPB <- drawBackground(efbSagP, sag)

# See that no two points are in the same cell (n is always 1):
#terra::cells(rast(sag), vect(efbSagPB, geom = c("x", "y"))) %>%
  #as.data.frame() %>% count(cell) %>% summarize(max(n))

################################################################################

############ EXTRACT VALUES OF EXPLANATORY VARIABLES FOR EACH POINT ############
# For all datasets -> presence-background, presence-only, and presence-absence
# At each correlation threshold

# Want land cover as character classes
landcoverClasses <- FedData::pal_nlcd() %>% dplyr::select(ID, Class)

extractEx <- function(i, stack) {
  bind_cols(i, terra::extract(stack, i[, 1:2])) %>%
    dplyr::select(-ID) %>%
    dplyr::rename(ID = layer) %>% # this is land cover
    left_join(landcoverClasses, by = "ID") %>%
    rename(land = Class) %>%
    dplyr::select(-ID)
}

# findCorrelation threshold 0.7
# Presence-background datasets
michPB0 <- extractEx(efbMichPB, stackMich)
michSpecPB0 <- extractEx(efbMichSpecPB, stackMich)
michObsPB0 <- extractEx(efbMichObsPB, stackMich)
sagPB0 <- extractEx(efbSagPB, stackSag)

# Presence-only datasets
michP0 <- extractEx(efbMichP, stackMich)
michSpecP0 <- extractEx(efbMichSpecP, stackMich)
michObsP0 <- extractEx(efbMichObsP, stackMich)
sagP0 <- extractEx(efbSagP, stackSag)

# Presence-absence
sagPA0 <- extractEx(efbSagPA, stackSag)

# findCorrelation threshold 0.5 (Michigan only)
# Presence-background datasets
michPB_B0 <- extractEx(efbMichPB, stackMichB)
michSpecPB_B0 <- extractEx(efbMichSpecPB, stackMichB)
michObsPB_B0 <- extractEx(efbMichObsPB, stackMichB)

# Presence-only datasets
michP_B0 <- extractEx(efbMichP, stackMichB)
michSpecP_B0 <- extractEx(efbMichSpecP, stackMichB)
michObsP_B0 <- extractEx(efbMichObsP, stackMichB)

# Saginaw Bay has only one environmental dataset because the 0.5 and 0.7 
# thresholds did not make a difference in removing correlated variables

rm(list = ls(pattern = "efb|ratio"))

################################################################################

###################### STANDARDIZE EXPLANATORY VARIABLES #######################
# By coercing all means to 0 and all standard deviations to 1

scaleVars <- function(i) {
  
  land <- i %>% dplyr::select(c(id, land))
  
  i2 <- i %>% dplyr::select(-land)
  
  r <- ncol(i2)
  
  i2 %>%
    dplyr::mutate(across(5:all_of(r),
                         ~ c(scale(.)))) %>% #c() makes output a data frame, not a series of matrices
    dplyr::left_join(land, by = "id")
}

# findCorrelation threshold 0.7
datList <- list(michPB0, michSpecPB0, michObsPB0, sagPB0,
                michP0, michSpecP0, michObsP0, sagP0, sagPA0)
datList2 <- lapply(datList, scaleVars)
names(datList2) <- c("michPB", "michSpecPB", "michObsPB", "sagPB", 
                     "michP", "michSpecP", "michObsP", "sagP", "sagPA")
list2env(datList2, envir = .GlobalEnv)

# Write out
write_delim(michPB, "output/michPB.txt", delim = "\t")
write_delim(michSpecPB, "output/michSpecPB.txt", delim = "\t")
write_delim(michObsPB, "output/michObsPB.txt", delim = "\t")
write_delim(sagPB, "output/sagPB.txt", delim = "\t")
write_delim(michP, "output/michP.txt", delim = "\t")
write_delim(michSpecP, "output/michSpecP.txt", delim = "\t")
write_delim(michObsP, "output/michObsP.txt", delim = "\t")
write_delim(sagP, "output/sagP.txt", delim = "\t")
write_delim(sagPA, "output/sagPA.txt", delim = "\t")

# findCorrelation threshold 0.5
datList_B <- list(michPB_B0, michSpecPB_B0, michObsPB_B0,
                  michP_B0, michSpecP_B0, michObsP_B0)
datList2_B <- lapply(datList_B, scaleVars)
names(datList2_B) <- c("michPB_B", "michSpecPB_B", "michObsPB_B",
                       "michP_B", "michSpecP_B", "michObsP_B")
list2env(datList2_B, envir = .GlobalEnv)

# Write out
write_delim(michPB_B, "output/michPB_B.txt", delim = "\t")
write_delim(michSpecPB_B, "output/michSpecPB_B.txt", delim = "\t")
write_delim(michObsPB_B, "output/michObsPB_B.txt", delim = "\t")
write_delim(michP_B, "output/michP_B.txt", delim = "\t")
write_delim(michSpecP_B, "output/michSpecP_B.txt", delim = "\t")
write_delim(michObsP_B, "output/michObsP_B.txt", delim = "\t")

################################################################################


############## PREPARE DATA FOR PREDICTIVE MAPPING IN LATER STEPS ##############

# For predictive mapping later on, we want variables across the entire study areas,
# not just at point locations

# Need mich and sag as SpatVectors
mich2 <- terra::subset(stackMich, 1) %>% as.polygons(values = FALSE) # can be any of the layers
sag2 <- terra::subset(stackSag, 1) %>% as.polygons(values = FALSE) # can be any of the layers

# Slightly modified scaleVars() function to work with these
scaleVars2 <- function(i) {
  
  i2 <- i %>% mutate(id = row_number()) %>% dplyr::select(-ID)
  
  land <- i2 %>% dplyr::select(c(id, land))
  
  i3 <- i2 %>% dplyr::select(-land)
  
  r <- ncol(i3) - 3 #for x and y and id
  
  i3 %>%
    dplyr::mutate(across(1:all_of(r),
                         ~ c(scale(.)))) %>%
    dplyr::left_join(land, by = "id")
}

# Michigan (first set, 0.7 correlation threshold)
fullMich <- terra::extract(stackMich, mich2, xy = TRUE) %>%
  rename(land = layer)
fullMichScaled <- scaleVars2(fullMich)
  
# Michigan (0.5 correlation threshold)
fullMichB <- terra::extract(stackMichB, mich2, xy = TRUE) %>%
  rename(land = layer)
fullMichBScaled <- scaleVars2(fullMichB)

# Saginaw Bay (0.7 correlation threshold)
fullSag <- terra::extract(stackSag, sag2, xy = TRUE) %>%
  rename(land = layer)
fullSagScaled <- scaleVars2(fullSag)

# Write out
write_delim(fullMichScaled, "output/fullMichScaled.txt", delim = "\t")
write_delim(fullMichBScaled, "output/fullMichBScaled.txt", delim = "\t")
write_delim(fullSagScaled, "output/fullSagScaled.txt", delim = "\t")

################################################################################

rm(list = ls())
