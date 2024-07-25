# Data-centric SDM Step 1) Download explanatory data from online sources
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

########################### DOWNLOAD RESPONSE DATA #############################

# Michigan occurrences were downloaded from the Global Biodiversity Information Facility:
# https://www.gbif.org/dataset/71454d8a-6e9c-49f5-bf37-353f9ad2e2b9
# For more information, see:
# Hansen, S. E., B. C. Cahill, R. A. Hackett, M. J. Monfils, R. T. Goebel, S. Asencio, and A. Monfils. 2022. 
# Aggregated occurrence records of invasive European frog-bit (Hydrocharis morsus-ranae L.) across North America.
# Biodiversity Data Journal 10: e77492. https://doi.org/10.3897/BDJ.10.e77492

# The downloaded data set is available in the data folder.
# Cite the data paper if you use the data set.

# Saginaw Bay occurrences will be extracted from this same dataset

################################################################################

############################## DEFINE STUDY AREAS ##############################

# Two study areas are being assessed: Michigan and Saginaw Bay
# To define each study area, we need European frog-bit occurrence points

# Study area = Saginaw Bay, Michigan
# We will use a subset of data collected during a field study in 2020
sagPoints <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                        sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::filter(bibliographicCitation == 
                  "Monfils, A. (2020). [European frog-bit occurrences in Saginaw Bay]. Unpublished data. Central Michigan University Herbarium." &
                  decimalLatitude > 43) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

sagPoints2 <- terra::vect(sagPoints %>% 
                            dplyr::select(x = decimalLongitude, y = decimalLatitude),
                          geom = c("x", "y"))
crs(sagPoints2) <- "+proj=longlat"

# Define a concave hull around points to establish the sampled area
sagHull <- concaveman::concaveman(sagPoints2 %>% 
                                    sf::st_as_sf(),
                                  concavity = 1) %>%
  terra::vect() # concavity = 1 follows the natural sampling path

# Define study area as 1 km buffer around concave hull
# then clipped to 1 km around land so most of the study area is not in the open water
# Michigan land area
data("state_boundaries_wgs84")
michLand0 <- subset(state_boundaries_wgs84, NAME == "Michigan" & TYPE == "Land",
                    select = Shape) %>%
  terra::vect()
crs(michLand0) <- "+proj=longlat"

# Create 1km buffer around land
michLand <- terra::buffer(michLand0, width = 1000) #unit is meters
crs(michLand) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

# Create study area
sag <- terra::buffer(sagHull, width = 1000) %>%
  terra::crop(michLand)
crs(sag) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

# Note that points further into the open water than 1km will later be excluded
# because they are not within the defined study area

# Visualize:
plot(sag)
plot(sagPoints2, add = TRUE)

# One more sagPoints layer, excluding the ones that fell outside the study area
# because of cropping to the land (this is not the layer for modeling yet, just mapping)
# Note this is not the modeling data yet, just mapping
# So points may not perfectly represent the actual points model
sagPoints3 <- sagPoints2 %>%
  terra::crop(sag)

# Study area = Eastern Michigan
# The spatial scale of explanatory variables will be 30 arc seconds, or about 1 km
# so we will limit our uncertainty to that amount
# We will also exclude the points used in the Saginaw Bay study area

h <- 1000 # threshold (in meters) for excluding records

michPoints <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                          sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::filter(occurrenceStatus == "PRESENT" & stateProvince == "Michigan" &
                  !is.na(decimalLongitude) & !is.na(decimalLatitude) &
                  (is.na(coordinateUncertaintyInMeters)
                   | coordinateUncertaintyInMeters <= h) &
                  decimalLongitude > -85.5) %>%
  anti_join(sagPoints) %>%
  distinct(decimalLongitude, decimalLatitude) %>%
  dplyr::select(x = decimalLongitude, y = decimalLatitude)

michPoints2 <- terra::vect(michPoints, geom = c("x", "y"))
crs(michPoints2) <- "+proj=longlat"

# Define a concave hull around points to establish the sampled area
michHull <- concaveman::concaveman(michPoints2 %>% 
                                     sf::st_as_sf()) %>%
  terra::vect()

# Define study area as 1 km buffer around concave hull, clipped to 1 km around land
# already have "michLand" which is a 1km buffer around Michigan land area

# Create study area
mich <- terra::buffer(michHull, width = 1000) %>%
  terra::crop(michLand)
crs(mich) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

################################################################################

# Write out the two extents for future use
terra::writeVector(sag, "output/sag.shp", filetype = "ESRI Shapefile")
terra::writeVector(mich, "output/mich.shp", filetype = "ESRI Shapefile")

rm(list = ls(pattern = "Points|Hull|Land|Box|Plot|state"))

######################## DOWNLOAD EXPLANATORY VARIABLES ########################
geodata_path("output")

# WorldClim bioclimatic variables
wc0 <- geodata::worldclim_global(var = "bio", res = 0.5, version = "2.1")

# WorldClim solar radiation
solar0 <- geodata::worldclim_global(var = "srad", res = 0.5, version = "2.1")

# WorldClim wind speed
wind_speed0 <- geodata::worldclim_global(var = "wind", res = 0.5, version = "2.1")

# WorldClim water vapor pressure
vapor0 <- geodata::worldclim_global(var = "vapr", res = 0.5, version = "2.1") 

# population density in 2020
pop0 <- geodata::population(year = 2020, res = 0.5) 

# elevation
elev0 <- geodata::elevation_global(res = 0.5) 

# nitrogen at surface of soil (0-5 cm)
nitrogen0 <- geodata::soil_world(var = "nitrogen", depth = 5)

# phosphorus fertilizer application from SEDAC
phosphorus0 <- download.file(url = "https://daac.ornl.gov/bundle/Global_Phosphorus_Dist_Map_1223.zip",
                             destfile = "output/phosphorus0.zip")

unzip("output/phosphorus0.zip", exdir = "output/phosphorus0")

# land cover (for each study area individually)
landcoverM0 <- FedData::get_nlcd(template = mich, label = "nlcd",
                                year = 2019, dataset = "landcover",
                                extraction.dir = "output/michLandcover")

landcoverS0 <- FedData::get_nlcd(template = sag, label = "nlcd",
                                 year = 2019, dataset = "landcover",
                                 extraction.dir = "output/sagLandcover")
# note that both landcover datasets will be resampled later

# Michigan boat access sites
boat_launch0 <- fromJSON("https://services3.arcgis.com/Jdnp1TjADvSDxMAX/ArcGIS/rest/services/dnrParksAndRecreation/FeatureServer/1/query?where=1%3D1&outFields=*&outSR=4326&f=json") %>%
  purrr::pluck("features") %>%
  tidyr::unnest(cols = c(attributes, geometry)) %>%
  dplyr::filter(!is.na(long) & !is.na(lat))

boat_launch_spatial_points <- sf::st_as_sf(boat_launch0, coords = c("long", "lat"),
                                            crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")

# distanceFromPoints() requires Raster objects
emptyRaster <- terra::rast(ext(wc0), resolution = res(wc0))
crs(emptyRaster) <- crs(wc0)

michRaster <- terra::rasterize(mich, emptyRaster) %>%
  terra::crop(mich, mask = T)

sagRaster <- terra::rasterize(sag, emptyRaster) %>%
  terra::crop(sag, mask = T)

plot(michRaster)
plot(sagRaster)

sagRaster

# boat launches in Michigan
boat_launch_distM0 <- distanceFromPoints(object = raster::raster(michRaster),
                                         xy = boat_launch_spatial_points)

terra::writeRaster(boat_launch_distM0, "output/boat_launch_distM0.tif",
                   filetype = "GTiff")

# boat launches in Saginaw Bay
boat_launch_distS0 <- distanceFromPoints(object = raster::raster(sagRaster),
                                         xy = boat_launch_spatial_points)

terra::writeRaster(boat_launch_distS0, "output/boat_launch_distS0.tif",
                   filetype = "GTiff")

################################################################################

# Clear workspace
rm(list = ls())
