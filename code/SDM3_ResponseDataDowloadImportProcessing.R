# Data-centric SDM Step 3) Download, import, and process response data
## Sara Hansen
## Modified July 25, 2024

library(conflicted)
library(tidyverse)
library(geosphere)
library(terra)
library(tidyterra)
library(USA.state.boundaries)
library(ggspatial)
library(ggpubr)

########################### READ IN AND PROCESS DATA ###########################

# Import raster stacks
stackMich <- terra::rast("output/stackMich.tif")
stackSag <- terra::rast("output/stackSag.tif")

# Import European frog-bit data
# Saginaw Bay
efbSag0 <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                        sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::filter(bibliographicCitation == 
                  "Monfils, A. (2020). [European frog-bit occurrences in Saginaw Bay]. Unpublished data. Central Michigan University Herbarium." &
                  decimalLatitude > 43) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  mutate(occurrence = recode(occurrenceStatus,
                             "PRESENT" = as.integer(1),
                             "ABSENT" = as.integer(0)),
         id = row_number()) %>%
  dplyr::select(occurrenceID, x = decimalLongitude, y = decimalLatitude, 
                occurrence, id)

# Michigan
h <- 1000 # threshold (in meters) for excluding records

efbMichP0 <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                         sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::filter(occurrenceStatus == "PRESENT" & stateProvince == "Michigan" &
                  !is.na(decimalLongitude) & !is.na(decimalLatitude) &
                  (is.na(coordinateUncertaintyInMeters)
                   | coordinateUncertaintyInMeters <= h) &
                  decimalLongitude > -85.5) %>%
  anti_join(efbSag0, by = "occurrenceID") %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  mutate(occurrence = as.integer(1)) %>%
  dplyr::select(basisOfRecord, occurrenceID, x = decimalLongitude,
                y = decimalLatitude, occurrence)

efbSag <- efbSag0 %>%
  dplyr::select(-occurrenceID)
rm(efbSag0)

# Michigan specimens only
efbMichSpecP0 <- efbMichP0 %>%
  dplyr::filter(basisOfRecord == "PRESERVED_SPECIMEN")

# Michigan observations only
efbMichObsP0 <- efbMichP0 %>%
  dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION")

# In order to avoid counting occurrences from the same immediate area multiple times,
# we will use dismo::gridSample() to randomly select one occurrence per cell

set.seed(345)
efbMichP_gs <- dismo::gridSample(efbMichP0 %>% 
                                   dplyr::select(x, y),
                                 r = raster::raster(stackMich), n = 1) %>% # n = 461
  mutate(occurrence = as.integer(1),
         id = row_number())

set.seed(345)
efbMichSpecP_gs <- dismo::gridSample(efbMichSpecP0 %>%
                                       dplyr::select(x, y),
                                     r = raster::raster(stackMich), n = 1) %>% # n = 34
  mutate(occurrence = as.integer(1),
         id = row_number())

set.seed(345)
efbMichObsP_gs <- dismo::gridSample(efbMichObsP0 %>%
                                      dplyr::select(x, y),
                                    r = raster::raster(stackMich), n = 1) %>% # n = 456
  mutate(occurrence = as.integer(1),
         id = row_number())

# Saginaw Bay presence-absence
# Crop Saginaw Bay data to the extent established to account for edges
efbSagPA <- bind_cols(efbSag, terra::extract(stackSag, efbSag[, 1:2])) %>%
  dplyr::filter(!is.na(boat) | !is.na(elev) | !is.na(nit) | 
                  !is.na(pop) | !is.na(wc10) | !is.na(layer)) %>%
  dplyr::select(x, y, occurrence, id)
# These points are not an issue, they are just missing data because of the way
# the study area was cropped

# Plot:
#sag <- sf::read_sf("sag.shp")
#ggplot() +
  #geom_sf(data = sag) +
  #geom_spatraster(data = stackSag[[1]]) +
  #geom_point(data = efbSag, aes(x = x, y = y))

# See that the points with missing data fall outside the study area
#ggplot() +
  #geom_sf(data = sag) +
  #geom_spatraster(data = stackSag[[1]]) +
  #geom_point(data = efbSag %>% anti_join(efbSagPA), aes(x = x, y = y))

#rm(sag)

set.seed(543)
efbSagPA_gs0 <- dismo::gridSample(efbSagPA %>%
                                   dplyr::select(x, y),
                                 r = raster::raster(stackSag), n = 1) # 115 observations
# And then just join this to the original data to get the identifiers and
# see which ones are presence and absence
efbSagPA_gs0.2 <- efbSagPA_gs0 %>%
  left_join(efbSag, by = c("x", "y"))

efbSagPA_gs0.2 %>% count(occurrence) # close to even: 50 absences, 65 presences
# To make this more comparable to other datasets which will use a 1:1 ratio,
# randomly same presences equal to absences

efbSagPA_gs <- efbSagPA_gs0.2 %>%
  dplyr::filter(occurrence == 1) %>%
  sample_frac(1) %>% # shuffle
  head(50) %>% # and take top 50
  rbind(efbSagPA_gs0.2 %>% dplyr::filter(occurrence == 0)) %>%
  sample_frac(1)
# Sample size is 100, which matches the Saginaw Bay presence data size entirely by chance

plot(stackSag[[1]])
plot(vect(efbSagPA_gs, geom = c("x", "y")), add = TRUE)

# Saginaw Bay presence
efbSagP <- efbSagPA %>%
  dplyr::filter(occurrence == 1)

set.seed(345)
efbSagP_gs <- dismo::gridSample(efbSagP %>%
                                  dplyr::select(x, y),
                                r = raster::raster(stackSag), n = 1) %>% #100 observations
  left_join(efbSagP, by = c("x", "y"))

################################################################################

############################# PLOT FOR MANUSCRIPT ##############################
mich <- sf::read_sf("output/mich.shp")

sag <- sf::read_sf("output/sag.shp")

data("state_boundaries_wgs84")
michLand0 <- subset(state_boundaries_wgs84, NAME == "Michigan" & TYPE == "Land",
                    select = Shape) %>%
  terra::vect()
crs(michLand0) <- "+proj=longlat"

boundingBox <- terra::ext(c(-85, -82.2, 41.5, 46.5)) %>% 
  terra::vect()
crs(boundingBox) <- crs(mich)

michPlot <- ggplot() +
  geom_sf(data = michLand0 %>% terra::crop(boundingBox), fill = "#CDCBCF") +
  geom_sf(data = boundingBox, fill = NA, linewidth = 2/.pt) +
  geom_sf(data = mich, fill = "#FFFFFF", linewidth = 1/.pt) +
  geom_point(data = efbMichP0, aes(x = x, y = y), color = "#28560b", size = 0.5) +
  scale_x_continuous(breaks = c(-85, -82.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(42, 46), expand = c(0, 0)) +
  annotation_scale(location = "tr", pad_y = unit(0.4, "cm"), pad_x= unit(0.3, "cm")) +
  annotation_north_arrow(location = "tr", pad_y = unit(0.8, "cm"), pad_x = unit(0.2, "cm"),
                         height = unit(0.8, "cm"), width = unit(0.6, "cm"),
                         style=north_arrow_orienteering(text_size = 6)) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_blank())

michPlotSmall <- ggplot() +
  geom_sf(data = michLand0, fill = "#CDCBCF") +
  geom_sf(data = boundingBox, fill = NA, linewidth = 2/.pt) +
  theme_void()

boundingBox2 <- terra::ext(c(-84, -83.35, 43.55, 44.05)) %>% terra::vect()
crs(boundingBox2) <- crs(sag)

sagPlot <- ggplot() +
  geom_sf(data = michLand0 %>% terra::crop(boundingBox2), fill = "#CDCBCF") +
  geom_sf(data = boundingBox2, fill = NA, linewidth = 2/.pt) +
  geom_sf(data = sag, fill = "#FFFFFF", linewidth = 1/.pt) +
  geom_point(data = efbSagPA, aes(x = x, y = y), color = "#28560b", size = 0.5) +
  scale_x_continuous(breaks = c(-83.9, -83.4), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(43.6, 44.0), expand = c(0, 0)) +
  annotation_scale(location = "tr", pad_y = unit(0.3, "cm"), pad_x= unit(0.3, "cm")) +
  annotation_north_arrow(location = "tr", pad_y = unit(0.7, "cm"), pad_x = unit(0.3, "cm"),
                         height = unit(0.8, "cm"), width = unit(0.6, "cm"),
                         style=north_arrow_orienteering(text_size = 6)) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_blank())

michPlotSmall2 <- ggplot() +
  geom_sf(data = michLand0, fill = "#CDCBCF") +
  geom_sf(data = boundingBox2, fill = NA, linewidth = 2/.pt) +
  theme_void()

# Combine maps (Appendix S5 in manuscript)
studyAreasPlot <- ggarrange(michPlotSmall, michPlotSmall2, michPlot, sagPlot, 
                            ncol = 2, nrow = 2, labels = c("A", "B", "", ""),
                            heights = c(1, 2))

studyAreasPlot
#ggsave("figures/studyAreasPlot.pdf", width = 7.25, height = 7.25, units = "in")
#ggsave("figures/studyAreasPlot.png", width = 7.25, height = 7.25, units = "in", dpi = 600)

################################################################################

################################################################################

# The environment should have five datasets labeled "_gs"
# Write each one out to a csv file

write_delim(efbMichP_gs, "output/efbMichP.txt", delim = "\t")
write_delim(efbMichSpecP_gs, "output/efbMichSpecP.txt", delim = "\t")
write_delim(efbMichObsP_gs, "output/efbMichObsP.txt", delim = "\t")
write_delim(efbSagPA_gs, "output/efbSagPA.txt", delim = "\t")
write_delim(efbSagP_gs, "output/efbSagP.txt", delim = "\t")

rm(list = ls())
