# Data-centric SDM Step 7) Data exploration and visualization for manuscript
## Sara Hansen
## Modified July 25, 2024

library(conflicted)
library(tidyverse)
library(terra)
library(sf)
library(USA.state.boundaries)
library(tidyterra)

# 0.7 findCorrelation threshold data
michPB <- read.delim("output/michPB.txt", header = TRUE,
                       sep = "\t", quote = "", fileEncoding = "UTF-8")
michSpecPB <- read.delim("output/michSpecPB.txt", header = TRUE,
                         sep = "\t", quote = "", fileEncoding = "UTF-8")
michObsPB <- read.delim("output/michObsPB.txt", header = TRUE,
                        sep = "\t", quote = "", fileEncoding = "UTF-8")
sagPB <- read.delim("output/sagPB.txt", header = TRUE,
                           sep = "\t", quote = "", fileEncoding = "UTF-8")
sagPA <- read.delim("output/sagPA.txt", header = TRUE,
                          sep = "\t", quote = "", fileEncoding = "UTF-8")

# 0.5 findCorrelation threshold data
michPB_B <- read.delim("output/michPB_B.txt", header = TRUE,
                     sep = "\t", quote = "", fileEncoding = "UTF-8")
michSpecPB_B <- read.delim("output/michSpecPB_B.txt", header = TRUE,
                         sep = "\t", quote = "", fileEncoding = "UTF-8")
michObsPB_B <- read.delim("output/michObsPB_B.txt", header = TRUE,
                        sep = "\t", quote = "", fileEncoding = "UTF-8")


# Data for mapping:
# Saginaw Bay (just for the anti_join on the Michigan EFB layer)
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

efbMich <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
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
                y = decimalLatitude, occurrence, recordedBy, year)

rm(efbSag0)

# Other layers for plots
data("state_boundaries_wgs84")
baseMich <- subset(state_boundaries_wgs84, NAME == "Michigan" & TYPE == "Land",
                   select = Shape) %>%
  terra::vect()
crs(baseMich) <- "+proj=longlat"

michBoundingBox <- terra::ext(c(-85, -82.2, 41.5, 46.5)) %>% terra::vect()
crs(michBoundingBox) <- crs(baseMich)

# Saginaw Bay
sagBoundingBox <- terra::ext(c(-84, -83.35, 43.55, 44.05)) %>% terra::vect()
crs(sagBoundingBox) <- crs(baseMich)

# Mapping function
pointsMap <- function(i, absenceType = NA, pointSize, baseMap, boundingBox) {
  
  if(str_detect(deparse(substitute(boundingBox)), "mich")) {
    a = -85
    b = -82.5
    d = 42
    e = 46 }
  
  else if(str_detect(deparse(substitute(boundingBox)), "sag")) {
    a = -83.9
    b = -83.4
    d = 43.6
    e = 44.0 }

ggplot() +
  geom_sf(data = baseMap %>% terra::crop(boundingBox), fill = "#FFFFFF") +
  geom_sf(data = boundingBox, fill = NA, linewidth = 2/.pt) +
  geom_point(data = i %>%
               mutate(occurrence = case_when(occurrence == "0" ~ as.character(absenceType),
                                             occurrence == "1" ~ "presence")) %>%
               mutate(occurrence = factor(occurrence, levels = c("presence", as.character(absenceType)))),
             aes(x = x, y = y, color = occurrence),
             size = pointSize) +
  scale_x_continuous(breaks = c(a, b), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(d, e), expand = c(0, 0)) +
  scale_color_manual(values = c("#28560b", "#f6999c")) +
  labs(color = "Occurrence status") +
  theme_classic() +
  theme(axis.text = element_text(size = 9),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) }

#################### STATISTICS AND FIGURES FOR MANUSCRIPT #####################

# Visualize distribution of absences vs. backgrounds point in Saginaw Bay
sagPAmap <- pointsMap(sagPA, absenceType = "absence", pointSize = 1, 
                      baseMap = baseMich, boundingBox = sagBoundingBox)
sagPAmap

sagPBmap <- pointsMap(sagPB, absenceType = "background", pointSize = 1, 
                      baseMap = baseMich, boundingBox = sagBoundingBox)

sagPBmap
# Greater distribution of background points compared with absences

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Visualize distribution of specimens and observations in Michigan
michSpecPmap <- pointsMap(efbMich %>% 
                            dplyr::filter(basisOfRecord == "PRESERVED_SPECIMEN"),
                          pointSize = 1, baseMap = baseMich,
                          boundingBox = michBoundingBox)

michObsPmap <- pointsMap(efbMich %>% 
                           dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION"),
                         pointSize = 1, baseMap = baseMich,
                         boundingBox = michBoundingBox)

# Combine specimen and observation plots (Appendix S6 in manuscript)
dataSourceMaps <- ggpubr::ggarrange(michSpecPmap + ggtitle("Specimens"), 
                                    michObsPmap + ggtitle("Observations"),
                                    ncol = 2, legend = "none")
dataSourceMaps
#ggsave("figures/dataSourceMaps.pdf", width = 7.25, height = 7.25, units = "in")
#ggsave("figures/dataSourceMaps.png", width = 7.25, height = 7.25, units = "in", dpi = 600)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Who observed specimens in each year?
efbMich %>%
  dplyr::filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
  count(recordedBy, year)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# What is the size of each study area? (square kilometers)
stackMich <- terra::rast("output/stackMich.tif")
mich <- terra::subset(stackMich, 1) %>% as.polygons(values = FALSE)
sum(expanse(mich, unit = "km"))

stackSag <- terra::rast("output/stackSag.tif")
sag <- terra::subset(stackSag, 1) %>% as.polygons(values = FALSE)
sum(expanse(sag, unit = "km"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Where are observations with missing values for the large Michigan datasets
# Use the pointsMap function, but understand that "occurrence" actually just means
# whether there are any missing values for explanatory variables
# "presence" = nothing missing, "missing" = something is missing
michPB_missing <- michPB %>%
  dplyr::select(-occurrence) %>%
  drop_na() %>%
  mutate(occurrence = 1) %>%
  bind_rows(michPB %>%
              anti_join(michPB %>%
                          dplyr::select(-occurrence) %>%
                          drop_na(),
                        by = c("x", "y")) %>%
              mutate(occurrence = 0))

pointsMap(michPB_missing, absenceType = "missing", pointSize = 1,
          baseMap = baseMich, boundingBox = michBoundingBox)

# Zoom in on specific parts:
northBoundingBox <- terra::ext(c(-84, -83, 43, 46.5)) %>% terra::vect()
crs(northBoundingBox) <- crs(baseMich)

# northern half:
ggplot() +
  geom_sf(data = baseMich %>% terra::crop(northBoundingBox), fill = "#FFFFFF") +
  geom_sf(data = northBoundingBox, fill = NA, linewidth = 2/.pt) +
  geom_point(data = michPB_missing %>%
               mutate(occurrence = case_when(occurrence == "0" ~ "missing",
                                             occurrence == "1" ~ "not missing")) %>%
               dplyr::filter(x > -84 & y > 43),
             aes(x = x, y = y, color = occurrence), size = 1) +
  scale_color_manual(values = c("#28560b", "#f6999c")) +
  labs(color = "Missing status") +
  theme_classic() +
  theme(axis.text = element_text(size = 9),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) 

# southeastern area
southBoundingBox <- terra::ext(c(-84, -82.5, 41.6, 43)) %>% terra::vect()
crs(southBoundingBox) <- crs(baseMich)

ggplot() +
  geom_sf(data = baseMich %>% terra::crop(southBoundingBox), fill = "#FFFFFF") +
  geom_sf(data = southBoundingBox, fill = NA, linewidth = 2/.pt) +
  geom_point(data = michPB_missing %>%
               mutate(occurrence = case_when(occurrence == "0" ~ "missing",
                                             occurrence == "1" ~ "not missing")) %>%
               dplyr::filter(x > -84 & y < 43),
             aes(x = x, y = y, color = occurrence), size = 1) +
  scale_color_manual(values = c("#28560b", "#f6999c")) +
  labs(color = "Missing status") +
  theme_classic() +
  theme(axis.text = element_text(size = 9),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) 

# Saginaw Bay area only (kind of central part of Michigan study area)
ggplot() +
  geom_sf(data = baseMich %>% terra::crop(sagBoundingBox), fill = "#FFFFFF") +
  geom_sf(data = sagBoundingBox, fill = NA, linewidth = 2/.pt) +
  geom_point(data = michPB_missing %>%
               mutate(occurrence = case_when(occurrence == "0" ~ "missing",
                                             occurrence == "1" ~ "not missing")) %>%
               dplyr::filter(x > -84 & x < -83.35 & y > 43.55 & y < 44.05),
             aes(x = x, y = y, color = occurrence), size = 1) +
  scale_color_manual(values = c("#28560b", "#f6999c")) +
  labs(color = "Missing status") +
  theme_classic() +
  theme(axis.text = element_text(size = 9),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) 

# Not seeing a lot of evidence that missing values are necessarily restricted to
# any specific area or type of landscape

# What about missing values relative to land cover?
# (land cover is not usually missing)

michPB_missing %>%
  count(occurrence)

michPB_missing %>%
  count(occurrence, land) %>%
  mutate(prop = case_when(occurrence == 0 ~ n / 139,
                          occurrence == 1 ~ n / 783))

# Larger proportion of points with missing values fall in Open Water land class

michObsPB_missing <- michObsPB %>%
  dplyr::select(-occurrence) %>%
  drop_na() %>%
  mutate(occurrence = 1) %>%
  bind_rows(michObsPB %>%
              anti_join(michObsPB %>%
                          dplyr::select(-occurrence) %>%
                          drop_na(),
                        by = c("x", "y")) %>%
              mutate(occurrence = 0))

michObsPB_missing %>%
  count(occurrence)

michObsPB_missing %>%
  count(occurrence, land) %>%
  mutate(prop = case_when(occurrence == 0 ~ n / 137,
                          occurrence == 1 ~ n / 775))

# The same is true for Michigan observations

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Which variables were lost with the stricter findCorrelation threshold?
setdiff(colnames(michPB), colnames(michPB_B))
# solar, wc18, wind
# None of which contributed a lot in the final XGBoost models for Michigan all and observations data
# which used the 0.7 threshold
# So it makes sense that the threshold performed similarly

# Note that the final XGBoost models for Michigan all data and observation data
# used the row removal method, so the variables included in the data were all modeled
# However, specimen data used the column removal method

summary(michSpecPB)
summary(michSpecPB_B)
# Only nitrogen has missing values though, so those variables wouldn't be affected

################################################################################

rm(list = ls())

