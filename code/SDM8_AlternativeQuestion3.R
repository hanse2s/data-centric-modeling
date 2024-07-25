# Data-centric SDM Step 8) Alternative options for assessing Question 3b: effect of scale
## Sara Hansen
## Modified July 25, 2023

library(conflicted)
library(tidyverse)
library(terra)
library(sf)
library(tidyterra)
library(mgcv)
library(earth)
library(dismo)
library(rJava)
library(FedData)
library(gbm) 
library(randomForest)
library(party)
library(maxnet)
library(neuralnet)
library(xgboost)
library(caret)
library(Metrics)

################################################################################

# The Michigan study area is both larger and has a larger sample size
# than the Saginaw Bay study area -> is one of these factors causing differences?

# This alternative code will provide options for redefining the Michigan scale
# to analyze model performance for Question 3. The two ways of redefining the scale are:
# 1) the same area as Saginaw Bay ("EqA")
# 2) the same sample size as Saginaw Bay ("EqN")

# Note that the correlation thresholds (0.5 and 0.7) results in the same dataset
# for Saginaw Bay, so this automatically controlled for (Question 2)
# Data will be from all sources, presence-background (controlling for Questions 4 and 5)
# Question 1 (method of removing missing values) will be controlled for by selecting better one

################################################################################

############# 1) Michigan scale = same area as Saginaw Bay ("EqA") #############

# Data preparation

# Define study area
stackSag <- terra::rast("output/stackSag.tif") #also need this for explanatory variables
sag <- terra::subset(stackSag, 1) %>% as.polygons(values = FALSE)

michEqA <- sag
michEqA_rast <- terra::subset(stackSag, 1) %>% raster::raster() # can be any of the layers

# Extract the Michigan presence clusters that fall within the new area

sagPoints <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                        sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::filter(bibliographicCitation == 
                  "Monfils, A. (2020). [European frog-bit occurrences in Saginaw Bay]. Unpublished data. Central Michigan University Herbarium." &
                  decimalLatitude > 43) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

h <- 1000 # threshold (in meters) for excluding records

efbMichP0 <- read.delim("data/MichiganEFBOccurrence.txt", header = TRUE,
                         sep = "\t", quote = "", fileEncoding = "UTF-8") %>%
  dplyr::filter(occurrenceStatus == "PRESENT" & stateProvince == "Michigan" &
                  !is.na(decimalLongitude) & !is.na(decimalLatitude) &
                  (is.na(coordinateUncertaintyInMeters)
                   | coordinateUncertaintyInMeters <= h) &
                  decimalLongitude > -85.5) %>%
  anti_join(sagPoints) %>%
  distinct(decimalLongitude, decimalLatitude) %>%
  dplyr::select(x = decimalLongitude, y = decimalLatitude) %>%
  terra::vect(geom = c("x", "y"),
              crs = "+proj=longlat +ellps=WGS84 +datum=WGS84") %>%
  terra::crop(michEqA) %>%
  as.data.frame(geom = "XY")

ggplot() +
  geom_sf(data = michEqA) +
  geom_point(data = efbMichP, aes(x = x, y = y))

# Remove any points that will have only missing values because of where they fall
# in the grid cells relative to the study area
efbMichP0.2 <- bind_cols(efbMichP0, terra::extract(stackSag, efbMichP0[, 1:2])) %>%
  dplyr::select(-ID) %>%
  dplyr::filter(!is.na(boat) | !is.na(elev) | !is.na(nit) | 
                  !is.na(pop) | !is.na(wc10) | !is.na(layer)) %>%
  dplyr::select(x, y)

count(efbMichP0) == count(efbMichP0.2) # TRUE, so cropping is fine 

# Grid sample to thin occurrences (reduce bias)
set.seed(345)
efbMichP <- dismo::gridSample(efbMichP0.2 %>% 
                                   dplyr::select(x, y),
                                 r = raster::raster(stackSag), n = 1) %>% # n = 75
  mutate(occurrence = as.integer(1),
         id = row_number())

# Get background points
# Define background function (as in Step 4)
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

# Run it
efbMichPB <- drawBackground(efbMichP, extentOb = michEqA_rast)

# See that no two points are in the same cell (n is always 1):
terra::cells(rast(michEqA_rast), vect(efbMichPB, geom = c("x", "y"))) %>%
  as.data.frame() %>% count(cell) %>% summarize(max(n))

ggplot() +
  geom_sf(data = michEqA) +
  geom_point(data = efbMichPB, aes(x = x, y = y))

# Define function to extract explanatory variables (as in Step 4)
landcoverClasses <- FedData::pal_nlcd() %>% dplyr::select(ID, Class)

extractEx <- function(i, stack) {
  bind_cols(i, terra::extract(stack, i[, 1:2])) %>%
    dplyr::select(-ID) %>%
    dplyr::rename(ID = layer) %>% # this is land cover
    left_join(landcoverClasses, by = "ID") %>%
    rename(land = Class) %>%
    dplyr::select(-ID)
}

# Run it
michPB0 <- extractEx(efbMichPB, stackSag)

# Define function to scale variables (as in Step 4)
scaleVars <- function(i) {
  
  land <- i %>% dplyr::select(c(id, land))
  
  i2 <- i %>% dplyr::select(-land)
  
  r <- ncol(i2)
  
  i2 %>%
    dplyr::mutate(across(5:all_of(r),
                         ~ c(scale(.)))) %>% #c() makes output a data frame, not a series of matrices
    dplyr::left_join(land, by = "id")
}

# Run it
michPB <- scaleVars(michPB0)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Models

# Define modeling functions exactly as in Step 5
# Function to build presence-absence (or presence-background) models
buildPAmodels <- function(i, gam = TRUE, brtRate = 0.01, bagFraction = 0.75) {
  
  landcoverClasses <- FedData::pal_nlcd() %>% dplyr::select(ID, Class)
  #recoding land cover to numeric for MaxEnt and Maxnet
  
  ### Data prep 
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    dplyr::mutate(land = as.factor(land)) %>%
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
    dplyr::mutate(land = as.factor(land)) %>%
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
  
  
  # MaxEnt requires that factors are stored as numbers (but they remain factors)
  # ANN and XGBoost require no factors, so we can use one-hot-encoding or "dummy variables"
  # (which other models are able to do on their own)
  
  lapply(datList, FUN = function(dat) {
    if("land" %in% colnames(dat)) {
      maxPred <- dat %>% dplyr::select(-occurrence) %>%
        dplyr::left_join(landcoverClasses, by = c("land" = "Class")) %>%
        dplyr::mutate(land = as.factor(ID)) %>% dplyr::select(-ID)
      
      annDat <- model.matrix(~ ., dat) %>% as.data.frame() %>% dplyr::select(-1)
      colnames(annDat) <- str_replace_all(colnames(annDat), " |\\)|\\)|,|\\/", "")
      
      xgDat <- annDat %>% as.matrix() 
      
      gamFormExclude <- c(which(colnames(dat) == "occurrence"), which(colnames(dat) == "land"))
      gamForm <- as.formula(paste0("occurrence ~", paste(
        paste0("s(", colnames(dat)[-gamFormExclude], ", k = 8)"), collapse = " + "), " + land")) } 
    
    else {
      maxPred <- dat %>% dplyr::select(-occurrence)
      
      annDat <- dat %>% as.data.frame()
      
      xgDat <- annDat %>% as.matrix() 
      
      gamForm <- as.formula(paste0("occurrence ~", paste(
        paste0("s(", colnames(dat)[-1], ", k = 8)"), collapse = " + "))) }
    
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
  
  landcoverClasses <- FedData::pal_nlcd() %>% dplyr::select(ID, Class)
  
  # Datasets with NA rows removed
  set.seed(789)
  rowrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>%
    dplyr::mutate(land = as.factor(land)) %>%
    drop_na() %>%
    dplyr::select(-c(x, y, id))
  
  rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  if("land" %in% colnames(rowrm)) {
    rowrm_test1 <- rowrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group) %>%
      dplyr::filter(land %in% unique(rowrm_train1$land))
    rowrm_test2 <- rowrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group) %>%
      dplyr::filter(land %in% unique(rowrm_train2$land))
    rowrm_test3 <- rowrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group) %>%
      dplyr::filter(land %in% unique(rowrm_train3$land))
    rowrm_test4 <- rowrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group) %>%
      dplyr::filter(land %in% unique(rowrm_train4$land))}
  else {
    rowrm_test1 <- rowrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group)
    rowrm_test2 <- rowrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group)
    rowrm_test3 <- rowrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group)
    rowrm_test4 <- rowrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group)
  }
  
  # Datasets with NA columns removed
  set.seed(987)
  colrm <- i %>%
    sample_frac(1) %>%
    dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
    dplyr::mutate(land = as.factor(land)) %>%
    dplyr::select(where(~ !any(is.na(.)))) %>% 
    dplyr::select(-c(x, y, id))
  
  colrm_train1 <- colrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
  colrm_train2 <- colrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
  colrm_train3 <- colrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
  colrm_train4 <- colrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)
  
  if("land" %in% colnames(colrm)) {
    colrm_test1 <- colrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group) %>%
      dplyr::filter(land %in% unique(colrm_train1$land))
    colrm_test2 <- colrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group) %>%
      dplyr::filter(land %in% unique(colrm_train2$land))
    colrm_test3 <- colrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group) %>%
      dplyr::filter(land %in% unique(colrm_train3$land))
    colrm_test4 <- colrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group) %>%
      dplyr::filter(land %in% unique(colrm_train4$land))}
  else {
    colrm_test1 <- colrm %>% dplyr::filter(group == 4) %>% dplyr::select(-group)
    colrm_test2 <- colrm %>% dplyr::filter(group == 3) %>% dplyr::select(-group)
    colrm_test3 <- colrm %>% dplyr::filter(group == 2) %>% dplyr::select(-group)
    colrm_test4 <- colrm %>% dplyr::filter(group == 1) %>% dplyr::select(-group)
  }
  
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
           
           if("land" %in% names(testDat)) {
             
             maxTestDat <- testDat %>%
               dplyr::left_join(landcoverClasses, by = c("land" = "Class")) %>%
               dplyr::mutate(land = as.factor(ID)) %>% dplyr::select(-ID)
             
             annTestDat <- model.matrix(~., testDat) %>% as.data.frame %>% dplyr::select(-1)
             
             colnames(annTestDat) <- str_replace_all(colnames(annTestDat), " |\\)|\\)|,|\\/", "") 
             
             xgTestDat <- annTestDat %>% as.matrix() }
           
           else {
             maxTestDat <- testDat
             
             annTestDat <- testDat
             
             xgTestDat <- testDat %>% as.matrix()}
           
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
    dplyr::select(-c(x, y, occurrence, id, land))
  
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
    dplyr::select(-c(x, y, occurrence, id, land))
  
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
    dplyr::select(-c(x, y, id, land))
  
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
    dplyr::select(-c(x, y, id, land))
  
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
                                    str_detect(evalName, "sag") ~ "sag"),
                  dataSource = case_when(str_detect(evalName, "Spec") ~ "specimen",
                                         str_detect(evalName, "Obs") ~ "observation",
                                         TRUE ~ "all"),
                  dataType = case_when(str_detect(evalName, "PB") ~ "presence-background",
                                       str_detect(evalName, "P") &
                                         !str_detect(evalName, "PB|PA") ~ "presence-only",
                                       str_detect(evalName, "PA") ~ "presence-absence",
                                       str_detect(evalName, "Ab") ~ "abundance"),
                  modelClass = case_when(str_detect(modelType, "BIOCLIM|DOMAIN") ~ "Envelope",
                                         str_detect(modelType, "GLM|GAM|MARS|MaxEnt|MaxNet") ~ "Regression",
                                         str_detect(modelType, "BRT|RF|Cforest|ANN|XGBoost") ~ "Machine learning"))
}

# Build presence-background models
michPBmodels0 <- buildPAmodels(i = michPB)
# the last BRT model failed but we still have the 3/4 needed

michPBmodels <- unlist(michPBmodels0, recursive = FALSE) %>%
  purrr::discard(is.null)

michPBmodels[["colrm_train2.ann"]] <- NULL
michPBmodels[["colrm_train4.ann"]] <- NULL

# Evaluate
michPBeval <- evaluatePAmodels(i = michPB, modList = michPBmodels, threshold = 0.5)
michPBmetrics <- convertMetrics(eval = michPBeval, modList = michPBmodels)

# Build presence-only models
michPmodels0 <- buildPmodels(michPB %>% dplyr::filter(occurrence == 1))
michPmodels <- unlist(michPmodels0, recursive = FALSE)

# Evaluate
michPeval <- evaluatePmodels(michPB %>% dplyr::filter(occurrence == 1), michPmodels, 0.5)
michPmetrics <- convertMetrics(michPeval, michPmodels)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Metrics
sagPBmetrics <- read.delim("output/sagPBmetrics.txt", header = TRUE, sep = "\t",
                             fileEncoding = "UTF-8", row.names = NULL)
sagPmetrics <- read.delim("output/sagPmetrics.txt", header = TRUE, sep = "\t",
                            fileEncoding = "UTF-8", row.names = NULL)
metricsEqA <- bind_rows(michPBmetrics, michPmetrics, sagPBmetrics, sagPmetrics) %>%
  rename(evalThreshold = threshold)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Analysis (as in Step 6)

# Only one findCorrelation threshold so we don't need to worry about that
# But we do want to select the better performing method of removing missing values (removalType)

betterMAEgroup <- metricsEqA %>%
  group_by(modelClass, modelType, scale, dataType, dataSource, removalType) %>%
  mutate(count = n()) %>%
  #exclude any where one of the groups has fewer than 3
  dplyr::filter(count >= 3) %>%
  dplyr::select(-count) %>%
  mutate(removalType = case_when(removalType == "row remove" ~ "rowrm",
                                 removalType == "col remove" ~ "colrm")) %>%
  summarize(mean = mean(mae)) %>%
  pivot_wider(names_from = removalType, values_from = mean) %>%
  mutate(betterMAE = case_when(rowrm < colrm | (is.na(colrm) & !is.na(rowrm)) ~ "row remove",
                               rowrm > colrm | (is.na(rowrm) & !is.na(colrm)) ~ "col remove")) %>%
  dplyr::select(-c(rowrm, colrm))

metricsEqA2 <- metricsEqA %>%
  right_join(betterMAEgroup, by = c("modelClass", "modelType", "scale", "dataType", "dataSource",
                                    "removalType" = "betterMAE"))

rm(betterMAEgroup)

# We're really just asking Question 3 again (effect of scale)
question3 <- bind_rows(
  metricsEqA2 %>%
    dplyr::filter(dataType != "presence-absence" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(mae ~ scale, paired = FALSE, p.adjust.method = "bonferroni") %>%
    left_join(metricsEqA2 %>%
                group_by(modelType, scale) %>%
                summarize(meanMAE = mean(mae)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanMAE"),
              by = "modelType"),
  metricsEqA2 %>%
    dplyr::filter(dataType != "presence-absence" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(sensitivity ~ scale, paired = FALSE, p.adjust.method = "bonferroni") %>%
    left_join(metricsEqA2 %>%
                group_by(modelType, scale) %>%
                summarize(meanSens = mean(sensitivity)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanSens"),
              by = "modelType"),
  metricsEqA2 %>%
    dplyr::filter(dataType == "presence-background" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(specificity ~ scale, paired = FALSE) %>%
    left_join(metricsEqA2 %>%
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

# Supplementary table will not be read out, but can be created below
suppTable3 <- question3 %>%
  mutate(across(.cols = c(mich, sag),
                ~ round(., 4))) %>%
  mutate(sig = str_to_sentence(sig)) %>%
  dplyr::select("Model type" = modelType, "Metric" = metric, 
                "Group 1 n" = n1, "Group 2 n" = n2,
                "Wilcox statistic" = statistic, "P" = p, "Significance" = sig,
                "Mean value for Michigan" = mich, 
                "Mean value for Saginaw Bay" = sag)

# Visualize differences in each metric between scales
question3Long <- question3 %>%
  pivot_longer(cols = c(mich, sag),
               names_to = "scale",
               values_to = "meanValue") %>%
  mutate(meanValue = case_when(metric == "Mean absolute error (MAE)" ~ 1 - meanValue,
                               TRUE ~ meanValue),
         metric = case_when(metric == "Mean absolute error (MAE)" ~ "1 - mean absolute error (MAE)",
                            TRUE ~ metric))

# Plot (will not be written out)
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
                     expand = c(0.025, 0.025)) +
  labs(x = "Mean value", y = "Model type",
       color = "Scale") +
  scale_color_manual(values = c("#37406B", "#96BBDA"),
                     labels = c("Michigan", "Saginaw Bay")) +
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
        legend.key.width = unit(0.25, "in"))

plot3

# Takeaway when Michigan area = Saginaw Bay area:
# They perform very similarly all around
# GLM MAE and MaxNet specificity are better for Michigan
# DOMAIN is better for Saginaw Bay

################################################################################

rm(list = ls())

######### 2) Michigan scale = same sample size as Saginaw Bay ("EqN") ##########

stackMich <- terra::rast("output/stackMich.tif") #findCorrelation threshold = 0.7
stackMichB <- terra::rast("output/stackMichB.tif") #findCorrelation threshold = 0.5

mich <- terra::subset(stackMich, 1) %>% raster::raster() # can be any of the layers

# Sample from Michigan data to get same number as Saginaw Bay presences
# We can draw from the existing data "efbMichP" to simplify things

efbMichP_2.0 <- read.delim("output/efbMichP.txt", header = TRUE,
                       sep = "\t", quote = "", fileEncoding = "UTF-8")

sagPB <- read.delim("output/sagPB.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")
num <- sagPB %>% dplyr::filter(occurrence == 1) %>% count() %>% pull()

set.seed(345)
efbMichP_2 <- efbMichP_2.0 %>%
  sample_frac(1) %>%
  head(n = num) # n = 100

# Draw background points
efbMichPB_2 <- drawBackground(efbMichP_2, extentOb = mich)

# Extract variables at both findCorrelationThresholds
michPB_2_0 <- extractEx(efbMichPB_2, stackMich)
michPB_2_B_0 <- extractEx(efbMichPB_2, stackMichB)

# Scale variables
michPB_2 <- scaleVars(michPB_2_0)
michPB_2_B <- scaleVars(michPB_2_B_0)

# And get the presence-only datsets
michP_2 <- michPB_2 %>% dplyr::filter(occurrence == 1)
michP_2_B <- michPB_2_B %>% dplyr::filter(occurrence == 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# findCorrelation threshold = 0.7
# Build presence-background models
michPB_2models0 <- buildPAmodels(michPB_2)

michPB_2models <- unlist(michPB_2models0, recursive = FALSE)

# Evaluate
michPB_2eval <- evaluatePAmodels(i = michPB_2, modList = michPB_2models, threshold = 0.5)
michPB_2metrics <- convertMetrics(eval = michPB_2eval, modList = michPB_2models)

# Build presence-only models
michP_2models0 <- buildPmodels(michP_2)
michP_2models <- unlist(michP_2models0, recursive = FALSE)

# Evaluate
michP_2eval <- evaluatePmodels(michP_2, michP_2models, 0.5)
michP_2metrics <- convertMetrics(michP_2eval, michP_2models)

# findCorrelation threshold = 0.5
# Build presence-background models
michPB_2_Bmodels0 <- buildPAmodels(michPB_2_B)

michPB_2_Bmodels <- unlist(michPB_2_Bmodels0, recursive = FALSE)

# Evaluate
michPB_2_Beval <- evaluatePAmodels(i = michPB_2_B, modList = michPB_2_Bmodels, threshold = 0.5)
michPB_2_Bmetrics <- convertMetrics(eval = michPB_2_Beval, modList = michPB_2_Bmodels)

# Build presence-only models
michP_2_Bmodels0 <- buildPmodels(michP_2_B)
michP_2_Bmodels <- unlist(michP_2_Bmodels0, recursive = FALSE)

# Evaluate
michP_2_Beval <- evaluatePmodels(michP_2_B, michP_2_Bmodels, 0.5)
michP_2_Bmetrics <- convertMetrics(michP_2_Beval, michP_2_Bmodels)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

### Metrics
metricsEqN <- bind_rows(michPB_2metrics, michP_2metrics, sagPBmetrics, sagPmetrics) %>%
  mutate(varRemovalThreshold = 0.7) %>%
  bind_rows(michPB_2_Bmetrics %>% mutate(varRemovalThreshold = 0.5)) %>%
  bind_rows(michP_2_Bmetrics %>% mutate(varRemovalThreshold = 0.5)) %>%
  rename(evalThreshold = threshold)

# Choose the better method of removing missing values in each case
betterMAEgroup <- metricsEqN %>%
  group_by(modelClass, modelType, scale, dataType, dataSource, removalType, varRemovalThreshold) %>%
  mutate(count = n()) %>%
  #exclude any where one of the groups has fewer than 3
  dplyr::filter(count >= 3) %>%
  dplyr::select(-count) %>%
  mutate(removalType = case_when(removalType == "row remove" ~ "rowrm",
                                 removalType == "col remove" ~ "colrm")) %>%
  summarize(mean = mean(mae)) %>%
  pivot_wider(names_from = removalType, values_from = mean) %>%
  mutate(betterMAE = case_when(rowrm < colrm | (is.na(colrm) & !is.na(rowrm)) ~ "row remove",
                               rowrm > colrm | (is.na(rowrm) & !is.na(colrm)) ~ "col remove")) %>%
  dplyr::select(-c(rowrm, colrm))

metricsEqN2 <- metricsEqN %>%
  right_join(betterMAEgroup, by = c("modelClass", "modelType", "scale", "dataType", "dataSource",
                                    "varRemovalThreshold", "removalType" = "betterMAE"))

rm(betterMAEgroup)

# And the better findCorrelation threshold in each case
betterMAEgroup2 <-  metricsEqN2 %>%
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

metricsEqN3 <- metricsEqN2 %>%
  right_join(betterMAEgroup2, by = c("modelClass", "modelType", "scale", "dataType", "dataSource",
                                     "varRemovalThreshold" = "betterMAE"))

rm(betterMAEgroup2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Analysis (again, it's just Question 3)
question3_2 <- bind_rows(
  metricsEqN3 %>%
    dplyr::filter(dataType != "presence-absence" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(mae ~ scale, paired = FALSE, p.adjust.method = "bonferroni") %>%
    left_join(metricsEqN3 %>%
                group_by(modelType, scale) %>%
                summarize(meanMAE = mean(mae)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanMAE"),
              by = "modelType"),
  metricsEqN3 %>%
    dplyr::filter(dataType != "presence-absence" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(sensitivity ~ scale, paired = FALSE, p.adjust.method = "bonferroni") %>%
    left_join(metricsEqN3 %>%
                group_by(modelType, scale) %>%
                summarize(meanSens = mean(sensitivity)) %>%
                pivot_wider(names_from = "scale",
                            values_from = "meanSens"),
              by = "modelType"),
  metricsEqN3 %>%
    dplyr::filter(dataType == "presence-background" &
                    dataSource == "all") %>%
    group_by(modelType) %>%
    rstatix::wilcox_test(specificity ~ scale, paired = FALSE) %>%
    left_join(metricsEqN3 %>%
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

# Supplementary table will not be read out, but can be created below
suppTable3_2 <- question3_2 %>%
  mutate(across(.cols = c(mich, sag),
                ~ round(., 4))) %>%
  mutate(sig = str_to_sentence(sig)) %>%
  dplyr::select("Model type" = modelType, "Metric" = metric, 
                "Group 1 n" = n1, "Group 2 n" = n2,
                "Wilcox statistic" = statistic, "P" = p, "Significance" = sig,
                "Mean value for Michigan" = mich, 
                "Mean value for Saginaw Bay" = sag)

# Visualize differences in each metric between scales
question3_2Long <- question3_2 %>%
  pivot_longer(cols = c(mich, sag),
               names_to = "scale",
               values_to = "meanValue") %>%
  mutate(meanValue = case_when(metric == "Mean absolute error (MAE)" ~ 1 - meanValue,
                               TRUE ~ meanValue),
         metric = case_when(metric == "Mean absolute error (MAE)" ~ "1 - mean absolute error (MAE)",
                            TRUE ~ metric))

# Plot (will not be written out)
plot3_2 <- ggplot(question3_2Long, aes(x = meanValue, 
                                   y = reorder(modelType, desc(modelType)))) +
  geom_line() +
  geom_point(aes(color = scale, group = metric),
             size = 4) +
  geom_text(data = question3_2Long %>%
              mutate(sig = case_when(sig == "s" ~ "*",
                                     sig == "ns" ~ "")) %>%
              arrange(metric, modelType, meanValue) %>%
              distinct(modelType, metric, sig, .keep_all = TRUE),
            aes(label = sig), position = position_nudge(x = -0.085, y = -0.1),
            size = 6) +
  facet_wrap(vars(metric)) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
                     expand = c(0.025, 0.025)) +
  labs(x = "Mean value", y = "Model type",
       color = "Scale") +
  scale_color_manual(values = c("#37406B", "#96BBDA"),
                     labels = c("Michigan", "Saginaw Bay")) +
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
        legend.key.width = unit(0.25, "in"))

plot3_2

# Takeaway when Michigan sample size = Saginaw Bay sample size:
# Opportunistic data over a large area (Michigan)
# often outperform systematic data over a smaller area (Saginaw Bay)

# Comparing the equal area and equal sample size methods,
# we can infer that the larger geographic area of Michigan is probably
# driving its better performance compared with Saginaw Bay

################################################################################

rm(list = ls())

