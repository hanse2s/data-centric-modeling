# Data-centric SDM Step 5) Fit and evaluate models
## Sara Hansen
## Modified Juy 25, 2024

library(conflicted)
library(tidyverse)
library(tools)
library(terra)
library(mgcv)
library(earth)
library(dismo)
# To use the dismo::maxent function, you need to have the MaxEnt software installed
# https://biodiversityinformatics.amnh.org/open_source/maxent/
# You should also make sure Java is installed and compatible with your version of R
# https://www.java.com/en/download/
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

############################### CREATE FUNCTIONS ###############################

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

################################################################################

#Build first set of models (correlation threshold 0.7)

############## BUILD PRESENCE-BACKGROUND / PRESENCE-ABSENCE MODELS #############

# Scale = Michigan, source = specimens + observations, type = presence + background

# Read in data
michPB <- read.delim("output/michPB.txt", header = TRUE,
                     sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michPBmodels0 <- buildPAmodels(i = michPB)

# Unnest the list for the evaluate functions
michPBmodels <- unlist(michPBmodels0, recursive = FALSE)

# Remove any models that failed to build fully
michPBmodels[["colrm_train1.ann"]] <- NULL
michPBmodels[["colrm_train2.ann"]] <- NULL
michPBmodels[["colrm_train3.ann"]] <- NULL
michPBmodels[["colrm_train4.ann"]] <- NULL

# Evaluate at threshold of 0.5
michPBeval <- evaluatePAmodels(i = michPB, modList = michPBmodels, threshold = 0.5)

# Convert metrics to a data frame for saving and analyzing 
michPBmetrics <- convertMetrics(eval = michPBeval, modList = michPBmodels)

# Write out metrics
write_delim(michPBmetrics, "output/michPBmetrics.txt", delim = "\t", quote = "none")

# You can write out the models too
#saveRDS(michPBmodels, "output/michPBmodels.rds")
#saveRDS(michPBeval, "output/michPBeval.rds")

rm(list = ls(pattern = "michPB"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Michigan, source = specimens, type = presence + background

# Read in data
michSpecPB <- read.delim("output/michSpecPB.txt", header = TRUE,
                         sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
# "Error: Model has more coefficients than data" means we have to build GAM separately
# Because at least one model is failing to build
michSpecPBmodels0 <- buildPAmodels(michSpecPB, gam = FALSE)

# GAM models cannot be built separately due to smoothing increasing
# the number of coefficients too much for this dataset size

michSpecPBmodels <- unlist(michSpecPBmodels0, recursive = FALSE) %>%
  purrr::discard(is.null) # note this may remove models other than GAM, if they failed

# Evaluate
michSpecPBeval <- evaluatePAmodels(michSpecPB, michSpecPBmodels, 0.5)
michSpecPBmetrics <- convertMetrics(michSpecPBeval, michSpecPBmodels)

# Write out metrics
write_delim(michSpecPBmetrics, "output/michSpecPBmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michSpecPBmodels, "output/michSpecPBmodels.rds")
#saveRDS(michSpecPBeval, "output/michSpecPBeval.rds")

rm(list = ls(pattern = "michSpecPB"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Michigan, source = observations, type = presence + background

# Read in data
michObsPB <- read.delim("output/michObsPB.txt", header = TRUE,
                        sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michObsPBmodels0 <- buildPAmodels(michObsPB)
michObsPBmodels <- unlist(michObsPBmodels0, recursive = FALSE)

michObsPBmodels[["colrm_train2.ann"]] <- NULL

# Evaluate
michObsPBeval <- evaluatePAmodels(michObsPB, michObsPBmodels, 0.5)
michObsPBmetrics <- convertMetrics(michObsPBeval, michObsPBmodels)

# Write out metrics
write_delim(michObsPBmetrics, "output/michObsPBmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michObsPBmodels, "output/michObsPBmodels.rds")
#saveRDS(michObsPBeval, "output/michObsPBeval.rds")

rm(list = ls(pattern = "michObsPB"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Saginaw Bay, type = presence + background

# Read in data
sagPB <- read.delim("output/sagPB.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
sagPBmodels0 <- buildPAmodels(sagPB, gam = FALSE)

# Build GAM models separately
set.seed(789)
rowrm <- sagPB %>%
  sample_frac(1) %>%
  dplyr::mutate(group = rep(c(1:4), length = n())) %>%
  dplyr::mutate(land = as.factor(land)) %>%
  drop_na() %>%
  dplyr::select(-c(x, y, id))

rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

set.seed(789)
colrm <- sagPB %>%
  sample_frac(1) %>%
  dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
  dplyr::mutate(land = as.factor(land)) %>%
  dplyr::select(where(~ !any(is.na(.)))) %>% 
  dplyr::select(-c(x, y, id))

colrm_train1 <- colrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
colrm_train2 <- colrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
colrm_train3 <- colrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
colrm_train4 <- colrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

lapply(list(rowrm_train1 = rowrm_train1,
            rowrm_train2 = rowrm_train2,
            rowrm_train3 = rowrm_train3,
            rowrm_train4 = rowrm_train4,
            colrm_train1 = colrm_train1,
            colrm_train2 = colrm_train2,
            colrm_train3 = colrm_train3,
            colrm_train4 = colrm_train4),
       FUN = function(dat) {dat %>% summarise_all(n_distinct)} )

# No changes to parameters, these just didn't want to run in the loop
set.seed(865)
rowrm_train1.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                               s(nit, k = 8) + s(pop, k = 8) + s(wc10, k = 8) + land,
                             family = binomial(link = "logit"),
                             data = rowrm_train1, method = "REML")

set.seed(865)
rowrm_train2.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                               s(nit, k = 8) + s(pop, k = 8) + s(wc10, k = 8) + land,
                             family = binomial(link = "logit"),
                             data = rowrm_train2, method = "REML")

set.seed(865)
rowrm_train3.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(nit, k = 8) + s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train3, method = "REML")

set.seed(865) 
rowrm_train4.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(nit, k = 8) + s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train4, method = "REML")

set.seed(865)
colrm_train1.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train1, method = "REML")

set.seed(865)
colrm_train2.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train2, method = "REML")

set.seed(865)
colrm_train3.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train3, method = "REML")

set.seed(865)
colrm_train4.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train4, method = "REML")
  
sagPBGAMmodels <- list(rowrm_train1.gam = rowrm_train1.gam, rowrm_train2.gam = rowrm_train2.gam,
                       rowrm_train3.gam = rowrm_train3.gam, rowrm_train4.gam = rowrm_train4.gam,
                       colrm_train1.gam = colrm_train1.gam, colrm_train2.gam = colrm_train2.gam,
                       colrm_train3.gam = colrm_train3.gam, colrm_train4.gam = colrm_train4.gam)

rm(list = ls(pattern = "rowrm|colrm"))

sagPBmodels <- unlist(sagPBmodels0, recursive = FALSE) %>% 
  purrr::discard(is.null) %>%
  append(sagPBGAMmodels)

sagPBmodels[["colrm_train1.ann"]] <- NULL
sagPBmodels[["colrm_train2.ann"]] <- NULL

# Evaluate
sagPBeval <- evaluatePAmodels(sagPB, sagPBmodels, 0.5)
sagPBmetrics <- convertMetrics(sagPBeval, sagPBmodels)

# Write out metrics
write_delim(sagPBmetrics, "output/sagPBmetrics.txt", delim = "\t", quote = "none")

#saveRDS(sagPBmodels, "output/sagPBmodels.rds")
#saveRDS(sagPBeval, "output/sagPBeval.rds")

rm(list = ls(pattern = "sagPB"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Saginaw Bay, type = presence + absence

# Read in data
sagPA <- read.delim("output/sagPA.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
# If BRT models do not build at all, set a learning rate smaller than 0.01
# The error will be "restart model with a smaller learning rate or smaller step size..."
# and/or try increasing bagFraction a little
sagPAmodels0 <- buildPAmodels(sagPA, gam = FALSE, brtRate = 0.001, bagFraction = 0.8)

# Build GAM models separately
set.seed(789)
rowrm <- sagPA %>%
  sample_frac(1) %>%
  dplyr::mutate(group = rep(c(1:4), length = n())) %>%
  dplyr::mutate(land = as.factor(land)) %>%
  drop_na() %>%
  dplyr::select(-c(x, y, id))

rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

set.seed(789)
colrm <- sagPA %>%
  sample_frac(1) %>%
  dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
  dplyr::mutate(land = as.factor(land)) %>%
  dplyr::select(where(~ !any(is.na(.)))) %>% 
  dplyr::select(-c(x, y, id))

colrm_train1 <- colrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
colrm_train2 <- colrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
colrm_train3 <- colrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
colrm_train4 <- colrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

lapply(list(rowrm_train1 = rowrm_train1,
            rowrm_train2 = rowrm_train2,
            rowrm_train3 = rowrm_train3,
            rowrm_train4 = rowrm_train4,
            colrm_train1 = colrm_train1,
            colrm_train2 = colrm_train2,
            colrm_train3 = colrm_train3,
            colrm_train4 = colrm_train4),
       FUN = function(dat) {dat %>% summarise_all(n_distinct)} )

# Only changes to parameters are elevation k value
set.seed(865)
rowrm_train3.gam = mgcv::gam(occurrence ~  s(boat, k = 8) +  s(elev, k = 7) +
                                 s(nit, k = 8) + s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train3, method = "REML")

set.seed(865)
rowrm_train4.gam = mgcv::gam(occurrence ~  s(boat, k = 8) +  s(elev, k = 7) +
                                 s(nit, k = 8) + s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train4, method = "REML")
set.seed(865)
colrm_train1.gam = mgcv::gam(occurrence ~  s(boat, k = 8) +  s(elev, k = 4) + 
                                 s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train1, method = "REML")

set.seed(865)
colrm_train2.gam = mgcv::gam(occurrence ~  s(boat, k = 8) +  s(elev, k = 6) + 
                                 s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train2, method = "REML")

set.seed(865)
colrm_train3.gam = mgcv::gam(occurrence ~  s(boat, k = 8) +  s(elev, k = 7) + 
                                 s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train3, method = "REML")

set.seed(865)
colrm_train4.gam = mgcv::gam(occurrence ~  s(boat, k = 8) +  s(elev, k = 7) + 
                                 s(pop, k = 8) + s(wc10, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train4, method = "REML")

sagPAGAMmodels <- list(rowrm_train3.gam = rowrm_train3.gam, rowrm_train4.gam = rowrm_train4.gam,
                       colrm_train1.gam = colrm_train1.gam, colrm_train2.gam = colrm_train2.gam,
                       colrm_train3.gam = colrm_train3.gam, colrm_train4.gam = colrm_train4.gam)

rm(list = ls(pattern = "rowrm|colrm"))

sagPAmodels <- unlist(sagPAmodels0, recursive = FALSE) %>%
  purrr::discard(is.null) %>%
  append(sagPAGAMmodels)

# Evaluate
sagPAeval <- evaluatePAmodels(sagPA, sagPAmodels, 0.5)
sagPAmetrics <- convertMetrics(sagPAeval, sagPAmodels)

# Write out metrics
write_delim(sagPAmetrics, "output/sagPAmetrics.txt", delim = "\t", quote = "none")

#saveRDS(sagPAmodels, "output/sagPAmodels.rds")
#saveRDS(sagPAeval, "output/sagPAeval.rds")

rm(list = ls(pattern = "sagPA"))

################################################################################

######################### BUILD PRESENCE-ONLY MODELS ############################

# Note that neither model will accept categorical variables,
# so land cover cannot be included as a predictor 

# Scale = Michigan, source = specimens + observations, type = presence
# Read in data
michP <- read.delim("output/michP.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")
 
# Build models
michPmodels0 <- buildPmodels(michP)
michPmodels <- unlist(michPmodels0, recursive = FALSE)

# Evaluate
michPeval <- evaluatePmodels(michP, michPmodels, 0.5)
michPmetrics <- convertMetrics(michPeval, michPmodels)

# Write out metrics
write_delim(michPmetrics, "output/michPmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michPmodels, "output/michPmodels.rds")
#saveRDS(michPeval, "output/michPeval.rds")

rm(list = ls(pattern = "michP"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Michigan, source = specimens, type = presence
# Read in data
michSpecP <- read.delim("output/michSpecP.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michSpecPmodels0 <- buildPmodels(michSpecP)
michSpecPmodels <- unlist(michSpecPmodels0, recursive = FALSE)

# Evaluate
michSpecPeval <- evaluatePmodels(michSpecP, michSpecPmodels, 0.5)
michSpecPmetrics <- convertMetrics(michSpecPeval, michSpecPmodels)

# Write out metrics
write_delim(michSpecPmetrics, "output/michSpecPmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michSpecPmodels, "output/michSpecPmodels.rds")
#saveRDS(michSpecPeval, "output/michSpecPeval.rds")

rm(list = ls(pattern = "michSpecP"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Michigan, source = observations, type = presence
# Read in data
michObsP <- read.delim("output/michObsP.txt", header = TRUE,
                        sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michObsPmodels0 <- buildPmodels(michObsP)
michObsPmodels <- unlist(michObsPmodels0, recursive = FALSE)

# Evaluate
michObsPeval <- evaluatePmodels(michObsP, michObsPmodels, 0.5)
michObsPmetrics <- convertMetrics(michObsPeval, michObsPmodels)

# Write out metrics
write_delim(michObsPmetrics, "output/michObsPmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michObsPmodels, "output/michObsPmodels.rds")
#saveRDS(michObsPeval, "output/michObsPeval.rds")

rm(list = ls(pattern = "michObsP"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Saginaw Bay, type = presence
# Read in data
sagP <- read.delim("output/sagP.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
sagPmodels0 <- buildPmodels(sagP)
sagPmodels <- unlist(sagPmodels0, recursive = FALSE)

# Evaluate
sagPeval <- evaluatePmodels(sagP, sagPmodels, 0.5)
sagPmetrics <- convertMetrics(sagPeval, sagPmodels)

# Write out metrics
write_delim(sagPmetrics, "output/sagPmetrics.txt", delim = "\t", quote = "none")

#saveRDS(sagPmodels, "output/sagPmodels.rds")
#saveRDS(sagPeval, "output/sagPeval.rds")

rm(list = ls(pattern = "sagP"))

################################################################################

#Build second set of models (correlation threshold 0.5)

############## BUILD PRESENCE-BACKGROUND / PRESENCE-ABSENCE MODELS #############

# Scale = Michigan, source = specimens + observations, type = presence-background
# Read in data
michPB_B <- read.delim("output/michPB_B.txt", header = TRUE,
                       sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michPB_Bmodels0 <- buildPAmodels(michPB_B)
michPB_Bmodels <- unlist(michPB_Bmodels0, recursive = FALSE)

michPB_Bmodels[["colrm_train1.ann"]] <- NULL

# Evaluate models
michPB_Beval <- evaluatePAmodels(michPB_B, michPB_Bmodels, threshold = 0.5)
michPB_Bmetrics <- convertMetrics(michPB_Beval, michPB_Bmodels)

# Write out metrics
write_delim(michPB_Bmetrics, "output/michPB_Bmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michPB_Bmodels, "output/michPB_Bmodels.rds")
#saveRDS(michPB_Beval, "output/michPB_Beval.rds")

rm(list = ls(pattern = "michPB_B"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Michigan, source = specimens, type = presence-background
# Read in data
michSpecPB_B <- read.delim("output/michSpecPB_B.txt", header = TRUE,
                       sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michSpecPB_Bmodels0 <- buildPAmodels(michSpecPB_B, gam = FALSE)

# Build GAM models separately
set.seed(789)
rowrm <- michSpecPB_B %>%
  sample_frac(1) %>%
  dplyr::mutate(group = rep(c(1:4), length = n())) %>%
  dplyr::mutate(land = as.factor(land)) %>%
  drop_na() %>%
  dplyr::select(-c(x, y, id))

rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

set.seed(789)
colrm <- michSpecPB_B %>%
  sample_frac(1) %>%
  dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
  dplyr::mutate(land = as.factor(land)) %>%
  dplyr::select(where(~ !any(is.na(.)))) %>% 
  dplyr::select(-c(x, y, id))

colrm_train1 <- colrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
colrm_train2 <- colrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
colrm_train3 <- colrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
colrm_train4 <- colrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

lapply(list(rowrm_train1 = rowrm_train1,
            rowrm_train2 = rowrm_train2,
            rowrm_train3 = rowrm_train3,
            rowrm_train4 = rowrm_train4,
            colrm_train1 = colrm_train1,
            colrm_train2 = colrm_train2,
            colrm_train3 = colrm_train3,
            colrm_train4 = colrm_train4),
       FUN = function(dat) {dat %>% summarise_all(n_distinct)} )

# Only changes is not modeling rowrm_train2 -> no changes to parameters here
set.seed(865)
rowrm_train1.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(nit, k = 8) + s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train1, method = "REML")

set.seed(865)
rowrm_train3.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(nit, k = 8) + s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train3, method = "REML")

set.seed(865)
rowrm_train4.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(nit, k = 8) + s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train4, method = "REML")

set.seed(865)
colrm_train1.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train1, method = "REML")

set.seed(865)
colrm_train2.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train2, method = "REML")

set.seed(865)
colrm_train3.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train3, method = "REML")

set.seed(865)
colrm_train4.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = colrm_train4, method = "REML")

michSpecPB_BGAMmodels <- list(rowrm_train1.gam = rowrm_train1.gam,
                              rowrm_train3.gam = rowrm_train3.gam, rowrm_train4.gam = rowrm_train4.gam,
                              colrm_train1.gam = colrm_train1.gam, colrm_train2.gam = colrm_train2.gam,
                              colrm_train3.gam = colrm_train3.gam, colrm_train4.gam = colrm_train4.gam)

rm(list = ls(pattern = "rowrm|colrm"))

michSpecPB_Bmodels <- unlist(michSpecPB_Bmodels0, recursive = FALSE) %>% 
  purrr::discard(is.null) %>%
  append(michSpecPB_BGAMmodels)

# Evaluate models
michSpecPB_Beval <- evaluatePAmodels(michSpecPB_B, michSpecPB_Bmodels, threshold = 0.5)
michSpecPB_Bmetrics <- convertMetrics(michSpecPB_Beval, michSpecPB_Bmodels)

# Write out metrics
write_delim(michSpecPB_Bmetrics, "output/michSpecPB_Bmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michSpecPB_Bmodels, "output/michSpecPB_Bmodels.rds")
#saveRDS(michSpecPB_Beval, "output/michSpecPB_Beval.rds")

rm(list = ls(pattern = "michSpecPB_B"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Michigan, source = observations, type = presence-background
# Read in data
michObsPB_B <- read.delim("output/michObsPB_B.txt", header = TRUE,
                           sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michObsPB_Bmodels0 <- buildPAmodels(michObsPB_B, gam = FALSE)

# Build GAM models separately
set.seed(789)
rowrm <- michObsPB_B %>%
  sample_frac(1) %>%
  dplyr::mutate(group = rep(c(1:4), length = n())) %>%
  dplyr::mutate(land = as.factor(land)) %>%
  drop_na() %>%
  dplyr::select(-c(x, y, id))

rowrm_train1 <- rowrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
rowrm_train2 <- rowrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
rowrm_train3 <- rowrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
rowrm_train4 <- rowrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

set.seed(789)
colrm <- michObsPB_B %>%
  sample_frac(1) %>%
  dplyr::mutate(group = rep(c(1:4), length = n())) %>% 
  dplyr::mutate(land = as.factor(land)) %>%
  dplyr::select(where(~ !any(is.na(.)))) %>% 
  dplyr::select(-c(x, y, id))

colrm_train1 <- colrm %>% dplyr::filter(group != 4) %>% dplyr::select(-group)
colrm_train2 <- colrm %>% dplyr::filter(group != 3) %>% dplyr::select(-group)
colrm_train3 <- colrm %>% dplyr::filter(group != 2) %>% dplyr::select(-group)
colrm_train4 <- colrm %>% dplyr::filter(group != 1) %>% dplyr::select(-group)

lapply(list(rowrm_train1 = rowrm_train1,
            rowrm_train2 = rowrm_train2,
            rowrm_train3 = rowrm_train3,
            rowrm_train4 = rowrm_train4,
            colrm_train1 = colrm_train1,
            colrm_train2 = colrm_train2,
            colrm_train3 = colrm_train3,
            colrm_train4 = colrm_train4),
       FUN = function(dat) {dat %>% summarise_all(n_distinct)} )

# Only adjustment here is not modeling rowrm_train2 -> parameters are the same
set.seed(865)
rowrm_train1.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(nit, k = 8) + s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train1, method = "REML")

set.seed(865)  
rowrm_train3.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(nit, k = 8) + s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train3, method = "REML")

set.seed(865)
rowrm_train4.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(nit, k = 8) + s(phos, k = 8) + s(pop, k = 8) + land,
                               family = binomial(link = "logit"),
                               data = rowrm_train4, method = "REML")

set.seed(865)
colrm_train1.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(pop, k = 8),
                               family = binomial(link = "logit"),
                               data = colrm_train1, method = "REML")

set.seed(865)  
colrm_train2.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(pop, k = 8),
                               family = binomial(link = "logit"),
                               data = colrm_train2, method = "REML")

set.seed(865)
colrm_train3.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(pop, k = 8),
                               family = binomial(link = "logit"),
                               data = colrm_train3, method = "REML")

set.seed(865) 
colrm_train4.gam = mgcv::gam(occurrence ~ s(boat, k = 8) + s(elev, k = 8) +
                                 s(pop, k = 8),
                               family = binomial(link = "logit"),
                               data = colrm_train4, method = "REML")
  
  
michObsPB_BGAMmodels <- list(rowrm_train1.gam = rowrm_train1.gam,
                             rowrm_train3.gam = rowrm_train3.gam, rowrm_train4.gam = rowrm_train4.gam,
                             colrm_train1.gam = colrm_train1.gam, colrm_train2.gam = colrm_train2.gam,
                             colrm_train3.gam = colrm_train3.gam, colrm_train4.gam = colrm_train4.gam)
    

rm(list = ls(pattern = "rowrm|colrm"))

michObsPB_Bmodels <- unlist(michObsPB_Bmodels0, recursive = FALSE) %>%
  purrr::discard(is.null) %>%
  append(michObsPB_BGAMmodels)

michObsPB_Bmodels[["colrm_train1.ann"]] <- NULL
michObsPB_Bmodels[["colrm_train2.ann"]] <- NULL
michObsPB_Bmodels[["colrm_train4.ann"]] <- NULL

# Evaluate models
michObsPB_Beval <- evaluatePAmodels(michObsPB_B, michObsPB_Bmodels, threshold = 0.5)
michObsPB_Bmetrics <- convertMetrics(michObsPB_Beval, michObsPB_Bmodels)

# Write out metrics
write_delim(michObsPB_Bmetrics, "output/michObsPB_Bmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michObsPB_Bmodels, "output/michObsPB_Bmodels.rds")
#saveRDS(michObsPB_Beval, "output/michObsPB_Beval.rds")

rm(list = ls(pattern = "michObsPB_B"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Saginaw Bay data has only one threshold, no "_B" models

################################################################################

######################### BUILD PRESENCE-ONLY MODELS ###########################

# Scale = Michigan, source = specimens + observations, type = presence
# Read in data
michP_B <- read.delim("output/michP_B.txt", header = TRUE,
                    sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michP_Bmodels0 <- buildPmodels(michP_B)
michP_Bmodels <- unlist(michP_Bmodels0, recursive = FALSE)
rm(michP_Bmodels0)

# Evaluate models
michP_Beval <- evaluatePmodels(michP_B, michP_Bmodels, 0.5)
michP_Bmetrics <- convertMetrics(michP_Beval, michP_Bmodels)

# Write out metrics
write_delim(michP_Bmetrics, "output/michP_Bmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michP_Bmodels, "output/michP_Bmodels.rds")
#saveRDS(michP_Beval, "output/michP_Beval.rds")

rm(list = ls(pattern = "michP_B"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Michigan, source = specimens, type = presence
# Read in data
michSpecP_B <- read.delim("output/michSpecP_B.txt", header = TRUE,
                      sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michSpecP_Bmodels0 <- buildPmodels(michSpecP_B)
michSpecP_Bmodels <- unlist(michSpecP_Bmodels0, recursive = FALSE)
rm(michSpecP_Bmodels0)

# Evaluate models
michSpecP_Beval <- evaluatePmodels(michSpecP_B, michSpecP_Bmodels, 0.5)
michSpecP_Bmetrics <- convertMetrics(michSpecP_Beval, michSpecP_Bmodels)

# Write out metrics
write_delim(michSpecP_Bmetrics, "output/michSpecP_Bmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michSpecP_Bmodels, "output/michSpecP_Bmodels.rds")
#saveRDS(michSpecP_Beval, "output/michSpecP_Beval.rds")

rm(list = ls(pattern = "michSpecP_B"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Scale = Michigan, source = observations, type = presence
# Read in data
michObsP_B <- read.delim("output/michObsP_B.txt", header = TRUE,
                          sep = "\t", quote = "", fileEncoding = "UTF-8-BOM")

# Build models
michObsP_Bmodels0 <- buildPmodels(michObsP_B)
michObsP_Bmodels <- unlist(michObsP_Bmodels0, recursive = FALSE)
rm(michObsP_Bmodels0)

# Evaluate models
michObsP_Beval <- evaluatePmodels(michObsP_B, michObsP_Bmodels, 0.5)
michObsP_Bmetrics <- convertMetrics(michObsP_Beval, michObsP_Bmodels)

# Write out metrics
write_delim(michObsP_Bmetrics, "output/michObsP_Bmetrics.txt", delim = "\t", quote = "none")

#saveRDS(michObsP_Bmodels, "output/michObsP_Bmodels.rds")
#saveRDS(michObsP_Beval, "output/michObsP_Beval.rds")

rm(list = ls(pattern = "michObsP_B"))

################################################################################

rm(list = ls())



