
# maybe keep ac model + rf for comparison (at least at the predictive level, e.g compare auc,thenfor extrapolation keep only rf model)
# to restart with (after 22/9)
# https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab15_SpatialRegression.html
# https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab10_PCA.html#Raster_data

# libraries ----
library(sp)
library(sdmvspecies)
library(tidyverse)
library(rgdal)
library(maptools)
library(rgeos)
library(rJava)
library(dismo)
library(sf)
library(raster)
library(virtualspecies)
library(spdplyr)
library(ENMeval)
library(lubridate)
library(amt)
library(osmextract)
library(rasterVis)
library(viridis)
library(caret)
library(CAST)
library(gridExtra)
library(knitr)
library(grid)
library(latticeExtra)
library(here)
library(nlme)
library(gstat)
library(tidymodels)
library(sfheaders)
library(sjmisc)

# data import ----

soci <- readRDS("data/derived_data/SocTbl/sampling1h/soc_indirect_wb.rds") %>% rename(indirect=layer) %>%  #rename(indirect=Indirect_interactions) #%>% dplyr::select(study,cells, withbetw)
  mutate(indirect = rescale(indirect))
socd <- readRDS("data/derived_data/SocTbl/sampling1h/soc_direct_wb.rds") %>%  rename(direct=layer) %>% dplyr::select(study,cells, direct) %>% 
  mutate(direct = rescale(direct))
#socwb <- readRDS("data/derived_data/SocTbl/sampling1h/soc_withbetw.rds") %>% rename(withbetw=WithinBetween_interactions) %>% dplyr::select(study,cells, withbetw)
soc_df <- left_join(soci, socd, by=c("study","cells"))
soc_df %>% group_by(study) %>% dplyr::summarise(min=min(direct, na.rm=T),
                                               max=max(direct, na.rm=T),
                                               n=n()) %>% 
  arrange(desc(n)) %>% View()


#soc_1 <- left_join(soci, socd, by=c("study","cells"))

# convert landcover to factor + regroup
soc_df <- soc_df %>% mutate(lc = case_when(landcover < 12 ~ "human",
                                                      landcover >= 12 & landcover < 18 ~ "agriculture",
                                                      landcover == 18 ~ "pasture",
                                                      landcover > 18 & landcover <= 22 ~ "agriculture",
                                                      landcover == 23 ~ "broad_leaved",
                                                      landcover == 24 ~ "coniferous",
                                                      landcover == 25 ~ "mixed",
                                                      landcover > 25 & landcover <= 34 ~ "natural_veg",
                                                      landcover > 34 ~ "water_marshes"), 
                            lc = factor(lc))
# create dummy variables from landcov
soc_df %>% 
  to_dummy(lc, suffix = "label") %>% 
  bind_cols(soc_df) -> soc_df 
nrow(soc_df)
# Collinearity ----
library(usdm)
library(stringi)


vif_df <- NULL
for (i in unique(soc_df$study)) {
  
  #i <- "Alb"
  soc1 <- soc_df %>% filter(study == i) %>% 
    
    dplyr::select(distance_forest,distance_edges,distance_urban,# the covariates
    distance_pasture,distance_crops, distance_river,human_print, enet_density,
    tcd, distance_roads,distance_paths, elevation,slope ,hli, twi)
  soc2 <- soc1[complete.cases(soc1),]
  soc2 <- soc2 %>% select_if(~ !any(is.na(.)))
  if (nrow(soc2) < 1) {next}
  #plot(grid_scale)
  vif_i <- vif(soc2)
  vif_i$study <- i
  vif_df <- bind_rows(vif_df, vif_i)
}

vif_df %>%  
  filter_all(all_vars(!is.infinite(.))) %>%  group_by(Variables) %>% 
  dplyr::summarise(vif_mean=median(VIF, na.rm=T),
                   vif_sd=sd(VIF, na.rm=T))

vif2 <- vifstep(soc2, th=10) 
sel3 <- dropLayer(grid_all, c("dist_sd","dist_median","Degree","freq_mean","n_ind_abs","n_ind_max_abs"))

# in consequence, remove distance_urban and distance_crops


cov_belgium <- stack("data/derived_data/dellicour_area.grd")
cov_belgium <- subset(cov_belgium,c("distance_forest","distance_edges","distance_urban",
                      "tcd", "distance_roads",
                      "distance_paths", "elevation", "slope" ,"hli", "twi"))
cov_belgium <- scale(cov_belgium) 

#cov_belgium <- approxNA(cov_belgium, rule = 2)
#table(is.na(values(cov_belgium)))
#plot(cov_belgium,colNA="red")

# loop across study areas -----

sample_size <- 600
m_res_df <- NULL
var_imp_df <- NULL
aoa_df <- NULL

# for (i in c("Geneva Basin_27", "Grimsoe_46","HainichNP_15", "Hertogenwald_11",
#             "kaszo_22","Koberg_47","kostelec_30" ,"LageK_26_1",  "lpp_23",
#             "Marche_14", "MT_35_3" ,"MT_35_4" , "Oasi_Arezzo_21","Sardinia_9","StHubert_17",
#             "sumava_31","tutrakan_12", "Wur_40")) {# # { c("Alb_39","Alt_41"))
for (i in unique(soc_df$study)) { 
  print(i)
  #i <- "HainichNP"
  names(soc_df)
  dat <- soc_df %>% filter(study == i) %>%
    #filter(season == "warm" & tod == "night") %>% 
    dplyr::select(x, y, direct, indirect, #withbetw, # the response
                  distance_forest,distance_edges,#distance_urban,# the covariates
                  distance_pasture,#distance_crops, 
                  distance_river, human_print, enet_density,
                  tcd, distance_roads,distance_paths, elevation,
                  slope ,hli, twi)
  names(dat)
  dat <- dat[complete.cases(dat[, 5:17]),]
  preProcValues <- preProcess(dat[,-c(1:4)], method = c("center", "scale")) # add "nzv", "YeoJohnson"
  dat <- predict(preProcValues, dat) 
  
  # create indices based on variogram model
  vario <- variogram(direct~1, data=dat, locations= ~x+y)
  vario_fit <-fit.variogram(vario, vgm("Sph"))
  #plot(vario)
  direct <- rasterFromXYZ(dat[,1:3]) 
  #plot(direct)
  
  dat_sf <- st_as_sf(dat,coords = c("x", "y"), crs = 32632 )

  vario_grid <- st_make_grid(dat_sf, cellsize = vario_fit$range[2],
                             crs = 32632) %>% st_sf() %>% mutate(index = row_number())
  vario_sf <- st_intersection(dat_sf, vario_grid) 
  
  dat_df <- sf_to_df(vario_sf, fill = T) %>% dplyr::select(-c(sfg_id, point_id))
 

  trainDat <- dat_df %>%
    sample_n(150) # to reduce processing time in the test phase
  
  predictors <- c("distance_forest","distance_edges",
                  "tcd", "distance_roads",
                  "distance_paths", "elevation", "slope" ,"hli", "twi",
                  "distance_pasture", "distance_river",
                  "human_print", "enet_density")
  
  m_res_all <- NULL
  var_imp_all <- NULL
  aoa_all <- NULL

  for (j in c("direct", "indirect", "withbetw")) {
    #j <- "indirect"
    response <- "direct"
    #if (response == "withbetw") {trainDat = na.omit(trainDat)}
    if(nrow(trainDat) < 50) {next}
    model <- NULL
    
    indices <- CreateSpacetimeFolds(trainDat,spacevar = "index")
    
    unregister_dopar <- function() {
      env <- foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
    }
    unregister_dopar()
    ffsmodel_LLO <- ffs(trainDat[,predictors],
                        trainDat[,response],
                        metric="Rsquared",
                        method="rf", #tuneGrid=data.frame("mtry"=2),
                        verbose=FALSE,#ntree=50,
                        trControl=trainControl(method="cv",
                                               index = indices$index))
   
    ffsmodel_LLO
    ffsmodel_LLO$selectedvars
    plot_ffs(ffsmodel_LLO)
    predictors_sp <- stack("D:/PROJECTS/xx_GIS/data/derived_data/stacks_2022/HainichNP.grd")
    prediction_ffs <- predict(predictors_sp,ffsmodel_LLO)
    spplot(prediction_ffs)
    
    
    model_LLO <- train(trainDat[,predictors],
                       trainDat[,response],
                       method="rf",
                       importance=TRUE,
                       trControl=trainControl(method="cv",
                                              index = indices$index))
    model_LLO$finalModel
    prediction <- predict(predictors_sp,model_LLO)
    spplot(prediction)
    plot(varImp(model_LLO))
    plot_ffs(ffsmodel_LLO)
    plot(model_LLO)
    AOA <- aoa(predictors_sp,model)
    
    spplot(prediction,main="prediction for the AOA \n(spatial CV error applied)")+
      spplot(AOA$AOA,col.regions=c("grey","transparent"))
    
    library(doParallel)
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)
    
    model$finalModel
    summary(model$finalModel)$r.squared
    model <- caret::train(trainDat[,predictors],
                   trainDat[,response],
                   method="rf",
                   #tuneLength=3,
                   importance=TRUE,
                   #na.action = na.omit,
                   #trControl=trainControl(method="LOOCV")
                   trControl = trainControl(method = "LOOCV",
                                            verboseIter = TRUE,
                                            allowParallel = TRUE)
                   )
    plot(model)
    prediction <- predict(predictors_sp,model)
    spplot(prediction)
    plot(varImp(model))
    # unregister_dopar <- function() {
    #   env <- foreach:::.foreachGlobals
    #   rm(list=ls(name=env), pos=env)
    # }
    
    stopCluster(cl)
   
    m_res <- model$results
    m_res$study <- i
    m_res$interaction <- j
    m_res_all <- rbind(m_res_all, m_res)
    
    var_imp <- data.frame(varImp(model,scale = F)$importance) %>% tibble::rownames_to_column(var = "cov")
    var_imp$study <- i
    var_imp$interaction <- j
    var_imp_all <- rbind(var_imp_all, var_imp)

    # predict to belgium + extract AOA
    predict <- raster::predict(cov_belgium, model)
    names(predict) <- paste("predictions",j, i, sep ="_" )
    AOA <- NULL
    AOA <- try(CAST::aoa(cov_belgium, model))
    #plot(AOA)
    if (class(AOA)=="try-error") {next}
    
    pred_aoa <- stack(predict, AOA)
    writeRaster(pred_aoa, paste("data/derived_data/belgium/",paste("predict", i, j, sep="_"),".grd", sep = ""), overwrite = T)
    
    areaAOA <- sum(AOA$AOA[] >= 1, na.rm = T)
    areaNAOA <- sum(AOA$AOA[] < 1, na.rm = T)
    aoa_i <- data.frame(areaAOA, areaNAOA)
    aoa_i$study <- i
    aoa_i$interaction <- j
    aoa_all <- rbind(aoa_all, aoa_i)
    model <- NULL
  }
    aoa_df <- rbind(aoa_df, aoa_all)
    m_res_df <- rbind(m_res_df, m_res_all)
    var_imp_df <- rbind(var_imp_df, var_imp_all)
}
  plot(pred_aoa)
write.csv(aoa_df,"tables/belgium/aoa_df.csv")
write.csv(m_res_df,"tables/belgium/m_res_df.csv")
write.csv(var_imp_df,"tables/belgium/var_imp_df.csv")
library(egg)
# variable importance
ggplot(var_imp_df, aes(x = cov, y = Overall)) + 
  geom_boxplot() + coord_flip() +
  facet_grid(.~interaction) +
  theme_article()

# aoa
ggplot(aoa_df,aes(x = areaAOA, y = reorder(study,areaAOA) )) +
  geom_histogram(stat = "identity") + 
  facet_grid(.~interaction) +
  theme_article()
 
  predict <- raster::predict(cov_belgium, model)

  plot(predict)
  AOA <- aoa(cov_belgium, model)
  plot(AOA$AOA)
  grid.arrange(spplot(AOA$DI,col.regions=viridis(100),main="DI"),
               spplot(predict, col.regions=viridis(100),
                      main="prediction for AOA \n(LOOCV error applies)")+ spplot(AOA$AOA,col.regions=c("grey","transparent")),ncol=2)
?aoa
  areaAOA <- sum(AOA$AOA[] >= 1, na.rm = T)
  areaNAOA <- sum(AOA$AOA[] < 1, na.rm = T)

  # Prediction
  cov <- stack(paste("D:/PROJECTS/xx_GIS/data/derived_data/",Swabian Alps (Alb),  ".grd", sep = "")
  cov <- cov_alb[[c("distance_forest","distance_edges","distance_urban",
                    "tcd","forest_type", "distance_roads",
                    "distance_paths", "elevation", "slope" ,"hli", "twi")]]
  observed <- dat %>% dplyr::select(x, y, direct) %>%  rasterFromXYZ() 
  observed2 <- resample(observed, cov)
  
  prediction <- predict(cov,model)
  
}



  #mutate(scaled_direct = rescale(direct,, to=c(1,10)),
  #       scaled_indirect = rescale(indirect, to=c(1,10)),
  #       scaled_wb = rescale(withbetw,to=c(1,10))) 



# response distribution
#png("D:/PROJECTS/04_SocialScape/figures/responses_distribution.png",width = 960, height = 480)
#par(mfrow=c(1,3))
hist(dat$direct, main="Direct_Interactions", xlab="") # poisson-like
hist(dat$indirect, main="Indirect_Interactions",xlab="") # poisson-like
hist(dat$withbetw, main="WithinBetween_Interactions",xlab="")
#dev.off()

# data preparation ----

trainDat <- dat %>%
  sample_n(100) # to reduce processing time in the test phase

predictors <- c("distance_forest","distance_edges","distance_urban",
                "tcd", "distance_roads",
                "distance_paths", "elevation", "slope" ,"hli", "twi")


print(model)
plot(varImp(model,scale = F),col="black")

# Predict
cov_alb <- stack("D:/PROJECTS/xx_GIS/data/derived_data/Swabian Alps (Alb).grd")
cov <- cov_alb[[c("distance_forest","distance_edges","distance_urban",
                  "tcd","forest_type", "distance_roads",
                  "distance_paths", "elevation", "slope" ,"hli", "twi")]]
observed <- dat %>% dplyr::select(x, y, direct) %>%  rasterFromXYZ() 
observed2 <- resample(observed, cov)

prediction <- predict(cov,model)
plot(prediction)
truediff <- abs(prediction-observed2)


spplot(stack(prediction,truediff),main=c("prediction","reference"))

AOA <- aoa(cov, model)
attributes(AOA)$aoa_stats
grid.arrange(
  spplot(truediff,col.regions=viridis(100),main="true prediction error"),
  spplot(AOA$DI,col.regions=viridis(100),main="DI"),
  spplot(prediction, col.regions=viridis(100),main="prediction for AOA")+ spplot(AOA$AOA,col.regions=c("grey","transparent")), ncol=3)





# I random forest -----


library('randomForest')
library('pROC') 
dfr <- rasterFromXYZ(dat)  #Convert first two columns as lon-lat and third as value  
plot(dfr)
# predictors
predictors <- dfr[[c("distance_forest","distance_edges","distance_urban",
                  "tcd","forest_type", "distance_roads",
                  "distance_paths", "elevation", "slope" ,"hli", "twi")]]

spplot(stretch(predictors,0,1),col.regions=viridis(100))


# Data splitting ----

datTrain <- dat %>% sample_n(200)
#datTest  <- dat[-trainIndex,] %>% sample_n(100)


predictors <- c("distance_forest","distance_edges","distance_urban",
                  "tcd", "distance_roads",
                  "distance_paths", "elevation", "slope" ,"hli", "twi")
response <- "direct"


print(model)
print(model_indirect)
print(model_wb)
plot(varImp(model_direct,scale = F),col="black")
plot(varImp(model_indirect,scale = F),col="black")
plot(varImp(model_wb,scale = F),col="black")



predict_rf <- raster::predict(cov, model)
predict_glm <- mask(predict_glm, ac_ras)
plot(predict_rf,
     col=viridis(100), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 )

AOA <- aoa(cov, model)

plot(AOA)
grid.arrange(spplot(AOA$DI,col.regions=viridis(100),main="DI with sampling locations (red)"),
             spplot(predict_rf, col.regions=viridis(100),main="prediction for AOA \n(LOOCV error applies)")+ spplot(AOA$AOA,col.regions=c("grey","transparent")),ncol=2)



# Model training
# based on Meyer: https://cran.r-project.org/web/packages/CAST/vignettes/AOA-tutorial.html#generate-predictors 
allDat <- dat[,c("clusters",
                 'distance_forest',"distance_edges", "distance_urban",
                 "tcd","forest_type","distance_roads","distance_paths","elevation","slope","hli","twi")] %>% 
  mutate(clusters=factor(clusters))
allDat <- allDat[complete.cases(allDat),]
set.seed(1)
inTrainingSet <- createDataPartition(allDat$clusters, p=.75, list=F)
trainDat <- Dat[inTrainingSet]
testDat <- Dat[-inTrainingSet]
# Train the model
model_random <- train(allDat[,c(
  'distance_forest',"distance_edges", "distance_urban",
  "tcd","forest_type","distance_roads",
  "distance_paths","elevation","slope","hli","twi")],
  allDat$clusters,
  method="rf",
  importance=TRUE,
  preProc = c("center", "scale"), ## Center and scale the predictors for the training
  trControl = trainControl(method="cv"),
  metric='Accuracy')

# fit a random forest model (using ranger)
# split into training and testing
set.seed(23489)


model_rf <- train(trainDat[,c('distance_forest',"distance_edges", "distance_urban",
                              "tcd","forest_type","distance_roads", "distance_paths","elevation","slope","hli","twi")],
                  trainDat$clusters,
                  method="rf",
                  importance=TRUE,
                  preProc = c("center", "scale"), ## Center and scale the predictors for the training
                  trControl = trainControl(method="cv"),
                  metric='Accuracy')
print(model_rf)
pROC::auc(pROC::roc(trainDat$clusters, fitted(model_rf)))
rf_fit <- train(clusters ~ ., 
                data = trainDat, 
                method = "ranger")
rf_fit
# predict the outcome on a test set
rf_pred <- predict(model_rf, testDat)
# compare predicted outcome and true outcome
confusionMatrix(rf_pred, testDat$clusters)


print(model_random)
plot(varImp(model_random,scale = F),col="black")
varImp(model_random, scale = FALSE)


predictors <- predictors[[c('distance_forest',"distance_edges", "distance_urban",
                            "tcd","forest_type","distance_roads",
                            "distance_paths","elevation","slope","hli","twi")]]
plot(predictors)
prediction <- predict(predictors,model_random)
par(mfrow=c(1,2))
plot(clusters, # what to plot
     col=viridis(6), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 # horizontal size of legend
)
plot(prediction, # what to plot
     col=viridis(6), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 # horizontal size of legend
)
plot(AOA, # what to plot
     col=viridis(10), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 # horizontal size of legend
)

AOA <- aoa(predictors, model_random)
AOA
plot(AOA)
attributes(AOA)$aoa_stats

grid.arrange(
  spplot(clusters,col.regions=viridis(100),main="Clusters"),
  spplot(prediction,col.regions=viridis(100),main="Model prediction"),
  spplot(AOA$DI,col.regions=viridis(100),main="Dissimilarity (DI)"),
  spplot(prediction, col.regions=viridis(100),main="prediction for AOA")+ spplot(AOA$AOA,col.regions=c("grey","transparent")), ncol=2, nrow=2)

# AOA for spatially clustered data? --> my case! see https://cran.r-project.org/web/packages/CAST/vignettes/AOA-tutorial.html
##  --> cross-validation should be based on a leave-cluster-out approach, and the AOA estimation based on distances to a nearest data point not located in the same spatial cluster.
plot(intensity)
rc <- clump(intensity, directions=4) 
freq(rc)
plot(r)

# cluster ID # from https://stackoverflow.com/questions/24465627/clump-raster-values-depending-on-class-attribute
clVal <- unique(intensity)
clstrID_df <- NULL
# loop over all unique class values
for (i in clVal) {
  
  # create & fill in class raster
  r.class <- setValues(raster(intensity), NA)
  r.class[intensity == i]<- 1
  
  # clump class raster
  clp <- clump(r.class)
  
  clstrID <- as.data.frame(clp, na.rm=T)  %>%  add_rownames("cells") %>% mutate(cells=as.numeric(cells)) %>% collect() %>%  data.frame()
  clstrID$class <- i
  clstrID_df <- rbind(clstrID_df, clstrID)
  
} 
clstrID_df <- clstrID_df %>% mutate(clstrID=paste(clumps, class, sep="_")) %>% arrange(cells)
range(clstrID_df$cells )
range(dat$cells)
library(tidyverse)
trainDat <- dat %>% arrange(cells)
trainDat <- trainDat[complete.cases(trainDat),]
trainDat <- left_join(trainDat, clstrID_df, by="cells")


folds <- CreateSpacetimeFolds(trainDat, spacevar="clstrID",k=10)
set.seed(15)
model <- train(trainDat[,c('distance_forest',"distance_edges", "distance_urban",
                           "tcd","dominant_leaf","forest_type","small_woody_forest","distance_roads",
                           "distance_paths","elevation","slope","hli","twi")],
               trainDat$raster_use_intensity,
               method="rf",
               importance=TRUE,
               tuneGrid = expand.grid(mtry = c(2:length(c('distance_forest',"distance_edges", "distance_urban",
                                                          "tcd","dominant_leaf","forest_type","small_woody_forest","distance_roads",
                                                          "distance_paths","elevation","slope","hli","twi")))),
               trControl = trainControl(method="cv",index=folds$index))

dataset <- trainDat %>% mutate(intensity=factor(intensity), landcover=factor(landcover)) 

control <- trainControl(method='cv', 
                        number=10)#, 

# Comparison prediction error with model error
model_rf$results
# Relationship between the DI and the performance measure
AOA_calib <- calibrate_aoa(AOA_spatial,model,window.size = 200,length.out =200, multiCV=TRUE,showPlot=FALSE)
AOA_calib$plot
spplot(AOA_calib$AOA$expected_RMSE,col.regions=viridis(100),main="expected RMSE")+
  spplot(AOA$AOA,col.regions=c("grey","transparent"))

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3,
                        search='grid')
#Metric compare model is Accuracy
dataset <- trainDat %>% mutate(intensity=factor(raster_use_intensity)) %>% dplyr::select(-raster_use_intensity)

metric <- "Accuracy"
set.seed(123)
library(randomForest)
install.packages("mlbench")
library(mlbench)
library(caret)
library(e1071)
data(Sonar)
str(Sonar)
#Number randomely variable selected is mtry
#mtry <- sqrt(ncol(x))
#tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(intensity~., 
                    data=dataset, 
                    method='rf', 
                    metric='Accuracy', 
                    importance=TRUE,
                    preProc = c("center", "scale"), 
                    #tuneGrid=tunegrid, 
                    trControl=control)
print(rf_default)
plot(rf_default)
plot(varImp(rf_default,scale = F),col="black")


f1 <- direct~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi






# 1 lm 
lm_direct <- lm(f1,data = dat)

# 2 lme 
lme_direct <- lme(f1,
                  random = ~ 1 | study,
                  data = dat)
# 3 glm 
glm_direct <- glm(f1, data = dat, family = "gaussian")

# residuals inspection
shapiro.test(resid(lm_direct)[1:500])
shapiro.test(resid(lme_direct)[1:500])            
shapiro.test(resid(glm_direct)[1:500])

library(broom.mixed)
df_lm <- augment(lm_direct)
df_lme <- augment(lme_direct)
df_glm <- augment(glm_direct)
ggplot(df_glm, aes(x = .fitted, y = .resid)) + geom_point()


# with autocovariance components

# 4 lme+ ----
write.csv(dat, "D:/PROJECTS/04_SocialScape/data/dat_alb.csv")
direct_Ratio <- lme(direct~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi,
                    random = ~ 1 | study,
                    correlation = corRatio(1, form = ~ x + y),
                    data = dat)
Vario_lme1 <- nlme::Variogram(lme_direct, form=~ x+y,
                             robust=T, maxDist=5000, resType="pearson")
Vario_lme2 <- nlme::Variogram(lme_direct, form=~ x+y,
                              robust=T, maxDist=5000, resType="normalized")
Vario_ratio1 <- nlme::Variogram(direct_Ratio, form=~ x+y,
                              robust=T, maxDist=5000, resType="pearson")
Vario_ratio2 <- nlme::Variogram(direct_Ratio, form=~ x+y,
                              robust=T, maxDist=5000, resType="normalized")
plot(Vario_lme_d_Ratio1E, smooth=F)
plot(Vario_lme2, smooth=F)
# compare with ratio corStructure

# 5 glm+ -----
coords <- as.matrix(cbind(dat$x, dat$y) )
ac <- autocov_dist(dat$direct, coords, nbs = 500, longlat =F,zero.policy=F,type="one")
dat$ac <- ac
ac_box <- autocov_dist(dat$box_direct, coords, nbs = 500, longlat =F,zero.policy=F,type="one")
dat$ac_box <- ac_box
f2 <- direct~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi+ac
glm_auto <- glm(f2, data = dat, family = 'gaussian')
f3 <- box_direct~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi+ac
glm_box <- glm(f3, data = dat, family = 'gaussian')

saveRDS(glm_auto, "tables/model_auto.rds")
summary(glm_auto)
summary(glm_box)
AIC(glm_box, glm_auto, glm_direct)

shapiro.test(glm_auto$residuals[1:100])
shapiro.test(glm_box$residuals[1:100])
png("D:/PROJECTS/04_SocialScape/figures/resid_distrib.png",width = 960, height = 480)
par(mfrow=c(1,3))
hist(glm_direct$residuals, main="GLM", xlab="residuals")
hist(glm_auto$residuals, main="GLM + autocovariate", xlab="residuals")
hist(glm_box$residuals, main="GLM + autocovariate + boxcox", xlab="residuals")
dev.off()

#extract residuals of spatial model
res2 <- glm_box$residuals
res2 <- data.frame(Residuals = res2, x = dat$x, y = dat$y)
dat_sp <- SpatialPointsDataFrame(cbind(dat$x, dat$y), dat)
lstw  <- nb2listw(knn2nb(knearneigh(dat_sp, k = 100)))
moran.test(residuals.glm(glm_box), lstw) 

df_box <- augment(glm_box) %>% dplyr::select(.fitted,.resid, box_direct) %>% mutate(model="box")
df_auto <- augment(glm_auto) %>% dplyr::select(.fitted,.resid, direct) %>% mutate(model="auto")
df_glm <- augment(glm_direct) %>% dplyr::select(.fitted,.resid, direct) %>% mutate(model="glm")
df <- bind_rows(df_auto, df_glm, df_box)
ggplot(df, aes(x = .fitted, y = .resid, group=model, color=model)) + geom_point(alpha=.5) + egg::theme_article()
#ggplot(df_glm, aes(x = .fitted, y = .resid)) + geom_point() 
ggsave("D:/PROJECTS/04_SocialScape/figures/resid_glm_auto.png")

# 6 spamm ----
# based on https://www.r-bloggers.com/2019/09/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/
## fit non-spatial model

spamm_d <- glm(direct ~ distance_forest+distance_edges+distance_urban+tcd+
                  distance_roads+distance_paths+elevation+slope+hli+twi, data=dat)
# plot residuals
dat$resid_glm <- resid(spamm_d)
# dat$resid_std <- rescale(dat$resid, 1, 10)
ggplot(dat, aes(x = x, y = y)) +
  geom_point() +
  scale_size_continuous(range = c(1,10))
# formal test
library(spaMM)
library(glmmTMB)
library("NLMR")
library(DHARMa)
library(ROI.plugin.glpk)

sims <- simulateResiduals(spamm_d)
testSpatialAutocorrelation(sims, x = dat$x, y = dat$y, plot = FALSE)

# fit the model
datsub <- dat %>% sample_n(500) #%>% na.omit()
spamm_direct <- fitme(direct ~ #landcov +
                   distance_forest+distance_edges+distance_urban+
                   tcd+distance_roads+
                   distance_paths+elevation+slope+hli+twi +
                   Matern(1 | x + y),
                 data = dat) # this take a bit of time

spamm_box <- fitme(box_direct ~ #landcov +
                        distance_forest+distance_edges+distance_urban+
                        tcd+distance_roads+
                        distance_paths+elevation+slope+hli+twi +
                        Matern(1 | x + y),
                      data = dat)
saveRDS(spamm_box, "tables/spamm_box.rds")
# to test 
glm_direct <- glm(f1, data = dat, family = gaussian(link=log))
# model summary
summary(spamm_direct)
summary(spamm_box)

dd <- dist(datsub[,c("x","y")])
mm <- MaternCorr(dd, nu = 2.72641461, rho =0.01181105)
plot(as.numeric(dd), as.numeric(mm), xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation")
sims <- simulateResiduals(spamm_box)
?simulateResiduals
plot(sims)
par(mfrow=c(1,2))
str(m_spamm)
?fitme







# Predict ----
cov_alb <- stack("D:/PROJECTS/xx_GIS/data/derived_data/Swabian Alps (Alb).grd")
cov <- cov_alb[[c("distance_forest","distance_edges","distance_urban",
           "tcd","forest_type", "distance_roads",
           "distance_paths", "elevation", "slope" ,"hli", "twi")]]
predict_glm <- raster::predict(cov, glm_direct, type="response")
predict_glm <- mask(predict_glm, ac_ras)
plot(predict_glm,
     col=viridis(100), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 )


# for glm_auto
# rasterize ac
ac <- data.frame(ac = ac, x = dat$x, y = dat$y)
coordinates(ac) <- ~ x + y
ac_ras <- rasterize(ac, cov, field = 'ac')
plot(ac_ras)
names(ac_ras) <- "ac"
cov_ac <- stack(cov, ac_ras)
predict_auto <- raster::predict(cov_ac, glm_auto, type="response")

png("D:/PROJECTS/04_SocialScape/figures/predict_glm_auto.png",width = 960, height = 480)
par(mfrow=c(1,2))
plot(predict_glm,
     col=viridis(100), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2,
     box=FALSE, axes=F,
     main="model without")
plot(predict_auto,
     col=viridis(100), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2,
     box=FALSE, axes=F,
     main="model with autocovariate")
dev.off()

# tidymodels ----
# norm_rec <- 
norm_rec <- recipe(direct + indirect + withbetw ~ distance_forest+distance_edges+distance_urban+tcd+#forest_type+
         distance_roads+distance_paths+elevation+slope+hli+twi, data = dat) %>%
  step_normalize(everything())
set.seed(57343)
folds <- vfold_cv(dat, repeats = 10)

folds <- 
  folds %>%
  mutate(recipes = map(splits, prepper, recipe = norm_rec))

library(pls)

get_var_explained <- function(recipe, ...) {
  
  # Extract the predictors and outcomes into their own matrices
  y_mat <- bake(recipe, new_data = NULL, composition = "matrix", all_outcomes())
  x_mat <- bake(recipe, new_data = NULL, composition = "matrix", all_predictors())
  
  # The pls package prefers the data in a data frame where the outcome
  # and predictors are in _matrices_. To make sure this is formatted
  # properly, use the `I()` function to inhibit `data.frame()` from making
  # all the individual columns. `pls_format` should have two columns.
  pls_format <- data.frame(
    endpoints = I(y_mat),
    measurements = I(x_mat)
  )
  # Fit the model
  mod <- plsr(endpoints ~ measurements, data = pls_format)
  
  # Get the proportion of the predictor variance that is explained
  # by the model for different number of components. 
  xve <- explvar(mod)/100 
  
  # To do the same for the outcome, it is more complex. This code 
  # was extracted from pls:::summary.mvr. 
  explained <- 
    drop(pls::R2(mod, estimate = "train", intercept = FALSE)$val) %>% 
    # transpose so that components are in rows
    t() %>% 
    as_tibble() %>%
    # Add the predictor proportions
    mutate(predictors = cumsum(xve) %>% as.vector(),
           components = seq_along(xve)) %>%
    # Put into a tidy format that is tall
    pivot_longer(
      cols = c(-components),
      names_to = "source",
      values_to = "proportion"
    )
}

folds <- 
  folds %>%
  mutate(var = map(recipes, get_var_explained),
         var = unname(var))

variance_data <- 
  bind_rows(folds[["var"]]) %>%
  filter(components <= 15) %>%
  group_by(components, source) %>%
  summarize(proportion = mean(proportion))

ggplot(variance_data, aes(x = components, y = proportion, col = source)) + 
  geom_line() + 
  geom_point() 

pls.model = plsr(direct + indirect + withbetw ~ distance_forest+distance_edges+distance_urban+tcd+#forest_type+
                   distance_roads+distance_paths+elevation+slope+hli+twi, data = dat, validation = "CV")
cv = RMSEP(pls.model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
# Rerun the model
pls.model = plsr(direct + indirect + withbetw ~ distance_forest+distance_edges+distance_urban+tcd+#forest_type+
                   distance_roads+distance_paths+elevation+slope+hli+twi, data = dat, ncomp = best.dims)
coefficients = coef(pls.model)
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
coefficients = sort(coefficients[, 1 , 1])
barplot(tail(coefficients, 10))
devtools::install_github("SValv/displayR")
install.packages("displayR")
library(displayR)
names(coefficients) = TidyLabels(Labels(dat)[-1])
coefficients = sort(coefficients, decreasing = TRUE)

# For multiple study areas
data <- soc_indirect %>% filter(study %in% c("Alb_39", "HainichNP_15"))
regressions <-
data %>%
  nest(data = c(-study)) %>% 
  mutate(
    fit = map(data, ~ lm(layer ~ distance_forest+distance_edges+distance_urban+tcd+#forest_type+
                           distance_roads+distance_paths+elevation+slope+hli+twi, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  ) %>% 
  #unnest(tidied) %>% 
  select(-data, -fit)

regressions %>% 
  select(tidied) %>% 
  unnest(tidied)

regressions %>% 
  select(glanced) %>% 
  unnest(glanced)

# Variance estimation
set.seed(27)
boots <- bootstraps(data, times = 1000, apparent = TRUE)
boots

fit_lm_on_bootstrap <- function(split) {
  lm(layer ~ distance_forest+distance_edges+distance_urban+tcd+#forest_type+
       distance_roads+distance_paths+elevation+slope+hli+twi, analysis(split))
}

boot_models <-
  boots %>% 
  mutate(model = map(splits, fit_lm_on_bootstrap),
         coef_info = map(model, tidy))

boot_coefs <- 
  boot_models %>% 
  unnest(coef_info)
percentile_intervals <- int_pctl(boot_models, coef_info)
percentile_intervals

ggplot(boot_coefs, aes(estimate)) +
  geom_histogram(bins = 30) +
  facet_wrap( ~ term, scales = "free") +
  geom_vline(aes(xintercept = .lower), data = percentile_intervals, col = "blue") +
  geom_vline(aes(xintercept = .upper), data = percentile_intervals, col = "blue")

boot_aug <- 
  boot_models %>% 
  sample_n(200) %>% 
  mutate(augmented = map(model, augment)) %>% 
  unnest(augmented)

boot_aug


# with Zuur book ----

E <- rstandard(lm_direct)
mydata <- data.frame(E,train_data$x, train_data$y)
library(sp)
coordinates(mydata) <- ~ train_data.x + train_data.y
bubble(mydata, "E", col=c("black", "grey"), main="residuals")
glsd <- gls(
  model=direct ~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi,
  data=train_data )
Vario.gls <- nlme::Variogram(glsd, form=~ x+y,
          robust=T, maxDist=10000, resType="pearson")
plot(Vario.gls, smooth=T)
?variogram
Vario1 <- gstat::variogram(E ~ 1, data=mydata)
Vario2 <- gstat::variogram(E ~ 1, data=mydata, alpha=c(0, 45, 90, 135))
plot(Vario2)
# Test multiple correlation structure 
B1A <- gls(model=direct ~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi,
           correlation=corSpher(form=~ x+y, nugget=T),
           data=train_data )
B1B <- gls(model=direct ~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi,
           correlation=corLin(form=~ x+y, nugget=T),
           data=train_data )
B1C <- gls(model=direct ~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi,
           correlation=corRatio(form=~ x+y, nugget=T),
           data=train_data )
B1D <- gls(model=direct ~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi,
           correlation=corGaus(form=~ x+y, nugget=T),
           data=train_data )
B1E <- gls(model=direct ~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi,
           correlation=corExp(form=~ x+y, nugget=T),
           data=train_data )
AIC(glsd);AIC(B1A);AIC(B1B);AIC(B1C);AIC(B1D);AIC(B1E)
anova(glsd, B1A)



# significance? from https://hughst.github.io/week-6/
#confint(m_spamm,"pred") # from https://kimura.univ-montp2.fr/~rousset/spaMM/spaMMintro.pdf
coefs <- as.data.frame(summary(m_spamm)$beta_table)
row <- row.names(coefs) %in% c('distance_forest',"distance_edges", "distance_urban",
                               "tcd","dominant_leaf","forest_type","small_woody_forest","distance_roads",
                               "distance_paths","elevation","slope","hli","twi")
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, 'Cond. SE']
data.frame(variables=c('distance_forest',"distance_edges", "distance_urban",
                       "tcd","dominant_leaf","forest_type","small_woody_forest","distance_roads",
                       "distance_paths","elevation","slope","hli","twi"),
           lower=lower, upper=upper)
# plotting effects
library(pdp)
pdep_effects(m_spamm,"elevation")
plot_effects2(m_spamm,focal_var="tcd")  # need source(functions.R)
plot_effects2(m_spamm,focal_var="elevation") 
points(corridor~tcd, data=dat, col="blue",pch=20)

# Mapping
filled.mapMM(lfit,add.map=TRUE,plot.axes=quote({axis(1);axis(2)}),
             decorations=quote(points(pred[,coordinates],pch=15,cex=0.3)),
             plot.title=title(main="Inferred prevalence, North Cameroon",
                              xlab="longitude",ylab="latitude"))



# Predictions 
# Create an empty raster with the same extent and resolution as the bioclimatic layers
latitude_raster <- longitude_raster <-raster(nrows = nrow(cov_alb),
                                             ncols = ncol(cov_alb),
                                             ext = extent(cov_alb))

# Change the values to be latitude and longitude respectively
longitude_raster[] <- coordinates(longitude_raster)[,1]
latitude_raster[] <- coordinates(latitude_raster)[,2]

# Now create a final prediction stack of the 4 variables we need
names(cov_alb)
pred_stack <- stack(cov_alb[[c("distance_forest","distance_edges","distance_urban",
                               "tcd","dominant_leaf","forest_type","small_woody_forest", "distance_roads",
                               "distance_paths", "elevation", "slope" ,"hli", "twi")]],
                    longitude_raster,
                    latitude_raster)

# Rename to ensure the names of the raster layers in the stack match those used in the model
names(pred_stack) <- c("distance_forest","distance_edges","distance_urban",
                       "tcd","dominant_leaf","forest_type","small_woody_forest", "distance_roads",
                       "distance_paths", "elevation", "slope" ,"hli", "twi", "x", "y")
plot(pred_stack)
# Make predictions using the stack of covariates and the model
predicted_prevalence_raster <- predict(pred_stack, m_spamm)
plot(predicted_prevalence_raster)
gps_sf <- st_as_sf(gps_utm, coords = c("x", "y"), crs = 32632, agr = "constant")
gps_sf2 = gps_sf %>%  st_sample(100000) %>% st_set_crs(32632) %>%   st_transform(3035) %>% as_Spatial()
points(gps_sf2, pch=".")
plot(soc_direct)

# Validation -
## v-folds cross-validation
# take 20% to act as validation set
set.seed(1)
validation_rows <- sample(1:nrow(dat), 120)
data_train <- dat[-validation_rows,]
data_valid <- dat[validation_rows,]

# Fit model using 80%
glm_mod_2_spatial_validation <- spaMM::fitme(social_direct ~ #landcov +
                                               distance_forest+distance_edges+distance_urban+
                                               tcd+dominant_leaf+forest_type+small_woody_forest+distance_roads+
                                               distance_paths+elevation+slope+hli+twi +
                                               Matern(1 | x + y),
                                             data = data_train, family = "binomial")
summary(glm_mod_2_spatial_validation)
# PRedict to validation rows and compare
predictions_validation <- predict(glm_mod_2_spatial_validation, data_valid)
ggplot() + geom_point(aes(as.vector(predictions_validation), data_valid$social_direct))

ETH_malaria_data <- read.csv("https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week1/Lab_files/Data/mal_data_eth_2009_no_dups.csv", header=T) # Case data
ETH_Adm_1 <- raster::getData("GADM", country="ETH", level=1) # Admin boundaries
# Calculate mse
??mse
library(ModelMetrics)
library(spaMM)
mse(predictions_validation, data_valid$social_direct)
# the effect of distance to forest
range(dat$x)
range(dat$y)
levels(dat$landcov)

newdat <- expand.grid(x = 5000, y = 5200, distance_forest= seq(0, 1000, length.out = 100), 
                      landcov = factor(1, levels=c("broad_leaved_forest", "non_irrigated_arable_land",
                                                   "coniferous_forest","natural_grasslands","pastures",
                                                   "fruit_trees_and_berry_plantations","transitional_woodland_shrub",
                                                   "mixed_forest")))

newdat$social_direct <- as.numeric(predict(m_spamm, newdat, re.form = NA)) # re.form = NA used to remove spatial effects
newdat$social_direct <- newdat$social_direct + mean(c(0,fixef(m_spamm)[2:8])) # to remove region effect
# get 95% confidence intervals around predictions
newdat <- cbind(newdat, get_intervals(m_spamm, newdata = newdat, intervals = "fixefVar", re.form = NA) + mean(c(0,fixef(m_spamm)[2:8])))


gg1 <- ggplot(dat, aes(x = distance_forest, y = social_direct)) +
  geom_point() +
  geom_path(data = newdat) +
  geom_ribbon(data = newdat, aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975), alpha = 0.2)

# now for landcover effect
newdat <- data.frame(x = 5000, y = 5200, distance_forest = mean(dat$distance_forest), 
                     landcov = c("broad_leaved_forest", "non_irrigated_arable_land",
                                 "coniferous_forest","natural_grasslands","pastures",
                                 "fruit_trees_and_berry_plantations","transitional_woodland_shrub",
                                 "mixed_forest")) # averaging out elevation effect
newdat$social_direct <- as.numeric(predict(m_spamm, newdat, re.form = NA))
# get 95% CI
newdat <- cbind(newdat,get_intervals(m_spamm, newdata = newdat, intervals = "fixefVar", re.form = NA))

gg2 <- ggplot(dat, aes(x = landcov, y = social_direct)) +
  geom_jitter() +
  geom_point(data = newdat, color = "red", size = 2) +
  geom_linerange(data = newdat, aes(ymin = fixefVar_0.025, ymax = fixefVar_0.975), color = "red")

# plot together
library(gridExtra)
grid.arrange(gg1, gg2, ncol = 2)

library(fields)
library(raster)
# derive a DEM
elev_m <- Tps(dat_samp[,c("x","y")], dat_samp$distance_forest)
r <- raster(xmn = 527343, xmx = 528343, ymn = 5361763, ymx = 5362763, resolution = 100)
elev <- raster::interpolate(r, elev_m)
plot(elev)
# for the region info use the limits given in ca20
pp <- SpatialPolygons(list(Polygons(list(Polygon(ca20$reg1)), ID = "reg1"),Polygons(list(Polygon(ca20$reg2)), ID = "reg2"), Polygons(list(Polygon(ca20$reg3)), ID = "reg3")))
region <- rasterize(pp, r)

# predict at any location
newdat <- expand.grid(x = seq(4960, 5960, length.out = 50), y = seq(4830, 5710, length.out = 50))
newdat$elevation <- extract(elev, newdat[,1:2])
newdat$region <- factor(extract(region, newdat[,1:2]))
# remove NAs
newdat <- na.omit(newdat)
# predict
newdat$calcium <- as.numeric(predict(m_spamm, newdat))

(gg_spamm <- ggplot(newdat,aes(x=x, y=y, fill = calcium)) +
    geom_raster() +
    scale_fill_viridis())




crs(stack_ras) <- CRS('+init=EPSG:32632')
CBS_bb <- bb(stack_ras)

# read Microsoft Bing satellite and OpenCycleMap OSM layers
CBS_osm1 <- read_osm(CBS_bb, type="bing")
CBS_osm2 <- read_osm(CBS_bb, type="opencyclemap")
class(CBS_osm1)
crs=32632
# plot OSM raster data
contour(stack_ras[[1]])
qtm(CBS_osm1)  + tm_raster(stack_ras[[1]], palette = pal20)
qtm(CBS_osm2)
pal20 <- c("#003200", "#3C9600", "#006E00")
tm_shape(land) +
  tm_raster("cover", palette = pal20, title = "Global Land Cover") + 
  tm_layout(scale=.8, legend.position = c("left","bottom"))


removeCollinearity(stack_ras, plot = T, multicollinearity.cutoff = 0.7)
vif(stack_ras)
writeRaster(stack_ras, paste("D:/PROJECTS/04_SocialScape/data/derived_data/",i), overwrite=T)
# II random forest -----
# based on https://topepo.github.io/caret/pre-processing.html
## Predictors Pre-processing
predictors_df <- dat[,c('landcover','distance_forest',"distance_edges", "distance_urban",
                     "tcd","forest_type","distance_roads",
                     "distance_paths","elevation","slope","hli","twi")]
str(allDat)
nzv <- nearZeroVar(predictors_df, saveMetrics= TRUE)
nzv[nzv$nzv,][1:10,]

# correlated predictors
descrCor <-  cor(predictors_df)
highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .999) # 0
highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
predictors <- predictors [,-highlyCorDescr]
descrCor2 <- cor(predictors_df )
summary(descrCor2[upper.tri(descrCor2)])

# Model training
# based on Meyer: https://cran.r-project.org/web/packages/CAST/vignettes/AOA-tutorial.html#generate-predictors 
allDat <- dat[,c("clusters",
                 'distance_forest',"distance_edges", "distance_urban",
                   "tcd","forest_type","distance_roads","distance_paths","elevation","slope","hli","twi")] %>% 
  mutate(clusters=factor(clusters))
allDat <- allDat[complete.cases(allDat),]
set.seed(1)
inTrainingSet <- createDataPartition(allDat$clusters, p=.75, list=F)
trainDat <- Dat[inTrainingSet]
testDat <- Dat[-inTrainingSet]
# Train the model
model_random <- train(allDat[,c(
  'distance_forest',"distance_edges", "distance_urban",
                                  "tcd","forest_type","distance_roads",
                                  "distance_paths","elevation","slope","hli","twi")],
                      allDat$clusters,
                      method="rf",
                      importance=TRUE,
                      preProc = c("center", "scale"), ## Center and scale the predictors for the training
                      trControl = trainControl(method="cv"),
                      metric='Accuracy')

# fit a random forest model (using ranger)
# split into training and testing
set.seed(23489)
train_index <- sample(1:nrow(allDat), 0.8 * nrow(allDat))
trainDat <- allDat[train_index, ]
testDat <- allDat[-train_index, ]

model_rf <- train(trainDat[,c('distance_forest',"distance_edges", "distance_urban",
                            "tcd","forest_type","distance_roads", "distance_paths","elevation","slope","hli","twi")],
                  trainDat$clusters,
                  method="rf",
                  importance=TRUE,
                  preProc = c("center", "scale"), ## Center and scale the predictors for the training
                  trControl = trainControl(method="cv"),
                  metric='Accuracy')
print(model_rf)
pROC::auc(pROC::roc(trainDat$clusters, fitted(model_rf)))
rf_fit <- train(clusters ~ ., 
                data = trainDat, 
                method = "ranger")
rf_fit
# predict the outcome on a test set
rf_pred <- predict(model_rf, testDat)
# compare predicted outcome and true outcome
confusionMatrix(rf_pred, testDat$clusters)


print(model_random)
plot(varImp(model_random,scale = F),col="black")
varImp(model_random, scale = FALSE)


predictors <- predictors[[c('distance_forest',"distance_edges", "distance_urban",
                             "tcd","forest_type","distance_roads",
                             "distance_paths","elevation","slope","hli","twi")]]
plot(predictors)
prediction <- predict(predictors,model_random)
par(mfrow=c(1,2))
plot(clusters, # what to plot
     col=viridis(6), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 # horizontal size of legend
)
plot(prediction, # what to plot
     col=viridis(6), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 # horizontal size of legend
)
plot(AOA, # what to plot
     col=viridis(10), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 # horizontal size of legend
)

AOA <- aoa(predictors, model_random)
AOA
plot(AOA)
attributes(AOA)$aoa_stats

grid.arrange(
  spplot(clusters,col.regions=viridis(100),main="Clusters"),
  spplot(prediction,col.regions=viridis(100),main="Model prediction"),
  spplot(AOA$DI,col.regions=viridis(100),main="Dissimilarity (DI)"),
  spplot(prediction, col.regions=viridis(100),main="prediction for AOA")+ spplot(AOA$AOA,col.regions=c("grey","transparent")), ncol=2, nrow=2)

# AOA for spatially clustered data? --> my case! see https://cran.r-project.org/web/packages/CAST/vignettes/AOA-tutorial.html
##  --> cross-validation should be based on a leave-cluster-out approach, and the AOA estimation based on distances to a nearest data point not located in the same spatial cluster.
plot(intensity)
rc <- clump(intensity, directions=4) 
freq(rc)
plot(r)

# cluster ID # from https://stackoverflow.com/questions/24465627/clump-raster-values-depending-on-class-attribute
clVal <- unique(intensity)
clstrID_df <- NULL
# loop over all unique class values
for (i in clVal) {
  
  # create & fill in class raster
  r.class <- setValues(raster(intensity), NA)
  r.class[intensity == i]<- 1
  
  # clump class raster
  clp <- clump(r.class)
  
  clstrID <- as.data.frame(clp, na.rm=T)  %>%  add_rownames("cells") %>% mutate(cells=as.numeric(cells)) %>% collect() %>%  data.frame()
  clstrID$class <- i
  clstrID_df <- rbind(clstrID_df, clstrID)
  
} 
clstrID_df <- clstrID_df %>% mutate(clstrID=paste(clumps, class, sep="_")) %>% arrange(cells)
range(clstrID_df$cells )
range(dat$cells)
library(tidyverse)
trainDat <- dat %>% arrange(cells)
trainDat <- trainDat[complete.cases(trainDat),]
trainDat <- left_join(trainDat, clstrID_df, by="cells")


folds <- CreateSpacetimeFolds(trainDat, spacevar="clstrID",k=10)
set.seed(15)
model <- train(trainDat[,c('distance_forest',"distance_edges", "distance_urban",
                           "tcd","dominant_leaf","forest_type","small_woody_forest","distance_roads",
                           "distance_paths","elevation","slope","hli","twi")],
               trainDat$raster_use_intensity,
               method="rf",
               importance=TRUE,
               tuneGrid = expand.grid(mtry = c(2:length(c('distance_forest',"distance_edges", "distance_urban",
                                                          "tcd","dominant_leaf","forest_type","small_woody_forest","distance_roads",
                                                          "distance_paths","elevation","slope","hli","twi")))),
               trControl = trainControl(method="cv",index=folds$index))

dataset <- trainDat %>% mutate(intensity=factor(intensity), landcover=factor(landcover)) 

control <- trainControl(method='cv', 
                        number=10)#, 

# Comparison prediction error with model error
model_rf$results
# Relationship between the DI and the performance measure
AOA_calib <- calibrate_aoa(AOA_spatial,model,window.size = 200,length.out =200, multiCV=TRUE,showPlot=FALSE)
AOA_calib$plot
spplot(AOA_calib$AOA$expected_RMSE,col.regions=viridis(100),main="expected RMSE")+
  spplot(AOA$AOA,col.regions=c("grey","transparent"))

control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3,
                        search='grid')
#Metric compare model is Accuracy
dataset <- trainDat %>% mutate(intensity=factor(raster_use_intensity)) %>% dplyr::select(-raster_use_intensity)

metric <- "Accuracy"
set.seed(123)
library(randomForest)
install.packages("mlbench")
library(mlbench)
library(caret)
library(e1071)
data(Sonar)
str(Sonar)
#Number randomely variable selected is mtry
#mtry <- sqrt(ncol(x))
#tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(intensity~., 
                    data=dataset, 
                    method='rf', 
                    metric='Accuracy', 
                    importance=TRUE,
                    preProc = c("center", "scale"), 
                    #tuneGrid=tunegrid, 
                    trControl=control)
print(rf_default)
plot(rf_default)
plot(varImp(rf_default,scale = F),col="black")



# III mlr ----
# based on https://geocompr.robinlovelace.net/spatial-cv.html

dat <- dat %>% mutate(landcover=factor(landcover), intensity=factor(intensity))
fit <- glm(intensity ~ distance_forest+distance_edges+distance_urban+tcd+dominant_leaf+forest_type+
             distance_roads+distance_paths+elevation+slope+hli+twi, data=dat, family="binomial")
summary(fit)
pred_glm = predict(object = fit, type = "response")

# making the prediction
ta <- cov_alb[[c("distance_forest","distance_edges","distance_urban","tcd","dominant_leaf",
               "forest_type", "distance_roads","distance_paths", "elevation", "slope","hli","twi")]]
pred = raster::predict(ta, model = fit, type = "response")
plot(pred)
pROC::auc(pROC::roc(dat$intensity, fitted(fit)))

# mlr package
install.packages("mlr")
library(mlr)
# coordinates needed for the spatial partitioning
# select response and predictors to use in the modeling
data = dat %>% dplyr::select(-human_print, -small_woody_forest, -cells, -landcov, -landcover) %>% 
  sample_n(300)
# drawing random sample with a subjective minimum distance in order to possibly alleviate spatial autocorrelation of your model residuals 

coords = data[, c("x", "y")]
data <- data %>% dplyr::select(-x, -y)
str(data)
# create task
task = makeClassifTask(data = data, target = "intensity",
                       positive = 1, coordinates = coords)
# learner
learner.rf = makeLearner("classif.randomForest", predict.type = "prob", fix.factors.prediction = TRUE)
# Train the learner
resampling = makeResampleDesc("RepCV", fold = 5, reps = 5)

set.seed(123)
out = resample(learner = learner.rf, task = task,
               resampling = resampling, measures = list(auc))

mean(out$measures.test$auc, na.rm=T)

# OLD CODE ----

# from https://rpubs.com/rowlandw17/456209
dat <- soc_df %>% filter(study %in% c("Alb_39")) %>%
  mutate(scaled_direct = rescale(direct),
         scaled_indirect = rescale(indirect),
         scaled_wb = rescale(withbetw))# %>% 
 # sample_n(1000) 


# extract residuals
res <- m1$residuals
# extract coordinates corresponding to each residual and combine
res <- data.frame(Residuals = res, x = dat$x, y = dat$y)
#plot
plot <- ggplot(res, aes(x = x, y = y)) + geom_point(aes(colour = Residuals, size = Residuals)) + 
  geom_point(shape = 1, aes(size = Residuals, colour = sign) ,colour = "black") + 
  scale_colour_gradient2(low = "#E8F3EB", mid = "#FF1C55",
                         high = "#560000", midpoint = 8, space = "Lab",
                         na.value = "grey50", guide = "colourbar")
plot
# Convert data to spatial points dataframe
dat <- SpatialPointsDataFrame(cbind(dat$x, dat$y), dat)
library(spdep)
# Construct weights matrix in weights list form using the 10 nearest neighbors 
lstw  <- nb2listw(knn2nb(knearneigh(dat, k = 100)))
# Compute Moran's I using residuals of model
moran.test(residuals.glm(m1), lstw) # p-value < 0.05 therefore we can reject null hypothesis that there is no spatial autocorrelation in the residuals of the model.

# define cell coordinates 
coords <- as.matrix(cbind(dat$x, dat$y) )

# construct autcovariate - increase neigbourhood dist (nbs) by increments of 0.1 till no cells have zero neighbours

ac <- autocov_dist(dat$direct, coords, nbs = 450, longlat =F,zero.policy=T,type="one")
# combine with cell coordinates
#ac <- data.frame(ac = ac, x = dat$x, y = dat$y)

# convert to spatial points df
#coordinates(ac) <- ~ x + y

# rasterize and specify ac as classifying field
#ac_ras <- rasterize(ac, agr, field = 'ac')

dat$ac <- ac
f2 <- scaled_direct~ distance_forest+distance_edges+distance_urban+tcd+distance_roads+distance_paths+elevation+slope+hli+twi+ac
m2 <- glm(f2, data = dat, family = 'gaussian')
summary(m2)
AIC(m1, m2)

#extract residuals of spatial model
res2 <- m2$residuals

# add coords
res2 <- data.frame(Residuals = res2, x = dat$x, y = dat$y)

#plot2
plot2 <- ggplot(res2, aes(x = x, y = y)) + geom_point(aes(colour = Residuals, size = Residuals)) + 
  geom_point(shape = 1, aes(size = Residuals, colour = sign) ,colour = "black") + 
  scale_colour_gradient2(low = "#E8F3EB", mid = "#FF1C55",
                         high = "#560000", midpoint = 8, space = "Lab",
                         na.value = "grey50", guide = "colourbar")
plot2
# Compute Moran's I using residuals of updated model and weight list
moran.test(residuals.glm(m2), lstw) 

# tidymodels ----
# set new model
set_new_model("gls_tidy")
set_model_mode(model = "gls_tidy", mode = "regression")
set_model_engine(
  "gls_tidy", 
  mode = "regression", 
  eng = "gls"
)

set_dependency("gls_tidy", eng = "gls", pkg = "nlme")
set_model_arg(
  model = "gls_tidy",
  eng = "gls",
  parsnip = "correlation",
  original = "correlation",
  func = list(pkg = "nlme", fun = "bar"),
  has_submodel = FALSE
)
set_model_arg(
  model = "gls_tidy",
  eng = "gls",
  parsnip = "formula",
  original = "model",
  func = list(pkg = "nlme", fun = "bar"),
  has_submodel = FALSE
)
show_model_info("gls_tidy")

gls_tidy <-
  function(mode = "regression", model=NULL, correlation = NULL) {
    # Check for correct mode
    if (mode  != "regression") {
      rlang::abort("`mode` should be 'regression'")
    }
    
    # Capture the arguments in quosures
    args <- list(correlation = rlang::enquo(correlation))
    
    # Save some empty slots for future parts of the specification
    new_model_spec(
      "gls_tidy",
      args = args,
      eng_args = NULL,
      mode = mode,
      method = NULL,
      engine = NULL
    )
  }
set_fit(
  model = "gls_tidy",
  eng = "gls",
  mode = "regression",
  value = list(
    interface = "formula",
    protect = c("formula", "data"),
    func = c(pkg = "nlme", fun = "gls"),
    defaults = list()
  )
)
set_encoding(
  model = "gls_tidy",
  eng = "gls",
  mode = "regression",
  options = list(
    predictor_indicators = "traditional",
    compute_intercept = TRUE,
    remove_intercept = TRUE,
    allow_sparse_x = FALSE
  )
)

set_pred(
  model = "gls_tidy",
  eng = "gls",
  mode = "regression",
  type = "raw",
  value = list(
    pre = NULL,
    post = NULL,
    func = c(fun = "predict"),
    args = list(object = expr(object$fit), newdata = expr(new_data))
  )
)

set_pred(
  model = "gls_tidy",
  eng = "gls",
  mode = "regression",
  type = "numeric",
  value = list(
    pre = NULL,
    post = NULL,
    func = c(fun = "predict"),
    args =
      list(
        object = expr(object$fit),
        newdata = expr(new_data),
        type = "response"
      )
  )
)

show_model_info("gls_tidy")

gls_tidy(model=NULL,correlation = corGaus(1, form = ~ x + y)) %>%
  translate(engine = "gls")

intensity <- raster("data/derived_data/raster_use_intensity.tif")
corridor <- raster("data/derived_data/raster_corridor.tif")
soc_direct <- raster("data/derived_data/raster_social_direct.tif")
soc_indir <- raster("data/derived_data/raster_social_indirect.tif")
clusters <- raster("data/derived_data/clusters.tif")
plot(clusters)
names(clusters)
table(getValues(soc_direct))
soc_direct[soc_direct < 3] <- 0
soc_direct[soc_direct == 3] <- 1
names(soc_direct) <- "social_direct"

table(getValues(corridor))
corridor[corridor < 5] <- 0
corridor[corridor == 5] <- 1
names(corridor) <- "corridor"

table(getValues(intensity))
intensity[intensity < 3] <- 0
intensity[intensity == 3] <- 1
names(intensity) <- "intensity"

plot(corridor)
# import covariates rasterstack
cov_alb <- stack("D:/PROJECTS/xx_GIS/data/derived_data/ Swabian Alps (Alb).grd")
predictors <- resample(cov_alb, clusters, method="ngb")

plot(predictors)
df_legend <- read.csv("D:/PROJECTS/xx_GIS/data/copernicus/corine/legend_corine.csv") %>% dplyr::select(landcover=ID, landcov=landcover)

# extract
#cov_alb <- projectRaster(cov_alb, crs = crs(soc_direct))
p <- rasterToPoints(clusters)
p1 <- p[,c("x","y")]
pcov <- raster::extract(predictors, p1,cellnumbers=T)
dat <- data.frame(cbind(p, pcov))
dat <- left_join(dat, df_legend, by="landcover")
nrow(dat)
#dat_samp <- sample_n(dat, 100)
B# I Multinominal log regression ----
install.packages("mlogit")
library(mlogit)
?mlogit
rDat <- dat[,c("clusters", 
               'distance_forest',"distance_edges", "distance_urban",
               "tcd","forest_type","distance_roads","distance_paths","elevation","slope","hli","twi")] %>% 
  mutate(clusters=factor(clusters) )
mdata3 <- mlogit.data(rDat, varying=NULL, choice="clusters", shape='wide')
model_1 <- mlogit(clusters~1|distance_forest+distance_edges+distance_urban+tcd+forest_type+
                    distance_roads+distance_paths+elevation+slope+hli+twi, data=mdata3, reflevel="1")
summary(model_1)

# 0 BoxCox transformation ----
library(MASS)

# run a linear model
hist(dat$scaled_direct)
# run the box-cox transformation
bc_direct <- boxcox(scaled_direct ~ 1, data=dat)
bc_indirect <- boxcox(scaled_indirect ~ 1, data=dat)
bc_wb <- boxcox(scaled_wb ~ 1, data=dat)

(lambda_direct <- bc_direct$x[which.max(bc_direct$y)])
(lambda_indirect <- bc_indirect$x[which.max(bc_indirect$y)])
(lambda_wb <- bc_wb$x[which.max(bc_wb$y)])


powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
  
  boxcoxTrans <- function(x, lam1, lam2 = NULL) {
    
    # if we set lambda2 to zero, it becomes the one parameter transformation
    lam2 <- ifelse(is.null(lam2), 0, lam2)
    
    if (lam1 == 0L) {
      log(y + lam2)
    } else {
      (((y + lam2)^lam1) - 1) / lam1
    }
  }
  
  switch(method
         , boxcox = boxcoxTrans(y, lambda1, lambda2)
         , tukey = y^lambda1
  )
}


(trans_direct <- bc_direct$x[which.max(bc_direct$y)])
dat$box_direct <- (dat$scaled_direct^trans_direct-1)/trans_direct
mnew_direct <- lm(box_direct ~ 1, data=dat) #
hist(dat$direct)

(trans_indirect <- bc_indirect$x[which.max(bc_indirect$y)])
dat$box_indirect <- (dat$scaled_indirect^trans_indirect-1)/trans_indirect
mnew_indirect <- lm(box_indirect ~ 1, data=dat) #
hist(dat$box_indirect)

(trans_wb <- bc_wb$x[which.max(bc_wb$y)])
dat$box_wb <- (dat$scaled_wb^trans_wb-1)/trans_wb

mnew_wb <- lm(box_wb ~ 1, data=dat) #
hist(dat$scaled_wb)
shapiro.test(dat$box_direct[1:100])
# QQ-plot
op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m$residuals); qqline(m$residuals)
qqnorm(mnew_direct$residuals); qqline(mnew_direct$residuals)
qqnorm(mnew_indirect$residuals); qqline(mnew_indirect$residuals)
par(op)
names(dat)
dat <- dat %>% relocate(x, y,study, cells, direct, indirect, withbetw,
                        scaled_direct, scaled_indirect, scaled_wb, 
                        box_direct, box_indirect, box_wb)

#--------transform into binary data ?-------------------
#
#plot(r)
# direct <- rasterFromXYZ(dat[,1:3])  #Convert first two columns as lon-lat and third as value 
# indirect <- rasterFromXYZ(dat[,c(1:2,4)])  #Convert first two columns as lon-lat and third as value 
# direct <- sdmvspecies::rescale(direct)
# indirect <- sdmvspecies::rescale(indirect)
# 
# pp1 <- quantile(direct, c(0, 0.25, 0.75))
# ix <- findInterval(getValues(direct), pp1)
# classified_dir <- setValues(direct, ix)
# 
# pp2 <- quantile(indirect, c(0, 0.25, 0.75))
# ix2 <- findInterval(getValues(indirect), pp2)
# classified_ind <- setValues(indirect, ix2)
# 
# 
# #direct[ direct < .8] <- 0
# #direct[ direct >= .8] <- 1
# #indirect[ indirect < .8] <- 0
# #indirect[ indirect >= .8] <- 1
# #par(mfrow=c(1,2))
# #plot(classified_dir);plot(classified_ind)
# 
# diff <- overlay(classified_dir,classified_ind, fun=function(a,b) return(a==b))
# par(mfrow=c(1,3))
# plot(classified_dir);plot(classified_ind)
# plot(diff,
#      col=c('#FFE4E1','#228B22'),
#      legend=FALSE,
#      axes=FALSE)
# legend("top", legend=c("Agree", "Disagree"),
#        col=c("#228B22", "#FFE4E1"), pch = 15, cex=0.8)
#------------------------------------
