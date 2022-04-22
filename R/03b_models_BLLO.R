# spatially buffered leave-one-out validation approach ----
# based on http://portal.uni-freiburg.de/biometrie/mitarbeiter/dormann/roberts-et-al-2017-ecography.pdf
# and http://www.ecography.org/appendix/ecog-02881

# https://davidrroberts.wordpress.com/2016/03/11/spatial-leave-one-out-sloo-cross-validation/


# functions -----

# Spatially Buffered Leave-One-out Function
sloo.fun.lm <- function(dat, x, y, resp, ss, rad, modform){
  # dat = complete data file
  # x_coord(y_coord) = names of x(y) coordinates in data
  # truth = true value of response
  # ss = sample size to process (number of LOO runs)
  # rad = radius for the spatial buffer around test point
  # modform = formula for the GLM
  
  # Select the testing points
  test <- dat[sample(nrow(dat),ss),]
  #rad <- 2000
  #modform <- as.formula(direct ~ distance_forest + distance_edges +# the covariates
  #                        distance_pasture + distance_river + human_print + enet_density +
  #                        tcd + distance_roads + distance_paths + elevation + slope + hli + twi)
  var_all <- NULL
  # Vector of predictions
  for(i in 1:nrow(test)){
    if(i==1){p <- c()}
    # Training data (test point & buffer removed)
    #train <- dat[sqrt((dat[,"x"]-test[i,"x"])^2 +(dat[,"y"]-test[i,"y"])^2)>rad,]
    train <- dat[sqrt((dat[,x]-test[i,x])^2 +(dat[,y]-test[i,y])^2)>rad,]
    # Build the model
    m <- lm(modform, data=train)
    
    # Predict on test point
    p <- c(p, predict(m, test[i,], type="response"))
    
    var.coef <- data.frame(coef = summary(m)$coefficients[,4]) %>% rownames_to_column("variable") %>% 
      mutate(study = study, rep = i,rsq = summary(m)$r.squared)
    var_all <- rbind(var_all, var.coef)
    
  }
  

  # Residuals
  p.res <- test[,resp]-p
  # RMSE
  p.rmse <- sqrt(mean(p.res^2))
  
  list(SampleRows=as.numeric(rownames(test)), Truth=test[,resp], Predictions=p, Residuals=p.res,RMSE=p.rmse, rsq = r.sq,
       res_df=var_all)
  
  }

sloo.fun.rf <- function(dat, x, y, resp, ss, rad, modform, cov){
  # dat = complete data file
  # x_coord(y_coord) = names of x(y) coordinates in data
  # truth = true value of response
  # ss = sample size to process (number of LOO runs)
  # rad = radius for the spatial buffer around test point
  # modform = formula for the GLM
  # cov = covariates rastr stack where to predict
  
  # Select the testing points
  #ss <- 10
  test <- dat[sample(nrow(dat),ss),]
  # Vector of predictions
  prediction_all <- stack()
  var_all <- NULL
  for(i in 1:nrow(test)){
    if(i==1){p <- c()}
    #i <- 2
    # Training data (test point & buffer removed)
    train <- dat[sqrt((dat[,x]-test[i,x])^2 +(dat[,y]-test[i,y])^2)>rad,]
    # Build the model
    m <- randomForest(modform, data=train)
    #m <- randomForest(rf_form_direct, data=dat)
    #print(m)
    prediction <- predict(predictors_sp,m)
    prediction_all <- stack(prediction_all, prediction)
    
    var_exp <- tail(m$rsq,1) * 100
    var_all <- cbind(var_all,var_exp)
    
    #var_import <- data.frame(importance(m)) %>% rownames_to_column("variables")
    #var_import$SampleRows <- i
    
    #print(var_import)
    # Predict on test point
    p <- c(p, predict(m, test[i,], type="response"))
    
  }
  
  predict_mean <- calc(prediction_all, mean)
  writeRaster(predict_mean, paste0("data/derived_data/predictions/predict_",study, ".grd"), overwrite = T)
  #plot(predict_mean)
  # variable importance
  #var_import_df <- bind_rows(var_import_df, var_import)
  # Residuals
  p.res <- test[,resp]-p
  # RMSE
  p.rmse <- sqrt(mean(p.res^2))
  # variable importance
  
  list(SampleRows=as.numeric(rownames(test)), Truth=test[,resp], Predictions=p, Residuals=p.res,RMSE=p.rmse, var_exp = var_all)
  #print(var_import_df)
}

# SLOO CV run with subset of the whole data
# Sample Size = 10, Buffer Radius = 5
# Linear GLM


# INFERENTIAL APPROACH -----
lin_glm_between <- as.formula(used ~ distance_forest + distance_edges +# the covariates
                               distance_pasture + distance_river + #human_print + enet_density +
                               tcd + distance_roads + distance_paths + elevation + slope + hli + twi)

# Test on 1 study
lin.mod <- glm(lin_glm_between, data=dat, family = "binomial")


all_sloo <- NULL
var_df <- NULL
all_var <- NULL
for (i in c("Famenne", "Bialowieza", "fanel","kaszo","FCNP_west",
            "kostelec","Alb", "Wur","Grimsoe","HainichNP" )) {#unique(soc_df$study)) { 
  print(i)
  study <- i 
  #i <- "HainichNP"
  dat <- soc_df %>% filter(study %in% i) %>%
    dplyr::select(x, y, direct, indirect, #withbetw, # the response
                  distance_forest,distance_edges,#distance_urban,# the covariates
                  distance_pasture,#distance_crops, 
                  distance_river, #human_print,enet_density,
                  tcd, distance_roads,distance_paths, elevation,
                  slope ,hli, twi)
  #dat <- dat[complete.cases(dat[, 5:17]),]
  preProcValues <- preProcess(dat[,-c(1:4)], method = c("center", "scale")) # add "nzv", "YeoJohnson"
  dat <- predict(preProcValues, dat) 
  
  # create indices based on variogram model
  vario <- variogram(direct~1, data=dat, locations= ~x+y)
  vario_fit <- fit.variogram(vario, vgm("Sph"))
  # sloo
  max.size <- nrow(dat)
  buffer <- vario_fit$range[2]
  if (buffer > 10000){buffer <- 1000} 
  # Linear GLM
  ss.test.lin <- data.frame(N=1:max.size)
  test.out.lin <- try(sloo.fun.lm(dat, "x", "y", "indirect", max.size, buffer, lin_glm_direct))
  if(class(test.out.lin) == "try-error") {next}
  ss.test.lin[c("Truth","Pred","Res")] <- with(test.out.lin, cbind(Truth,Predictions,Residuals))
  
# Add the cumulative RMSE by simply calculating RMSE in sequence
# Add the random RMSE values by iteratively sampling from the data
  for(j in 1:max.size){
  # Random RMSE
  ss.test.lin[j,"Ran_RMSE"] <- with(ss.test.lin[sample(1:max.size,j),], sqrt(mean((Truth-Pred)^2)))
  # Cumulative RMSE
  ss.test.lin[j,"Cum_RMSE"] <- sqrt(mean((ss.test.lin$Res[1:j])^2))
  }
  ss.test.lin$study <- i
 all_sloo <- bind_rows(all_sloo, ss.test.lin)
 all_var <- rbind(all_var, test.out.lin$res_df) 
}
# Note that both RMSEs start low/high and stabilise to the same value
# The same happens to the quadratic model (not shown)
head(ss.test.lin); tail(ss.test.lin)

ggplot(all_sloo, aes(x = N, y = Ran_RMSE)) + 
  geom_point(alpha=0.1) +
  geom_line(aes(x = N, y = Cum_RMSE, colour = "red")) +
  facet_wrap(~study, scales = "free") + theme(legend.position = "none")
ggsave("doc/figures/supp_rmse.png", , unit = "px",width = 3000, height = 2500)

all_var %>% group_by(study, rep) %>% slice(1) %>% 
  ungroup() %>% group_by(study) %>% 
  dplyr::summarise(mean_rsq = mean(rsq))

all_var %>% group_by(variable) %>% 
  #dplyr::summarise(mean_coef = mean(coef)) %>% 
  ggplot(aes(x = reorder(variable,coef), y = coef)) + geom_boxplot() +
  ylim(c(0, 0.1)) + 
  geom_hline(yintercept = 0.05, color = "red1") +
  geom_hline(yintercept = 0.01, color = "red3") +
  geom_hline(yintercept = 0.001, color = "red4") +
  coord_flip() +
  facet_wrap(~study) 

p <- ggplot() + geom_vline(xintercept = 5) + coord_flip()
ggplotly(p)
# PREDICTIVE -----
# Test for a single study area
library(doParallel)
cl<-makeCluster(4)
registerDoParallel(cl)
library(randomForest)
sloo.fun.rf(dat, "x", "y", "direct", 10, 2500, rf_form)
model <- sloo.fun.rf(dat, "x", "y", "direct", max.size, buffer, rf_form)
stopCluster(cl)
predict <- raster::predict(cov, model)

# formulas
rf_form_direct <- as.formula(direct ~ distance_forest + distance_edges +# the covariates
                        distance_pasture + distance_river + human_print + enet_density +
                      tcd + distance_roads + distance_paths + elevation + slope + hli + twi)
rf_form_indirect <- as.formula(indirect ~ distance_forest + distance_edges +# the covariates
                               distance_pasture + distance_river + human_print + enet_density,
                             tcd + distance_roads + distance_paths + elevation + slope + hli + twi)
#rf_form_wb <- as.formula(withbetw ~ distance_forest + distance_edges +# the covariates
#                               distance_pasture + distance_river + human_print + enet_density +
#                             tcd + distance_roads + distance_paths + elevation + slope + hli + twi)

max.size <- 10
rf_direct <- NULL
rf_indirect <- NULL
rf_wb <- NULL
test_direct_df <- NULL
test_indirect_df <- NULL
test_wb_df <- NULL

for (i in unique(soc_df$study)) { #unique(soc_df$study)
  study <- "HainichNP"
  dat_i <- soc_df %>% filter(study == i)
  if (nrow(dat_i) < 250) {next}
  #i <- "HainichNP"
  #for (j in c("warm", "cold")) {
    #j <- "warm"
    #dat_j <- dat_i %>% filter(season == j) 
    #for (k in c("night", "day")){
      dat <- dat_i %>% # filter(tod == k) %>% 
        dplyr::select(x, y, direct, indirect, #withbetw, # the response
                      distance_forest,distance_edges,#distance_urban,# the covariates
                      distance_pasture,#distance_crops, 
                      distance_river, human_print,
                      enet_density,
                      tcd, distance_roads,distance_paths, elevation,
                      slope ,hli, twi) 
      if (nrow(dat) > 500) {dat <- dat %>% sample_n(500)}
      dat <- dat[complete.cases(dat[, 5:15]),]
      preProcValues <- preProcess(dat[,-c(1:4)], method = c("center", "scale")) # add "nzv", "YeoJohnson"
      dat <- predict(preProcValues, dat) 
      dat_omit = na.omit(dat)
      max.size <- nrow(dat)
      # create indices based on variogram model
      vario <- variogram(direct~1, data=dat, locations= ~x+y)
      vario_fit <- fit.variogram(vario, vgm("Sph"))
      # sloo
      
      buffer <- vario_fit$range[2]
      if (buffer > 10000){buffer <- 1000} 
      # randomForest
      #
      cl<-makeCluster(4)
      registerDoParallel(cl)
      test_direct_rf <- data.frame(N=1:max.size)
      test_indirect_rf <- data.frame(N=1:max.size)
      #test_wb_rf <- data.frame(N=1:max.size)
      
      max.size = 10
      cov <- stack(paste0("D:/PROJECTS/xx_GIS/data/derived_data/stacks_2022/", i,".grd"))
      test_out_direct <- try(sloo.fun.rf(dat, "x", "y", "direct", max.size, buffer, rf_form_direct, cov))
      test_out_indirect <- try(sloo.fun.rf(dat, "x", "y", "indirect", max.size, buffer, rf_form_indirect))
      #test_out_wb <- try(sloo.fun.rf(dat_omit, "x", "y", "withbetw", max.size, buffer, rf_form_wb))
      #+if(class(test_out_wb) == "try-error") {next}
      test_direct_rf[c("Truth","Pred","Res")] <- with(test_out_direct, cbind(Truth,Predictions,Residuals))
      test_indirect_rf[c("Truth","Pred","Res")] <- with(test_out_indirect, cbind(Truth,Predictions,Residuals))
      #test_wb_rf[c("Truth","Pred","Res")] <- with(test_out_wb, cbind(Truth,Predictions,Residuals))
      
      # Add the cumulative RMSE by simply calculating RMSE in sequence
      # Add the random RMSE values by iteratively sampling from the data
      for(l in 1:max.size){
        # Random RMSE
        test_direct_rf[l,"Ran_RMSE"] <- with(test_direct_rf[sample(1:max.size,l),], sqrt(mean((Truth-Pred)^2)))
        test_indirect_rf[l,"Ran_RMSE"] <- with(test_indirect_rf[sample(1:max.size,l),], sqrt(mean((Truth-Pred)^2)))
        #test_wb_rf[l,"Ran_RMSE"] <- with(test_wb_rf[sample(1:max.size,l),], sqrt(mean((Truth-Pred)^2)))
        # Cumulative RMSE
        test_direct_rf[l,"Cum_RMSE"] <- sqrt(mean((test_direct_rf$Res[1:l])^2))
        test_indirect_rf[l,"Cum_RMSE"] <- sqrt(mean((test_indirect_rf$Res[1:l])^2))
        #test_wb_rf[l,"Cum_RMSE"] <- sqrt(mean((test_wb_rf$Res[1:l])^2))
        
        test_direct_rf$study <- i
        #test_direct_rf$tod <- k
        #test_direct_rf$season <- j
        
        test_indirect_rf$study <- i
        #test_indirect_rf$tod <- k
        #test_indirect_rf$season <- j
        
        #test_wb_rf$study <- i
        #test_wb_rf$tod <- k
        #test_wb_rf$season <- j
        
      }
      test_direct_df <- bind_rows(test_direct_df, test_direct_rf)
      test_indirect_df <- bind_rows(test_indirect_df, test_indirect_rf)
      #test_wb_df <- bind_rows(test_wb_df, test_wb_rf)
      rf_direct <- bind_rows(rf_direct, test_direct_df)
      rf_indirect <- bind_rows(rf_indirect, test_indirect_df)
      
      stopCluster(cl)
      #rf_wb <- bind_rows(rf_wb, test_wb_df)
    }
    
    
  #}
  
 
#}

rf_direct %>% group_by(study) %>% tally() 
rf_wb %>% filter(study=="HainichNP" & season == "cold")

ggplot(rf_direct, aes(x = N, y = Ran_RMSE)) + 
  geom_point(alpha=0.1) +
  geom_line(aes(x = N, y = Cum_RMSE, colour = "red")) +
  facet_wrap(~study, scales = "free") + theme(legend.position = "none")
