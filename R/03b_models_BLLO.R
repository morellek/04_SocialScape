# spatially buffered leave-one-out validation approach ----
# based on http://portal.uni-freiburg.de/biometrie/mitarbeiter/dormann/roberts-et-al-2017-ecography.pdf
# and http://www.ecography.org/appendix/ecog-02881

# https://davidrroberts.wordpress.com/2016/03/11/spatial-leave-one-out-sloo-cross-validation/
dat
lin.glm.form <- as.formula(direct~distance_forest+distance_edges+distance_urban)
lin.mod <- glm(lin.glm.form, data=dat)
summary(lin.mod)
plot(predict(lin.mod),                                # Draw plot using Base R
     dat$direct,
     xlab = "Predicted Values",
     ylab = "Observed Values")
abline(a = 0,                                        # Add straight line
       b = 1,
       col = "red",
       lwd = 2)
data_mod <- data.frame(Predicted = predict(lin.mod),  # Create data for ggplot2
                       Observed = dat$direct)
ggplot(data_mod,                                     # Draw plot using ggplot2 package
       aes(x = Predicted,
           y = Observed)) +
  geom_point() +
  geom_abline(intercept = 0,
              slope = 1,
              color = "red",
              size = 2)

# Spatially Buffered Leave-One-out Function
sloo.fun.lm <- function(dat, x_coord, y_coord, resp, ss, rad, modform){
  # dat = complete data file
  # x_coord(y_coord) = names of x(y) coordinates in data
  # truth = true value of response
  # ss = sample size to process (number of LOO runs)
  # rad = radius for the spatial buffer around test point
  # modform = formula for the GLM
  
  # Select the testing points
  test <- dat[sample(nrow(dat),ss),]
  # Vector of predictions
  for(i in 1:nrow(test)){
    if(i==1){p <- c()}
    # Training data (test point & buffer removed)
    train <- dat[sqrt((dat[,x_coord]-test[i,x_coord])^2 +(dat[,y_coord]-test[i,y_coord])^2)>rad,]
    # Build the model
    m <- glm(modform, data=train)
    # Predict on test point
    p <- c(p, predict(m, test[i,], type="response"))
  }
  # Residuals
  p.res <- test[,resp]-p
  # RMSE
  p.rmse <- sqrt(mean(p.res^2))
  list(SampleRows=as.numeric(rownames(test)), Truth=test[,resp], Predictions=p, Residuals=p.res,RMSE=p.rmse)
}

sloo.fun.rf <- function(dat, x, y, resp, ss, rad, modform){
  # dat = complete data file
  # x_coord(y_coord) = names of x(y) coordinates in data
  # truth = true value of response
  # ss = sample size to process (number of LOO runs)
  # rad = radius for the spatial buffer around test point
  # modform = formula for the GLM
  
  # Select the testing points
  #ss <- 10
  test <- dat[sample(nrow(dat),ss),]
  # Vector of predictions
  for(i in 1:nrow(test)){
    if(i==1){p <- c()}
    #i <- 2
    # Training data (test point & buffer removed)
    train <- dat[sqrt((dat[,x]-test[i,x])^2 +(dat[,y]-test[i,y])^2)>rad,]
    # Build the model
    m <- randomForest(modform, data=train)
    print(m)
    # Predict on test point
    p <- c(p, predict(m, test[i,], type="response"))
  }
  # Residuals
  p.res <- test[,resp]-p
  # RMSE
  p.rmse <- sqrt(mean(p.res^2))
  list(SampleRows=as.numeric(rownames(test)), Truth=test[,resp], Predictions=p, Residuals=p.res,RMSE=p.rmse)
}

# SLOO CV run with subset of the whole data
# Sample Size = 10, Buffer Radius = 5
# Linear GLM
sloo.fun(dat, "x", "y", "direct", 100, 2500, lin.glm.form)
# rmse across study areas -----
all_sloo <- NULL
for (i in unique(soc_df$study)) { 
  print(i)
  #i <- "cab_24"
  dat <- soc_df %>% filter(study %in% i) %>%
    dplyr::select(x, y, direct, indirect, withbetw, # the response
                  distance_forest,distance_edges,distance_urban,# the covariates
                  tcd, distance_roads,distance_paths, elevation,
                  slope ,hli, twi)
  dat <- dat[complete.cases(dat[, 6:15]),]
  preProcValues <- preProcess(dat[,-c(1:2)], method = c("center", "scale")) # add "nzv", "YeoJohnson"
  dat <- predict(preProcValues, dat) 
  
  # create indices based on variogram model
  vario <- variogram(direct~1, data=dat, locations= ~x+y)
  vario_fit <- fit.variogram(vario, vgm("Sph"))
  # sloo
  max.size <- nrow(dat)
  buffer <- vario_fit$range[2]
  if (buffer > 10000){buffer <- 10000} 
  # Linear GLM
  ss.test.lin <- data.frame(N=1:max.size)
  test.out.lin <- try(sloo.fun(dat, "x", "y", "direct", max.size, buffer, lin.glm.form))
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
}
# Note that both RMSEs start low/high and stabilise to the same value
# The same happens to the quadratic model (not shown)
head(ss.test.lin); tail(ss.test.lin)

ggplot(all_sloo, aes(x = N, y = Ran_RMSE)) + 
  geom_point(alpha=0.1) +
  geom_line(aes(x = N, y = Cum_RMSE, colour = "red")) +
  facet_wrap(~study, scales = "free") + theme(legend.position = "none")
ggsave("doc/figures/supp_rmse.png", , unit = "px",width = 3000, height = 2500)


# randomForest across study areas -----
mod_rf <- randomForest(y=dat$y[train], x=data.used[train,c(3,4,7,8,9,10,11)], ntree=250, nodesize=10)
m <- randomForest(direct~.,data=dat[,c("direct","distance_forest","distance_edges")])
m$mse
m$rsq
m$predicted
varImpPlot(m)

iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE,
                        proximity=TRUE)
rf_form <- as.formula(direct~distance_forest+distance_edges+distance_urban)

# Test for a single study area
cl<-makeCluster(4)
registerDoParallel(cl)
sloo.fun.rf(dat, "x", "y", "direct", 10, 2500, rf_form)
model <- sloo.fun.rf(dat, "x", "y", "direct", max.size, buffer, rf_form)
stopCluster(cl)
predict <- raster::predict(cov, model)

max.size <- 10
rf_sloo <- NULL

for (i in unique(soc_df$study)) { 
  print(i)
  #i <- "cab_24"
  dat <- soc_df %>% filter(study %in% i) %>%
    dplyr::select(x, y, direct, indirect, withbetw, # the response
                  distance_forest,distance_edges,distance_urban,# the covariates
                  tcd, distance_roads,distance_paths, elevation,
                  slope ,hli, twi)
  dat <- dat[complete.cases(dat[, 6:15]),]
  preProcValues <- preProcess(dat[,-c(1:2)], method = c("center", "scale")) # add "nzv", "YeoJohnson"
  dat <- predict(preProcValues, dat) 
  
  # create indices based on variogram model
  vario <- variogram(direct~1, data=dat, locations= ~x+y)
  vario_fit <- fit.variogram(vario, vgm("Sph"))
  # sloo
  
  buffer <- vario_fit$range[2]
  if (buffer > 10000){buffer <- 10000} 
  # randomForest
  ss.test.lin <- data.frame(N=1:max.size)
  cl<-makeCluster(4)
  registerDoParallel(cl)
  test_out_rf <- try(sloo.fun.rf(dat, "x", "y", "direct", max.size, buffer, rf_form))
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
  rf_sloo <- bind_rows(rf_sloo, ss.test.lin)
}