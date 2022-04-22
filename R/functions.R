# depedencies 
library(adehabitatLT)
library(raster)
library(sf)
library(sp)
library(lubridate)
library(tidyverse)
library(pupilr)
library(data.table)
library(spatsoc)
library(wildlifeDI)
library(here)
#remotes::install_github("burchill/pupilr")

#gps_traj <- readRDS(file = "D:/PROJECTS/02_Contacts/materials/Bastille_Rousseau/gps_ltraj.rds")

# data requirement: ltraj object
social <- read.csv("tables/social_all.csv")

socialGrid <- function(data=data,  res=res, crs=crs, #window=day...,type=c('absolut', 'relative'),
                       contact_threshold=contact_threshold) {
  output <- NULL
  output_list <- NULL
  # Step 1: build up the grid (based on resolution argurment) ----
  #data <- gps_traj # to be deleted at the end
  #res <- 200 # to be deleted at the end
  #crs <- 32632 # to be deleted at the end
  # data <- gps_traj
  tt <-SpatialPoints(ld(data)[,1:2])
  tt1<-apply(coordinates(tt), 2, min)
  tt2<-apply(coordinates(tt), 2, max)
  ras<-raster(xmn=floor(tt1[1]), ymn=floor(tt1[2]),xmx=ceiling(tt2[1]), ymx=ceiling(tt2[2]), res=res)
  ras[] <- 1:ncell(ras)
  names(ras) <- 'grid_id'
  
  spdf_grid <- as(ras,'SpatialPolygonsDataFrame')
  sf_grid <- st_as_sf(spdf_grid)
  st_crs(sf_grid) = crs
  # plot(sf_grid$geometry)
  
  # Steps 2: tracking data conversion ----
  gps_sf <- st_as_sf(ld(data), coords = c("x", "y"), crs = crs, agr = "constant")
  gps_sf <- st_join(gps_sf, sf_grid, join = st_intersects)
  gps_df <- as.data.frame(gps_sf) %>% mutate(x = unlist(purrr::map(geometry,1)),
                                             y = unlist(purrr::map(geometry,2)),
                                             yr_month=format(as.Date(date), "%Y-%m"),
                                             yr_week=format(as.Date(date), "%Y-%W"),
                                             yr_day=format(as.Date(date), "%Y-%m-%d")) %>% 
    dplyr::select(id=burst, date, x, y, grid_id,yr_month, yr_week,  yr_day)
  # Steps 3: metrics calculation ----
  # 1 - tot_ind in the specific study area
  tot_ind <- gps_df %>% distinct(id) %>% nrow()
  # 2 - tracking_period by individuals
  track_period <- gps_df %>% group_by(id) %>% 
    dplyr::summarise(start=min(date),
                     end=max(date),
                     period=as.numeric(difftime(end, start, units="days")))
  # n_ind 
  # = total number of distinct ind having visited a particular grid
  n_ind <- gps_df %>% dplyr::select(grid_id,id, date) %>%
    group_by(grid_id) %>%
    dplyr::summarise(#n_ind_rel=(n_distinct(id)/tot_ind)*100,
                     n_ind_abs=n_distinct(id),.groups='keep') %>% arrange(grid_id)
  n_ind_empty <- data.frame(grid_id=getValues(ras), n_ind_abs=0) #, n_ind_rel=0
  n_ind_df <- left_join_and_overwrite(n_ind_empty, n_ind, by="grid_id")# %>% #str()
    #mutate(n_ind_rel=replace_na(n_ind_rel,0),
    #       n_ind_abs=replace_na(n_ind_abs,0))
  
  # 2 - n_ind_simult
  # = max number of individual observed in a particular grid
  # time_max = time at which it took place
  # add timegroup variable by grid
  DT <- data.table(gps_df)
  DT[, datetime := as.POSIXct(date, tz = 'UTC')]
  DT2 <- NULL
  #contact_threshold <- '5 minutes'
  for (i in unique(DT$grid_id)) {
    #i = 2
    subDT <- DT[grid_id == i]
    #if (length(unique(subDT$animals_original_id))==1) {next}
    group_times(DT = subDT, datetime = 'datetime', threshold = contact_threshold)
    #nrow(subDT)
    DT2 <- rbind(DT2, subDT)
  }
  
  #DT2 %>% dplyr::select(grid_id,id, date, timegroup) %>% arrange(timegroup) %>% pivot_wider(id_cols=grid_id)
  
  n_ind_simult <- DT2 %>% dplyr::select(grid_id,id, date, timegroup) %>%
    group_by(grid_id, timegroup) %>%
    dplyr::summarise(time_max=mean(date),
                     n_ind=n_distinct(id),.groups='keep') %>%
    group_by(grid_id) %>% dplyr::summarise(#n_ind_max_rel=(max(n_ind)/tot_ind)*100,
                                           n_ind_simult_max=max(n_ind),
                                           #n_ind_mean_rel=(mean(n_ind)/tot_ind)*100,
                                           n_ind_simult_mean=mean(n_ind))#%>% ggplot(.,aes(x=n_ind))  + geom_histogram()
  #n_ind_simult %>% ggplot(.,aes(x=n_ind_mean_rel))  + geom_histogram()
  n_ind_simult_empty <- data.frame(grid_id=getValues(ras), n_ind_simult_max=NA_real_,  n_ind_simult_mean=NA_real_) # n_ind_mean_rel=0, n_ind_max_rel=0
  n_ind_simult_df <- left_join_and_overwrite(n_ind_simult_empty, n_ind_simult, by="grid_id") #%>% #str()
    #mutate(n_ind_max_rel=replace_na(n_ind_max_rel,0),
    #       n_ind_max_abs=replace_na(n_ind_max_abs,0),
    #       n_ind_mean_rel=replace_na(n_ind_mean_rel,0),
    #       n_ind_mean_abs=replace_na(n_ind_mean_abs,0))
  
  # duration
  duration <- DT2 %>% dplyr::select(grid_id,id, date) %>%
    group_by(grid_id) %>% arrange(id, date) %>% data.table()
  
  duration[, visit_id := rleid(grid_id)]
  
  
  duration_df <- duration %>% group_by(id, grid_id, visit_id) %>% 
    dplyr::summarise(grid_id=min(grid_id),
                     grid_id2=max(grid_id),
                     start=min(date),
                     end=max(date),
                     duration=difftime(end, start, units="hours")) %>%
    group_by(grid_id) %>% dplyr::summarise(dur_mean=as.numeric(mean(duration)))
  
  duration_empty <- data.frame(grid_id=getValues(ras), dur_mean=NA_real_)
  duration_df <- left_join_and_overwrite(duration_empty, duration_df, by="grid_id")# %>% #str()
    #mutate(dur_mean=replace_na(dur_mean,0))
  
  # frequency 
  #duration %>% filter(id=="alb_7091") %>% arrange(grid_id, date) %>% View()
  freq <- duration %>% group_by(id, grid_id, visit_id) %>% tally() %>% 
    group_by(id, grid_id) %>% tally()
  freq_days <- left_join(freq, track_period, by="id") %>% 
    mutate(frequency=n/period) %>% 
    group_by(grid_id) %>% dplyr::summarise(freq_mean=mean(frequency))
  #freq_days %>% arrange(desc(freq_mean)) %>% ggplot(., aes(x=freq_mean)) + geom_histogram()
  frequency_empty <- data.frame(grid_id=getValues(ras), freq_mean=0)
  frequency_df <- left_join_and_overwrite(frequency_empty, freq_days, by="grid_id")# %>% #str()
  
  # distance
  distance_all <- NULL
  for (j in unique(DT2$grid_id)) {
   # i <- 725
    DT_edge <- DT2[grid_id == j]
    edges <- edge_dist(
      DT=DT_edge,
      threshold = 1000000,
      id = 'id',
      coords = c('x', 'y'),
      timegroup = 'timegroup',
      returnDist = TRUE,
      fillNA = F
  )
    distance <- edges %>% group_by(timegroup) %>%
      dplyr::summarise(mean_dist=mean(distance)) %>% 
      dplyr::summarise(dist_mean=mean(mean_dist),
                       dist_sd=sd(mean_dist),
                       dist_median=median(mean_dist))
    distance$grid_id <- j
    distance_all <- rbind(distance_all, distance)
  }
  distance_empty <- data.frame(grid_id=getValues(ras), dist_mean=NA_real_, dist_sd=NA_real_,dist_median=NA_real_)
  distance_df <- left_join_and_overwrite(distance_empty, distance_all, by="grid_id")# %>% #str()
  #distance_all %>% arrange(desc(mean_dist2)) %>% ggplot(., aes(x=mean_dist2)) + geom_histogram()
  
  # Group composition change: gcc
  gcc_all <- NULL
  for (k in unique(DT2$grid_id)) {
   # i <- 1027
   DT3 <-  DT2 %>% filter(grid_id==k) %>% arrange(timegroup)
   #View(DT3)
   cols <- c(grp_0 = NA_real_, grp_1 = NA_real_)
   gcc <- DT3 %>% dplyr::select(grid_id, timegroup,id, date) %>% group_by(timegroup) %>%  arrange(date) %>% 
    dplyr::summarise(date=mean(date),
                     grid_id=max(grid_id),
                     n=n_distinct(id),
                     IDs=list(sort(id))) %>% #View() # as_tibble()  %>%
    mutate(IDs2 = map2(IDs, lag(IDs), intersect),
           n2=lengths(IDs2),
           n1n2=n-n2,
           grp=case_when(n1n2 > 1~1,
                         TRUE~as.numeric(n1n2))) %>% #View()
     group_by(grp) %>% tally() %>% 
     pivot_wider(names_from=grp, values_from=n,names_prefix="grp_") %>% 
     add_column(!!!cols[!names(cols) %in% names(.)])
     
     gcc$grid_id <- k
     
     gcc_all <- rbind(gcc_all, gcc)
     #group_by(IDs) %>% tally() %>% View()
    
  }
  gcc_all <- gcc_all %>% group_by(grid_id) %>%  mutate(grp_change=grp_1*(grp_1/(grp_1+grp_0))) %>% dplyr::select(grid_id, grp_change)
  #gcc_all %>% arrange(grp_change) %>% View()
  gcc_empty <- data.frame(grid_id=getValues(ras), grp_change= NA_real_)
  gcc_df <- left_join_and_overwrite(gcc_empty, gcc_all, by="grid_id")
  
  # Social: within-Between interactions
  wb_all <- NULL
  for (l in unique(DT2$grid_id)) {
    # l <- 773
    DT3 <-  DT2 %>% filter(grid_id==l) %>% arrange(timegroup)
    if (nrow(DT3)==0) {next}
    #View(DT3)
    wb <- DT3 %>% dplyr::select(grid_id, timegroup,id, date) %>% group_by(timegroup) %>%  arrange(date) %>% #View()
      dplyr::summarise(date=mean(date),
                       grid_id=max(grid_id),
                       n=n_distinct(id),
                       IDs=list(sort(id))) %>% 
      filter(n>1)
    if (nrow(wb)==0) {next}
    for (m in 1:nrow(wb)) {
      wbj <- wb[m,]
      #j <- unique(wb$IDs)
      inds <- wbj$IDs[[1]]
      socj <- social %>% filter(ind1 %in% inds & ind2 %in% inds) %>% 
        filter(!duplicated(paste0(pmax(ind1, ind2), pmin(ind1, ind2)))) %>% 
        dplyr::summarise(meanSI=mean(si,na.rm=T),
                         meanPX=mean(prox,na.rm=T),
                         meanCS=mean(cs, na.rm=T),
                         meanHAI=mean(hai, na.rm=T),
                         meanCR=mean(cr, na.rm=T))
      
      wbj$si <- socj$meanSI
      wbj$px <- socj$meanPX
      wbj$cs <- socj$meanCS
      wbj$cr <- socj$meanCR
      wbj$hai <- socj$meanHAI
      wb_all <- rbind(wb_all, wbj)
    }
    
  }
  if(is.null(wb_all)) 
  {wb_all <- data.frame(grid_id=getValues(ras),si=NA_real_,  px=NA_real_,cs=NA_real_,  cr=NA_real_,  hai=NA_real_ ) }

  wb_sum <- wb_all %>%  group_by(grid_id) %>% dplyr::summarise(SI=median(si),
                                                               PX=median(px),
                                                               CS=median(cs),
                                                               CR=median(cr),
                                                               HAI=median(hai))
  wb_empty <- data.frame(grid_id=getValues(ras),SI=NA_real_,  PX=NA_real_,CS=NA_real_,  CR=NA_real_,  HAI=NA_real_ ) # n_ind_mean_rel=0, n_ind_max_rel=0
  wb_df <- left_join_and_overwrite(wb_empty, wb_sum, by="grid_id") #%>% #str()
  
  # Merge all metrics
  output <- left_join(n_ind_df, n_ind_simult_df, by="grid_id")
  output <- left_join(output, duration_df, by="grid_id")
  output <- left_join(output, frequency_df, by="grid_id")
  output <- left_join(output, distance_df, by="grid_id")
  output <- left_join(output, gcc_df, by="grid_id")
  output <- left_join(output, wb_df, by="grid_id")
  
  # create Raster stack
  r_stack=stack()
  for (l in 2:ncol(output)) {
    r_l <- setValues(ras, output[,l])
    names(r_l) <- names(output)[l]
    #ls_stack[[j]] = r_j
    r_stack <- stack( r_stack , r_l )
  }
  
  output_list <- list(output, r_stack, ras)
  return(output_list)
  
}


#test <- socialGrid(data=gps_traj, res=200, crs=32632, contact_threshold='15 minutes')
#plot(test[[2]])
#output <-test
# r_stack=stack()
# for (j in 2:ncol(output)) {
#   r_j <- setValues(ras, output[,j])
#   names(r_j) <- names(output)[j]
#   #ls_stack[[j]] = r_j
#   r_stack <- stack( r_stack , r_j )
# }

DT2 %>% group_by(grid_id) %>% tally()
wb_all <- NULL
for (i in unique(DT2$grid_id)) {
  # i <- 2651
  DT3 <-  DT2 %>% filter(grid_id==i) %>% arrange(timegroup)
  if (nrow(DT3)==0) {next}
  #View(DT3)
  wb <- DT3 %>% dplyr::select(grid_id, timegroup,id, date) %>% group_by(timegroup) %>%  arrange(date) %>% #View()
    dplyr::summarise(date=mean(date),
                     grid_id=max(grid_id),
                     n=n_distinct(id),
                     IDs=list(sort(id))) %>%
    filter(n>1)
  if (nrow(wb)==0) {next}
  for (j in 1:nrow(wb)) {
    wbj <- wb[j,]
    #j <- unique(wb$IDs)
    inds <- wbj$IDs[[1]]
    socj <- social %>% filter(ind1 %in% inds & ind2 %in% inds) %>%
      filter(!duplicated(paste0(pmax(ind1, ind2), pmin(ind1, ind2)))) %>%
      dplyr::summarise(meanSI=mean(si,na.rm=T),
                       meanPX=mean(prox,na.rm=T),
                       meanCS=mean(cs, na.rm=T),
                       meanHAI=mean(hai, na.rm=T),
                       meanCR=mean(cr, na.rm=T))

   wbj$si <- socj$meanSI
   wbj$px <- socj$meanPX
   wbj$cs <- socj$meanCS
   wbj$cr <- socj$meanCR
   wbj$hai <- socj$meanHAI
   wb_all <- rbind(wb_all, wbj)
  }
}
# wb_sum <- wb_all %>%  group_by(grid_id) %>% dplyr::summarise(SI=median(si),
#                                                              PX=median(px),
#                                                              CS=median(cs),
#                                                              CR=median(cr),
#                                                              HAI=median(hai))
# wb_empty <- data.frame(grid_id=getValues(ras),SI=NA_real_,  PX=NA_real_,CS=NA_real_,  CR=NA_real_,  HAI=NA_real_ ) # n_ind_mean_rel=0, n_ind_max_rel=0
# wb_df <- left_join_and_overwrite(wb_empty, wb_sum, by="grid_id") #%>% #str()
# 
# library(virtualspecies)
# output <- wb_df
# r_stack=stack()
# for (l in 2:ncol(wb_df)) {
#   r_l <- setValues(ras, output[,l])
#   names(r_l) <- names(output)[l]
#   #ls_stack[[j]] = r_j
#   r_stack <- stack( r_stack , r_l )
# }
# vif(r_stack)
# removeCollinearity(r_stack, plot = T, multicollinearity.cutoff = 0.5)


# Unsupervised classification ----
# raster_scale <- scale(r_stack, center=TRUE, scale=TRUE)
# plot(raster_scale)
# table_test <- test %>% dplyr::select(-c(grid_id,dist_sd,dist_median)) %>% na.omit()
# str(table_test)
# library(mclust)
# dataBIC <- mclustBIC(table_test[,c(1:7)])
# print(summary(dataBIC))
# plot(dataBIC)
# ls <- Mclust(table_test[,c(1:7)], G=5, modelNames="VVE")
# str(ls)
# ls[["parameters"]][["mean"]]
# #plot(ls, what = "uncertainty")
# ModPred <- predict.Mclust(ls, table_test[,c(1:7)]) # prediction
# Pred_ras <- ras # establishing a rediction raster
# values(Pred_ras) <- NA # set everything to NA
# # set values of prediction raster to corresponding classification according to rowname
# values(Pred_ras)[as.numeric(rownames(table_test))] <- as.vector(ModPred$classification)
# plot(Pred_ras)
# colours <- rainbow(ls$G) # define 7 colours
# plot(Pred_ras, # what to plot
#      col = colours, # colours for groups
#      colNA = "white", # which colour to assign to NA values
#      legend.shrink=1, # vertical size of legend
#      legend.width=2 # horizontal size of legend
# )
# table(values(Pred_ras))
# 
# 
# mod1dr <- MclustDR(ls)
# summary(mod1dr)
# plot(mod1dr, what = "pairs")
# plot(mod1dr, what = "boundaries", ngrid = 100)
# 
# # Supervised classification
# library(rpart)
# # Train the model
# cart <- rpart(as.factor(table_test$n_ind_abs)~., data=table_test, method = 'class', minsplit = 5)
# # print(model.class)
# # Plot the trained classification tree
# plot(cart, uniform=TRUE, main="Classification Tree")
# text(cart, cex = 0.8)
# 
# #Paremeters associated to each cluster (i.e. center)
# ls$parameters$mean
# #Proportion of each cluster (how frequent it is spatially)
# ls$parameters$pro
# 
# socialClust<-function(ls, max.n.clust=8) {
#   nvar<-nrow(ls$parameters$mean)
#   coef<-data.frame()
#   G<-ls$G
#   gg<-data.frame(cbind(t(ls$parameters$mean), ls$parameters$pro))
#   coef<-rbind(coef, gg)
#   coef[,-ncol(coef)] <- sapply(coef[-ncol(coef)],function(x) as.numeric(as.character(x)))
#   names(coef)[nvar+1]<-"Prop"
#   clust<-Mclust(coef[,1:nvar], G=1:max.n.clust)
#   coef$clust<-clust$classification
#   out<-list(clust, coef)
# }
# ls_clust <- socialClust(ls)
# 
# ls_clust[[1]]$parameters$mean
# 
# pop_stack<-function(clust_stack) {
#   n.clust<-length(names(clust_stack[[1]]))-1
#   ls<-lapply(clust_stack, function(x) x[[1]])
#   fct<-function(rast, val) {
#     values(rast) <- ifelse(values(rast) %in% c(val), 1, 0) # just replace
#     return(rast)
#   }
#   ls1<-list()
#   for (i in 1:n.clust) {ls1[[i]]<-lapply(ls, function(x) fct(x, i)) }
#   
#   #Mosaic all together
#   x <- ls1[[1]]
#   names(x)[1:2] <- c('x', 'y')
#   x$fun <- max
#   x$na.rm <- TRUE
#   y <- do.call(mosaic, x)
#   out<-stack(y)
#   for (i in 2:n.clust) {
#     x <- ls1[[i]]
#     names(x)[1:2] <- c('x', 'y')
#     x$fun <- max
#     x$na.rm <- TRUE
#     y <- do.call(mosaic, x)
#     out[[i]]<-y
#   }
#   names(out)<-paste("Clust", c(1:n.clust), sep="_")
#   return(out)
# }
# 
# 
# 
# socialClust<-function(table, max.n.clust=8, modelname="EEV", vars=c("n_ind_rel", "n_ind_abs","n_ind_max_rel",
#                                                                   "n_ind_max_abs","n_ind_mean_rel","n_ind_mean_abs",
#                                                                   "dur_mean")) {
#   
#   id<-unique(table$ID)
#   ls<-list()
#   for (i in 1:length(id)) {
#     i=1
#     tt<-scale(table[table$ID==id[i],vars])
#     try(ls[[i]]<-Mclust(tt, G=2:max.n.clust, modelNames=modelname))
#     print(id[i])
#   }
#   return(ls)
# }
# 
# 
# 
# 
# sfc_as_cols <- function(x, names = c("x","y")) {
#   stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
#   ret <- sf::st_coordinates(x)
#   ret <- tibble::as_tibble(ret)
#   stopifnot(length(names) == ncol(ret))
#   x <- x[ , !names(x) %in% names]
#   ret <- setNames(ret,names)
#   dplyr::bind_cols(x,ret)
# }
# 
# plot_effects2 <- function (object, focal_var, newdata = object$data, effects = NULL, 
#           xlab = focal_var, ylab = NULL, rgb.args = col2rgb("blue"), 
#           add = FALSE, ylim = NULL, ...) 
# {
#   if (is.null(effects)) 
#     effects <- pdep_effects(object, newdata = newdata, focal_var = focal_var, 
#                             indiv = FALSE, ...)
#   if (object$family$family == "binomial") {
#     resp <- object$y/object$BinomialDen
#     #if (is.null(ylab)) 
#     #  ylab <- paste("frequency(", formula.HLfit(object, 
#     #                                            which = "")[[2]][[2]], ")")
#   }
#   else {
#     resp <- object$y
#     if (is.null(ylab)) 
#       ylab <- paste(formula.HLfit(object, which = "")[[2]])
#   }
#   if (is.null(ylim)) 
#     ylim <- stats::quantile(resp, c(0.025, 0.975))
#   rgb.args <- as.list(rgb.args)
#   if (is.null(rgb.args$maxColorValue)) 
#     rgb.args$maxColorValue <- 255
#   colpts <- do.call(rgb, rgb.args)
#   rgb.args$alpha <- rgb.args$maxColorValue * 0.1
#   colshd <- do.call(rgb, rgb.args)
#   new.x <- effects[, "focal_var"]
#   if (!add) 
#     plot(NULL, type = "l", ylab = ylab, xlab = xlab, 
#          xlim = range(new.x), ylim = ylim)
#   polygon(c(new.x, rev(new.x)), c(effects[, "low"], rev(effects[, 
#                                                                 "up"])), border = NA, col = colshd)
#   focal_class <- class(newdata[, focal_var])
#   if ("numeric" %in% focal_class) {
#     lines(new.x, effects[, "pointp"], lwd = 2, col = colpts)
#   }
#   else points(new.x, effects[, "pointp"], col = colpts, 
#               pch = 19)
#   graphics::rug(newdata[, focal_var], side = 1, col = colpts)
#   graphics::rug(resp, side = 2, col = colpts)
#   invisible(effects)
# }

gls2 <- function (formula, data = sys.frame(sys.parent()), correlation = NULL, 
                  weights = NULL, subset, method = c("REML", "ML"), 
                  na.action = na.fail, control = list(), verbose = FALSE) 
{
  Call <- match.call()
  controlvals <- glsControl()
  if (!missing(control)) 
    controlvals[names(control)] <- control
  if (!inherits(formula, "formula") || length(formula) != 
      3L) {
    stop("\nformula must be a formula of the form \"resp ~ pred\"")
  }
  method <- match.arg(method)
  REML <- method == "REML"
  groups <- if (!is.null(correlation)) 
    getGroupsFormula(correlation)
  glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
  formula <- terms(formula, data = data)
  mfArgs <- list(formula = asOneFormula(formula(glsSt), formula, 
                                        groups), data = data, na.action = na.action)
  if (!missing(subset)) {
    mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2L]]
  }
  mfArgs$drop.unused.levels <- TRUE
  dataMod <- do.call(formula.frame, mfArgs)
  origOrder <- row.names(dataMod)
  if (!is.null(groups)) {
    groups <- eval(substitute(~1 | GR, list(GR = groups[[2L]])))
    grps <- getGroups(dataMod, groups, level = length(getGroupsFormula(groups, 
                                                                       asList = TRUE)))
    ord <- order(grps)
    grps <- grps[ord]
    dataMod <- dataMod[ord, , drop = FALSE]
    revOrder <- match(origOrder, row.names(dataMod))
  }
  else grps <- NULL
  X <- formula.frame(formula, dataMod)
  contr <- lapply(X, function(el) if (inherits(el, "factor")) 
    contrasts(el))
  contr <- contr[!unlist(lapply(contr, is.null))]
  X <- formula.matrix(formula, X)
  if (ncol(X) == 0L) 
    stop("no coefficients to fit")
  y <- eval(formula[[2L]], dataMod)
  N <- nrow(X)
  p <- ncol(X)
  parAssign <- attr(X, "assign")
  fTerms <- terms(as.formula(formula), data = data)
  namTerms <- attr(fTerms, "term.labels")
  if (attr(fTerms, "intercept") > 0) {
    namTerms <- c("(Intercept)", namTerms)
  }
  namTerms <- factor(parAssign, labels = namTerms)
  parAssign <- split(order(parAssign), namTerms)
  fixedSigma <- (controlvals$sigma > 0)
  attr(glsSt, "conLin") <- list(Xy = array(c(X, y), c(N, 
                                                      ncol(X) + 1L), list(row.names(dataMod), c(colnames(X), 
                                                                                                deparse(formula[[2]])))), dims = list(N = N, p = p, REML = as.integer(REML)), 
                                logLik = 0, sigma = controlvals$sigma, fixedSigma = fixedSigma)
  glsEstControl <- controlvals["singular.ok"]
  glsSt <- Initialize(glsSt, dataMod, glsEstControl)
  parMap <- attr(glsSt, "pmap")
  numIter <- numIter0 <- 0L
  repeat {
    oldPars <- c(attr(glsSt, "glsFit")[["beta"]], 
                 coef(glsSt))
    if (length(coef(glsSt))) {
      optRes <- if (controlvals$opt == "nlminb") {
        nlminb(c(coef(glsSt)), function(glsPars) -logLik(glsSt, 
                                                         glsPars), control = list(trace = controlvals$msVerbose, 
                                                                                  iter.max = controlvals$msMaxIter))
      }
      else {
        optim(c(coef(glsSt)), function(glsPars) -logLik(glsSt, 
                                                        glsPars), method = controlvals$optimMethod, 
              control = list(trace = controlvals$msVerbose, 
                             maxit = controlvals$msMaxIter, reltol = if (numIter == 
                                                                         0L) controlvals$msTol else 100 * .Machine$double.eps))
      }
      coef(glsSt) <- optRes$par
    }
    else {
      optRes <- list(convergence = 0)
    }
    attr(glsSt, "glsFit") <- glsEstimate(glsSt, control = glsEstControl)
    if (!needUpdate(glsSt)) {
      if (optRes$convergence) 
        stop(optRes$message)
      break
    }
    numIter <- numIter + 1L
    glsSt <- update(glsSt, dataMod)
    aConv <- c(attr(glsSt, "glsFit")[["beta"]], 
               coef(glsSt))
    conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
    aConv <- c(beta = max(conv[1:p]))
    conv <- conv[-(1:p)]
    for (i in names(glsSt)) {
      if (any(parMap[, i])) {
        aConv <- c(aConv, max(conv[parMap[, i]]))
        names(aConv)[length(aConv)] <- i
      }
    }
    if (verbose) {
      cat("\nIteration:", numIter)
      cat("\nObjective:", format(optRes$value), "\n")
      print(glsSt)
      cat("\nConvergence:\n")
      print(aConv)
    }
    if (max(aConv) <= controlvals$tolerance) {
      break
    }
    if (numIter > controlvals$maxIter) {
      stop("maximum number of iterations reached without convergence")
    }
  }
  glsFit <- attr(glsSt, "glsFit")
  namBeta <- names(glsFit$beta)
  attr(glsSt, "fixedSigma") <- fixedSigma
  attr(parAssign, "varBetaFact") <- varBeta <- glsFit$sigma * 
    glsFit$varBeta * sqrt((N - REML * p)/(N - p))
  varBeta <- crossprod(varBeta)
  dimnames(varBeta) <- list(namBeta, namBeta)
  Fitted <- fitted(glsSt)
  if (!is.null(grps)) {
    grps <- grps[revOrder]
    Fitted <- Fitted[revOrder]
    Resid <- y[revOrder] - Fitted
    attr(Resid, "std") <- glsFit$sigma/varWeights(glsSt)[revOrder]
  }
  else {
    Resid <- y - Fitted
    attr(Resid, "std") <- glsFit$sigma/varWeights(glsSt)
  }
  names(Resid) <- names(Fitted) <- origOrder
  apVar <- if (controlvals$apVar) 
    glsApVar(glsSt, glsFit$sigma, .relStep = controlvals[[".relStep"]], 
             minAbsPar = controlvals[["minAbsParApVar"]], 
             natural = controlvals[["natural"]])
  else "Approximate variance-covariance matrix not available"
  dims <- attr(glsSt, "conLin")[["dims"]]
  dims[["p"]] <- p
  attr(glsSt, "conLin") <- NULL
  attr(glsSt, "glsFit") <- NULL
  attr(glsSt, "fixedSigma") <- fixedSigma
  grpDta <- inherits(data, "groupedData")
  structure(class = "gls", list(formulaStruct = glsSt, 
                                dims = dims, contrasts = contr, coefficients = glsFit[["beta"]], 
                                varBeta = varBeta, sigma = if (fixedSigma) controlvals$sigma else glsFit$sigma, 
                                apVar = apVar, logLik = glsFit$logLik, numIter = if (needUpdate(glsSt)) numIter else numIter0, 
                                groups = grps, call = Call, method = method, fitted = Fitted, 
                                residuals = Resid, parAssign = parAssign, na.action = attr(dataMod, 
                                                                                           "na.action")), namBetaFull = colnames(X), units = if (grpDta) 
                                                                                             attr(data, "units"), labels = if (grpDta) 
                                                                                               attr(data, "labels"))
}
