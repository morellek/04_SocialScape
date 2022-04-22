library(here)
here()
library(moveNT)
source('R/functions_MoveNT.R')
source('R/functions.R')
library(adehabitatLT)
library(raster)
library(sp)
library(mclust)
library(igraph)
library(sp)
library(sf)
library(tidyverse)
library(amt)
library(spatsoc)
library(data.table)

gps_data <- readRDS("data/gps1h.rds")

gps_data %>% group_by(short_name2) %>% tally() %>% View()
# 2 Social grid ----
for (i in unique(gps_data$short_name2)) { #c("StHubert_17","Hertogenwald_11","Famenne_8", "Marche_14","Thuin_18")) {
  #i <- "HainichNP_15"
  gps_i <- gps_data %>% filter(short_name2==i)
  #for (j in c("cold", "warm")) {
    #j <- "warm"
  #  gps_season <- gps_i %>% filter(season==j)
  #  for (k in c("day", "night")) {
      # k <- "night"
  #    gps_tod <- gps_season %>% filter(tod_==k)
      gps_traj <- as.ltraj(xy = gps_i[,c("x","y")], date = gps_i$acquisition_time,
                           id = gps_i$animals_original_id, typeII = TRUE)
      socgrid <- socialGrid(data=gps_traj, res=200, crs=32632, contact_threshold='5 minutes')
      soc_tbl <- socgrid[[1]] %>% dplyr::select(-grid_id) %>%   janitor::remove_empty(which = "rows") %>% rownames_to_column("grid_id")
      soc_ras <- socgrid[[2]]
      soc_grid <- socgrid[[3]]
      
      crs(soc_ras) <- "+init=epsg:32632"
      crs(soc_grid) <- "+init=epsg:32632"
      i2 <- str_replace(i, "[.]", "_")
      write.csv(soc_tbl,paste0("tables/SocTbl/sampling1h/2022_03_30/soc_tbl_",i2, ".csv", sep=""))
      writeRaster(soc_ras, paste("data/derived_data/SocGrid/sampling1h/2022_03_30/",i2, ".grd", sep=""), format="raster", overwrite=T)
      
  #  }
 # }
  writeRaster(soc_grid, paste("data/derived_data/SocGrid/sampling1h/2022_03_30/",i2,"_grid", ".grd", sep=""), format="raster", overwrite=T)
}

names(soc_ras)
grid <- stack(paste("data/derived_data/SocGrid/sampling1h/",i2,"_","warm","_","night",".grd", sep=""))
grid <- stack(paste("data/derived_data/SocGrid/sampling1h/",i2,"_grid",".tif", sep=""))
plot(grid)


par(mfrow=c(1,1))
herto <- stack("data/derived_data/SocGrid/sampling2h/full/Hertogenwald_11.tif")
# study areas ----
grids_list <- list.files("data/derived_data/SocGrid/sampling1h/2022_03_25", pattern="_grid.grd")
study_areas <- str_remove_all(grids_list, ".grd") %>% str_remove_all(., "_grid") %>% unique() 
#saI <- study_areas[!study_areas %in% c("DNP_34_1","fanel_33_1",  "LageK_26_1", "MT_35_3", "MT_35_4")]
#saII <- study_areas[study_areas %in% c("DNP_34_1","fanel_33_1",  "LageK_26_1", "MT_35_3", "MT_35_4")] %>% 
#  stri_replace_last(fixed = "_", ".") 
#study_areas <- c(saI, saII)
study_areas <- study_areas[!study_areas %in% "Marche_14"]

# 3 Collinearity ----
library(usdm)
library(stringi)
vif(soc_ras)
removeCollinearity(soc_ras, plot = T, multicollinearity.cutoff = 0.5)


  
vif_df <- NULL
for (i in study_areas) {
  grid <- stack(paste("data/derived_data/SocGrid/sampling2h/full/",i,".tif", sep=""))
  names(grid) <- c("n_ind_abs", "n_ind_simult_max", "n_ind_simult_mean","dur_mean", "freq_mean",  
                   "dist_mean", "dist_sd", "dist_median",  "grp_change","SI","PX","CS","CR", "HAI")
  #plot(grid)
  j <- cellStats(is.na(grid), sum)
  j <- j/ncell(grid)
  grid <- grid[[which(j<.99)]]
  #grid_scale <- scale(grid, center=TRUE, scale=TRUE)
  #plot(grid_scale)
  vif_i <- vif(grid)
  vif_i$short_name <- i
  vif_df <- bind_rows(vif_df, vif_i)
}

vif_df %>%  
  filter_all(all_vars(!is.infinite(.))) %>%  group_by(Variables) %>% 
  dplyr::summarise(vif_mean=mean(VIF, na.rm=T),
                   vif_sd=sd(VIF, na.rm=T))

vif2 <- vifstep(raster_scale, th=10) 
sel3 <- dropLayer(grid_all, c("dist_sd","dist_median","Degree","freq_mean","n_ind_abs","n_ind_max_abs"))

# or....
library(virtualspecies)
removeCollinearity(soc_ras, plot = T, multicollinearity.cutoff = 0.5)
png("D:/PROJECTS/04_SocialScape/figures/collinearity.png",
    width=800, height=500, bg="transparent")
removeCollinearity(grid_all, plot = T, multicollinearity.cutoff = 0.5)
dev.off()

# 4 create Social_direct/indirect stack ----
library(stringi)
#study_loop <- read.csv("tables/study_loop.csv")
#study_loop <- gps_coord %>% group_by(study_name, short_name2) %>% tally() %>% dplyr::select(-n) %>% 
#  mutate(short_name2=stri_replace_last(short_name2,fixed = ".", "_") )


#write.csv(study_loop, "tables/study_loop.csv")

mean_narm = function(x,...){mean(x,na.rm=TRUE)}

remove <- c( "_","\\.", "1", "2","3", "4","5", "6", "7", "8", "9", "0")
study_list <- unique(gps_data$short_name2)
study_loop <- study_list %>% str_remove_all(paste(remove, collapse = "|")) %>% unique()
list_grids <- data.frame(grid_path=list.files("D:/PROJECTS/04_SocialScape/data/derived_data/SocGrid/sampling1h/2022_03_25/",full.names=T ))
list_cov <- data.frame(cov_path=list.files("D:/PROJECTS/xx_GIS/data/derived_data/stacks_2022", pattern = ".grd",full.names=T ))

dat_i <- NULL
dat_d <- NULL
dat_w <- NULL

dati_df <- NULL
datd_df <- NULL
datw_df <- NULL


for (i in 1:length(study_loop)) {
  #i <- 23
  print(study_loop[i])
  #i <- 24
  if (study_loop[i]=="FCNPwest") {study_loop[i]="FCNP_west"}
  if (study_loop[i]=="FCNPeast") {study_loop[i]="FCNP_east"}
  if (study_loop[i]=="OasiArezzo") {study_loop[i]="Oasi_Arezzo"}
  #i <- 21
  #for (j in c("night", "day")) {
    #j <- "night"
    #for (k in c("warm", "cold")) {
      #k <- "warm"
      #i <- "_Alb_39_cold_day_"
      grid <- NULL
      
      grid_i <- list_grids %>% filter(str_detect(grid_path, regex(study_loop[i], ignore_case = TRUE)))
      #grid_j <- grid_i %>% filter(str_detect(grid_path, j))
      #grid_k <- grid_j %>% filter(str_detect(grid_path, k))
      grid <- stack(grid_i$grid_path[1])

      #grid <- stack(paste("D:/PROJECTS/04_SocialScape/data/derived_data/SocGrid/sampling1h/",grid_i,"_",k,"_", j,".grd", sep=""))
      #plot(grid)
      #names(grid) <- c("n_ind_abs", "n_ind_simult_max", "n_ind_simult_mean","dur_mean", "freq_mean",  
      #                 "dist_mean", "dist_sd", "dist_median",  "grp_change","SI","PX","CS","CR", "HAI")
      sel <- dropLayer(grid, c("dist_sd","dist_median"))
      sel_scale <- scale(sel, center=TRUE, scale=TRUE)
      
      soc_direct <- overlay(sel_scale[["n_ind_simult_max"]],sel_scale[["n_ind_simult_mean"]],
                            sel_scale[["dur_mean"]],sel_scale[["dist_mean"]], fun=mean_narm)
      names(soc_direct) <- "Direct_interactions"
      soc_indirect <- overlay(sel_scale[["n_ind_abs"]],sel_scale[["freq_mean"]],
                              sel_scale[["grp_change"]], fun=mean_narm)
      names(soc_indirect) <- "Indirect_interactions"
      soc_withbetw <- overlay(sel[["SI"]],sel[["PX"]],
                              sel[["CS"]],sel[["CR"]],sel[["HAI"]],
                              #sel_scale[["grp_change"]],
                              fun=mean_narm)
      plot(soc_direct)
      plot(soc_indirect)
      par(mfrow=c(1,2))
      plot(soc_direct)
      plot(soc_withbetw)
      #plot(soc_withbetw_sd)
      #plot(soc_withbetw_min)
      if (sum(!is.na(values(soc_withbetw))) < 100){next}
      # select only areas where cells contains withbetw values
      soc_direct <- mask(soc_direct, soc_withbetw)
      soc_dir_wb <- overlay(soc_direct, soc_withbetw, fun=function(x,y){(x*y)})
      #soc_indirect <- mask(soc_indirect, soc_withbetw)
      #soc_indir_wb <- overlay(soc_indirect, soc_withbetw, fun=function(x,y){(x*y)})
      #if (table(is.na(values(soc_dir_wb)))[1]/ table(is.na(values(soc_dir_wb)))[2] < .01){next}
      #plot(soc_dir_wb)
      #plot(soc_indir_wb)
      #plot(soc_withbetw)
      
      #l <- study_loop[study_loop$study==i,]$covariates
      #if (l=="Geneva Basin - CH/FR") {l <- "Geneva Basin"}
      #l <- "Swabian Alps (Alb)"
      
      if (study_loop[i]=="doupov") {study_loop[i]="doupov"}
      cov <- NULL
      cov <- list_cov %>% filter(str_detect(cov_path, regex(study_loop[i], ignore_case = TRUE)))
      if (nrow(cov) == 1) {
        names_list <- names(stack(cov$cov_path))
        cov <- stack(cov$cov_path)
        names(cov) <- names_list
        }
      if (nrow(cov) == 2) {
        name_diff <- setdiff(names(stack(cov$cov_path[2])),names(stack(cov$cov_path[1])))
        if (length(name_diff) == 0) {
          names_list <- names(stack(cov$cov_path[1]))
          cov <- merge(stack(cov$cov_path[1]),stack(cov$cov_path[2]))
          names(cov) <- names_list}
        if (length(name_diff) > 0) {cov_diff <- dropLayer(stack(cov$cov_path[2]), name_diff)
          names_list <- names(stack(cov$cov_path[1]))
          cov <- merge(stack(cov$cov_path[1]),cov_diff)
          names(cov) <- names_list
        }}
      
      if (nrow(cov) == 3) {
        #name_diff <- setdiff(names(stack(cov$cov_path[1])),names(stack(cov$cov_path[3])))
         {
          cov_diff1 <- dropLayer(stack(cov$cov_path[1]), "distance_urban")
          cov_diff2 <- dropLayer(stack(cov$cov_path[2]), "distance_urban")
          names_list <- names(cov_diff2)
          cov <- merge(stack(cov_diff1),stack(cov_diff2),stack(cov$cov_path[3]))
          names(cov) <- names_list
         } 
        if (study_loop[i]!="sumava") {
          names_list <- names(stack(cov$cov_path[1]))
          cov <- merge(stack(cov$cov_path[1]),stack(cov$cov_path[2]),stack(cov$cov_path[3]))
          names(cov) <- names_list
        }
        
        
      }
      
      if (nrow(cov) == 4) {
        names_list <- names(stack(cov$cov_path[1]))
        cov <- merge(stack(cov$cov_path[1]),stack(cov$cov_path[2]),stack(cov$cov_path[3]),stack(cov$cov_path[4]))
        names(cov) <- names_list
        }
      if (nrow(cov) == 5) {
        cov_diff4 <- dropLayer(stack(cov$cov_path[4]), "distance_roads")
        cov_diff3 <- dropLayer(stack(cov$cov_path[3]), "distance_roads")
        names_list <- names(stack(cov$cov_path[1]))
        cov <- merge(stack(cov$cov_path[1]),stack(cov$cov_path[2]),cov_diff4,cov_diff4,stack(cov$cov_path[5]))
        names(cov) <- names_list
        }
      #plot(cov)
      soc_direct <- projectRaster(soc_direct, 
                                 crs ='+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs' , 
                                 method = "bilinear")
      soc_indirect <- projectRaster(soc_indirect, crs ='+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs' , 
                                   method = "bilinear")
      soc_withbetw <- projectRaster(soc_withbetw, crs ='+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs' , 
                                    method = "bilinear")
      
      predictors <- resample(cov,soc_direct, method="ngb")
      #predictors <- mask(predictors, wb_direct)
      #plot(predictors[[1]])
     
      # extract
      #cov_alb <- projectRaster(cov_alb, crs = crs(soc_direct))
      p_d <- rasterToPoints(soc_direct)
      p1_d <- p_d[,c("x","y")]
      pcov_d <- raster::extract(predictors, p1_d,cellnumbers=T)
      dat_d <- data.frame(cbind(p_d, pcov_d))
      if(nrow(dat_d)==0) {dat_d[1,]<-NA }
      dat_d$study <- study_loop[i]
      #dat_d$season <- k
      #dat_d$tod <- j
      
      p_i <- rasterToPoints(soc_indirect)
      p1_i <- p_i[,c("x","y")]
      pcov_i <- raster::extract(predictors, p1_i,cellnumbers=T)
      dat_i <- data.frame(cbind(p_i, pcov_i))
      if(nrow(dat_i)==0) {dat_i[1,]<-NA }
      dat_i$study <- study_loop[i]
      #dat_i$season <- k
      #dat_i$tod <- j
      
      p_w <- rasterToPoints(soc_withbetw)
      p1_w <- p_w[,c("x","y")]
      pcov_w <- raster::extract(predictors, p1_w,cellnumbers=T)
      dat_w <- data.frame(cbind(p_w, pcov_w))
      if(nrow(dat_w)==0) {dat_w[1,]<-NA }
      dat_w$study <- study_loop[i]
      #dat_w$season <- k
      #dat_w$tod <- j
      
      dati_df <- bind_rows(dati_df, dat_i)
      datd_df <- bind_rows(datd_df, dat_d)
      datw_df <- bind_rows(datw_df, dat_w)
      
    }
    
 # }
  

#}
dati_df %>% filter(study == "HainichNP")
datd_df %>% group_by(study) %>% tally()
saveRDS(dati_df, "data/derived_data/SocTbl/sampling1h/soc_indirect_wb.rds")
saveRDS(datd_df, "data/derived_data/SocTbl/sampling1h/soc_direct_wb.rds")
saveRDS(datw_df, "data/derived_data/SocTbl/sampling1h/soc_withbetw.rds")

# now go to 03_models

# OLD CODE ----


# Next step, extract covaiates at the high values and make a data.frame for model section
# 1 > 0.85, 0 for the rest

plot(classified_r, # what to plot
     col=viridis(3), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=.5, # vertical size of legend
     legend.width=2, # horizontal size of legend
     box=FALSE, axes=F
)
library(rasterVis)
library(viridis)
png("D:/PROJECTS/04_SocialScape/figures/social_map.png",
    width=800, height=500, bg="transparent")
par(mfrow=c(1,3))
plot(soc_indirect, # what to plot
     col=viridis(100), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=.5, # vertical size of legend
     legend.width=2, # horizontal size of legend
     box=FALSE, axes=F,
     main="indirect interactions"
)
plot(soc_direct, # what to plot
     col=viridis(100), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=.5, # vertical size of legend
     legend.width=2, # horizontal size of legend
     box=FALSE, axes=F,
     main="direct interactions"
)
plot(soc_withbetw, # what to plot
     col=viridis(100), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=.5, # vertical size of legend
     legend.width=2, # horizontal size of legend
     box=FALSE, axes=F,
     main="within-between groups"
)
dev.off()

test <- stack(soc_indirect,soc_direct)
plot(test)
 
class_uncertainty %>% group_by(classification) %>% dplyr::summarise(mean_unc=mean(uncertainty))
write.csv(ls[["parameters"]][["mean"]], "D:/PROJECTS/04_SocialScape/tables/classification2.csv")



ModPred <- predict.Mclust(ls, table_test[,c(1:7)]) # prediction
Pred_ras <- grid[[1]] # establishing a rediction raster
values(Pred_ras) <- NA # set everything to NA
# set values of prediction raster to corresponding classification according to rowname
values(Pred_ras)[as.numeric(rownames(table_test))] <- as.vector(ModPred$classification)
# uncertainty raster
unc_ras <- grid[[1]] # establishing a rediction raster
values(unc_ras) <- NA # set everything to NA
# set values of prediction raster to corresponding classification according to rowname
values(unc_ras)[as.numeric(rownames(class_uncertainty))] <- as.vector(class_uncertainty$uncertainty)
names(unc_ras) <- "uncertainty"
names(Pred_ras) <- "class"
plot(grid[[1]])
plot(Pred_ras)
library(rasterVis)
library(viridis)
library(httr)
png("D:/PROJECTS/04_SocialScape/figures/cluster_prediction.png",
    width=800, height=500, bg="transparent")
plot(Pred_Unc, # what to plot
     col=viridis(6), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=.5, # vertical size of legend
     legend.width=2, # horizontal size of legend
     box=FALSE, axes=F
)
dev.off()
table(values(Pred_ras))
Pred_Unc <- stack(Pred_ras, unc_ras)
# Reproject to fit covirables projection (epsg 3035)
crs(Pred_Unc) <- "+init=epsg:32632"
clusters <- projectRaster(Pred_Unc, crs = '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs',
                          method="ngb")
writeRaster(clusters, "D:/PROJECTS/04_SocialScape/data/derived_data/clusters.tif",format="GTiff", overwrite=T)

# Cluster interpretation
clust <- data.frame(ls[["parameters"]][["mean"]]) %>% tibble::rownames_to_column( var = "metric") 
# quantile estimation
clust_quantile <- na.omit(as.data.frame(sel3)) %>% melt() %>% group_by(metric=variable) %>% summarise(q = list(quantile(value))) %>% 
  unnest_wider(q) %>% collect() %>% data.frame()
# clust_quantile <- clust %>% pivot_longer(cols=c(-metric))  %>% 
#   group_by(metric) %>% 
#   summarise(q = list(quantile(value))) %>% 
#   unnest_wider(q) %>% collect() %>% data.frame()
clust_rank <- clust %>% pivot_longer(cols=c(-metric)) %>%  group_by(metric) 
clustq <- left_join(clust_rank, clust_quantile, by="metric")
names(clustq) <- c("metric", "name",   "value",  "q0","q25","q50","q75", "q100")



clust_class <-clustq %>% group_by(metric) %>% 
  mutate(cluster=parse_number(name),
         nrank=case_when(value <=q25~1,
                         value > q25 & value <=q75~2,
                         value > q75~3))

type_use <- clust_class %>% filter(grepl("Weight",metric)) %>% 
  mutate(property=case_when(metric=="Weight" & nrank==1~1,#"low",
                            metric=="Weight" & nrank==2~2,#"medium",
                            metric=="Weight" & nrank==3~3)) %>% #"high")) %>% ungroup() %>% 
  dplyr::select(cluster, property) %>% mutate(type="use_intensity")

type_corridor <- clust_class %>% filter(metric=="Between"|metric=="Directionality"|metric=="Speed") %>%
  dplyr::select(cluster, nrank) %>% pivot_wider(names_from = metric, values_from=nrank) %>% 
  mutate(property=case_when(Between==2 & Speed==2 & Directionality==2~2, #"medium_slow",
                            Between==2 & Speed==3 & Directionality==2~2, #"medium_fast",
                            Between==3 & Speed==1 & Directionality==2~4, #"high_slow",
                            Between==3 & Speed==2 & Directionality==2~5,#"high_medium",
                            Between==3 & Speed==3 & Directionality==2~6)) %>% #"high_fast")) %>% 
  dplyr::select(cluster, property) %>% mutate(type="corridor")

type_social_direct <- clust_class %>% filter(metric=="dur_mean"|metric=="dist_mean"|metric=="n_ind_mean_abs") %>% 
  dplyr::select(cluster, nrank) %>% pivot_wider(names_from = metric, values_from=nrank) %>%
  mutate(property=case_when(n_ind_mean_abs==2 & dur_mean==1 & dist_mean==2~1,#"low",
                            n_ind_mean_abs==2 & dur_mean==2 & dist_mean==2~2,#"medium",
                            n_ind_mean_abs==1 & dur_mean==3 & dist_mean==2~2,#"medium",
                            n_ind_mean_abs==3 & dur_mean==2 & dist_mean==2~3, #"high",
                            n_ind_mean_abs==3 & dur_mean==3 & dist_mean==2~3,#"high",
                            n_ind_mean_abs==2 & dur_mean==3 & dist_mean==2~3)) %>% #"high")) %>% 
  dplyr::select(cluster, property) %>% mutate(type="social_direct")

type_social_indirect <- clust_class %>% filter(metric=="grp_change") %>% 
  mutate(property=case_when(nrank==1~1, #"low",
                            nrank==2~2,#"medium",
                            nrank==3~3)) %>% #"high")) %>% ungroup() %>% 
  dplyr::select(cluster, property) %>% mutate(type="social_indirect")

type_all <- bind_rows(type_use, type_corridor, type_social_direct, type_social_indirect)                            

# plot the different type

new <- list()
for (i in c("use_intensity", "corridor", "social_direct", "social_indirect")) {
  #i <- "corridor"
  xx_clus <- type_all %>% filter(type==i) #%>% mutate(srank=sub("^[^_]*_", "", property))
  rcl <- xx_clus %>% ungroup() %>%  dplyr::select(cluster, property) %>% data.frame() %>% as.matrix()
  #if (i=="social_direct") {rcl <- xx_clus %>% group_by(cluster) %>% dplyr::summarise(avg=round(mean(nrank),0)) %>% as.matrix()}
  xx_ras <- reclassify(Pred_ras, rcl)
  names(xx_ras) <- i
  new[[i]] <- xx_ras
}
stack_ras <- stack(new)
plot(stack_ras)
arg <- list(at=c(1,2,3,4,5,6), labels=c("Low","Medium","High","High_Slow","High_Medium","High_Fast"))
png("D:/PROJECTS/04_SocialScape/figures/cluster_property.png",
    width=800, height=500, bg="transparent")
par(mfrow=c(2,2))
plot(stack_ras[[1]], col=viridis(3), colNA = "white", legend.shrink=.5, legend.width=1,axis.args=list(at=c(1,2,3), labels=c("Low","Medium","High")),
     box=F, axes=F, main="Intensity of use")
plot(stack_ras[[2]], col=viridis(4), colNA = "white", legend.shrink=.5, legend.width=1,axis.args=list(at=c(2,4,5,6), labels=c("Medium","High_Slow","High_Medium","High_Fast")),
     box=F, axes=F, main="Corridor")
plot(stack_ras[[3]], col=viridis(3), colNA = "white", legend.shrink=.5, legend.width=1,axis.args=list(at=c(1,2,3), labels=c("Low","Medium","High")),
     box=F, axes=F, main="Direct social interactions")
plot(stack_ras[[4]], col=viridis(3), colNA = "white", legend.shrink=.5, legend.width=1,axis.args=list(at=c(1,2,3), labels=c("Low","Medium","High")),
     box=F, axes=F, main="Indirect social interactions")
dev.off()
par(mfrow=c(1,1))

# Reproject to fit covirables projection (epsg 3035)
crs(stack_ras) <- "+init=epsg:32632"
stack_ras <- projectRaster(stack_ras, crs = '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs',
                           method="ngb")

# save

writeRaster(stack_ras[[1]], "D:/PROJECTS/04_SocialScape/data/derived_data/raster_use_intensity.tif",format="GTiff", overwrite=T)
writeRaster(stack_ras[[2]], "D:/PROJECTS/04_SocialScape/data/derived_data/raster_corridor.tif",format="GTiff", overwrite=T)
writeRaster(stack_ras[[3]], "D:/PROJECTS/04_SocialScape/data/derived_data/raster_social_direct.tif",format="GTiff", overwrite=T)
writeRaster(stack_ras[[4]], "D:/PROJECTS/04_SocialScape/data/derived_data/raster_social_indirect.tif",format="GTiff", overwrite=T)

#writeRaster(stack(soc_direct,soc_indirect,soc_withbetw), "D:/PROJECTS/04_SocialScape/data/derived_data/Alb_Interactions_maps.grd")

# if working by 'class'

#pp_indirect <- quantile(soc_indirect, c(0, 0.15, 0.85))
#ix_indirect <- findInterval(getValues(soc_indirect), pp_indirect)
#class_indirect <- setValues(soc_indirect, ix_indirect)
#class_indirect[class_indirect<3] <- 0
#class_indirect[class_indirect==3] <- 1

#pp_direct <- quantile(soc_direct, c(0, 0.15, 0.85))
#ix_direct <- findInterval(getValues(soc_direct), pp_direct)
#class_direct <- setValues(soc_indirect, ix_direct)
#class_direct[class_direct<3] <- 0
#class_direct[class_direct==3] <- 1
