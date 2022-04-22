library(here)
library(moveNT)
source('D:/PROJECTS/04_SocialScape/R/functions_MoveNT.R')
source('D:/PROJECTS/04_SocialScape/R/functions.R')
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
devtools::install_github("mkearney/uslides")
library(uslides)
# 1 data import ----
load("D:/PROJECTS/02_Contacts/data/raw_data/gps_animals.Rda")
load("D:/PROJECTS/02_Contacts/data/raw_data/gps_animals2.Rda")


gps_animals <- gps_animals %>% 
  distinct(animals_original_id, acquisition_time, .keep_all = T) 
names(gps_animals)
gps_animals %>% group_by(short_name) %>% tally() %>% View()
# 2 GPS data preparation ----

# merge with study_areas_mcp to properly check number of animals per 'mcp'
sas <- st_read("D:/PROJECTS/xx_GIS/data/study_areas_mcp.shp")
st_crs(sas) <- 4326
names(sas) <- c("study_areas_id","research_groups_id","short_name","study_name","geometry")
sas_poly <- st_cast(sas, "POLYGON")
sas_poly$id <- 1:nrow(sas_poly)
sas_poly <- sas_poly %>% rownames_to_column() %>% 
  mutate(short_name2=paste(short_name,rowname, sep="_")) %>% dplyr::select(-rowname)
# make gps sf object
gps_sf <- st_as_sf(gps_animals,coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
gps_sas <- st_join(gps_sf, sas_poly, join = st_within)
coords <- st_coordinates(gps_sas)
gps_coord <- cbind(gps_sas, coords) %>% collect() %>% data.frame() %>% 
  group_by(short_name2) %>% mutate(n_ind=n_distinct(animals_id)) %>%
  filter(n_ind > 5) %>% ungroup()
# count number of individual by study area
gps_utm %>% group_by(short_name2) %>% mutate(n_ind=n_distinct(animals_id)) %>% 
  dplyr::select(short_name2, n_ind) %>% distinct(short_name2, n_ind) %>% arrange(n_ind) %>% View()
gps_coord %>% dplyr::select(short_name2, n_ind) %>% distinct(short_name2, n_ind) %>% arrange(n_ind) %>% View()


# make a track
gps_utm <- mk_track(gps_coord, .x = X, .y = Y, .t =acquisition_time,id=animals_original_id,
                    crs = sp::CRS("+proj=longlat +ellps=WGS84"), all_cols = T) %>%
  transform_coords(sp::CRS("+init=epsg:32632")) %>%
  mutate(x=as.numeric(x_), y=as.numeric(y_), acquisition_time=as.POSIXct(t_),
         animals_original_id=case_when(is.na(animals_original_id)~animals_original_name,
                                       TRUE ~ as.character(animals_original_id))) %>% 
  dplyr::select(short_name2, animals_id, animals_original_id,gps_sensors_id, acquisition_time,t_,
                x_, y_, start_time, end_time, sex, age=age_class_code_capture) %>% 
  arrange(animals_id, acquisition_time) %>% 
  time_of_day()


# Resample
# to have all data at 1 loc /2 hour (to keep some key study areas)
gps2h <- gps_utm %>% #filter(short_name2 %in% c("StHubert_17","Hertogenwald_11","Famenne_8","CondrozLg_3",
                           #                          "Couvin_6.1","Namur_16","Marche_14","Froidchapelle_10",
                           #                          "Couvin_6", "Thuin_18")) %>% 
  nest(data=!"animals_original_id") %>% 
  mutate(steps = purrr::map(data, function(x)
    x %>% amt::track_resample(rate = hours(2), tolerance = minutes(15)))) %>% 
  dplyr::select(animals_original_id, steps) %>%  
  unnest(cols = steps) %>% 
  time_of_day()

gps1h <- gps_utm %>% #filter(short_name2 %in% c("StHubert_17","Hertogenwald_11","Famenne_8","CondrozLg_3",
  #                          "Couvin_6.1","Namur_16","Marche_14","Froidchapelle_10",
  #                          "Couvin_6", "Thuin_18")) %>% 
  nest(data=!"animals_original_id") %>% 
  mutate(steps = purrr::map(data, function(x)
    x %>% amt::track_resample(rate = hours(1), tolerance = minutes(10)))) %>% 
  dplyr::select(animals_original_id, steps) %>%  
  unnest(cols = steps)  %>% 
  time_of_day()

# Add summer-spring and autum-winter period
gps1h <- gps1h %>% mutate(season=case_when(month(acquisition_time) %in% c(1:3, 10:12)~"cold",
                                                       month(acquisition_time) %in% c(4:9)~"warm"),
                          x=x_, y=y_)  
  
gps2h <- gps2h %>% mutate(season=case_when(month(acquisition_time) %in% c(1:3, 10:12)~"cold",
                                           month(acquisition_time) %in% c(4:9)~"warm"),
                          x=x_, y=y_)
gps1h %>% group_by(short_name2,season) %>% dplyr::summarise(start=min(acquisition_time),
                                                 end=max(acquisition_time)) %>% 
  View()

gps2h %>% group_by(short_name2,season) %>% dplyr::summarise(start=min(acquisition_time),
                                                            end=max(acquisition_time)) %>% 
  View()

gps_cold <- gps_belgium %>% filter(season=="cold")
gps_warm <- gps_belgium %>% filter(season=="warm")

saveRDS(gps_utm,"data/gps_utm.rds")
saveRDS(gps2h,"data/gps2h.rds")
saveRDS(gps1h,"data/gps1h.rds")



# 3 Social cohesion ----
# using WILDLIFE DI
# Convert gps data into ltraj object (adehabitat)
gps_utm <- readRDS("data/gps1h.rds") %>% st_as_sf(coords = c("x_", "y_"), crs = 32632, agr = "constant")
library(adehabitatHR)
library(wildlifeDI)
library(adehabitatLT)
library(rgeos)
library(tidyverse)
library(lubridate)

social_index <- NULL
con_sf_all <- NULL

setdiff(unique(gps_utm$short_name2),unique(social_indexI$study_area))
for (i in unique(gps_utm$short_name2)) { #c("Boxholm_49","Hertogenwald_11","Bialowieza_20")) { #setdiff(unique(gps_utm$short_name2),unique(social_indexI$study_area))
  #i <- "LageK_26.1"
  print(i)
  gps_i <-  gps_utm %>% filter(short_name2==i) %>% mutate(JDate=yday(acquisition_time)) %>%
    dplyr::select(animals_id, gps_sensors_id,animals_original_id, acquisition_time,x, y, sex, age, JDate) %>% 
    distinct(animals_original_id,acquisition_time, .keep_all= TRUE)
  #coordinates(gps_i) <- ~x+y
  raw_gps <- as.ltraj(st_coordinates(gps_i),
                      date=gps_i$acquisition_time,
                      id=gps_i$animals_original_id, typeII=TRUE,
                      proj4string=CRS(SRS_string = "EPSG:32632"))
  
  # regularization
  refda <- min(gps_i$acquisition_time)
  gps_NA <- setNA(raw_gps,refda,60,units="min")
  
  boar.traj <- sett0(gps_NA,refda,60,units="min")
  boar.traj <- na.omit(boar.traj)
  ######## test ##############
  # 
  # 
  # #boar2 <- boar.traj[burst(boar.traj) %in% c(pair$id1,pair$id2)]
  # gps_pair <- gps_i %>% filter(animals_original_id %in% c("HAN_038_039","HAN_027_028"))
  # pair_traj <- as.ltraj(st_coordinates(gps_pair),
  #                        date=gps_pair$acquisition_time,
  #                        id=gps_pair$animals_original_id, typeII=TRUE,
  #                        proj4string=CRS(SRS_string = "EPSG:32632"))
  # 
  # plt <- dcPlot(pair_traj,tc=15*60,dmax=1000)
  # doecons <- conProcess(pair_traj,dc=50,tc=10*60)
  # doephas <- conPhase(doecons, pc=10*60)
  # conSummary(doephas)
  # doepair <- conPairs(doephas)
  # doetemp <- conTemporal(doephas,units='mins')
  # 
  # doepair$hod <- as.POSIXlt(doepair$date)$hour + as.POSIXlt(doepair$date)$min / 60  #convert POSIX to hours
  # doetemp$hod <- as.POSIXlt(doetemp$start_time)$hour + as.POSIXlt(doetemp$start_time)$min / 60  #convert POSIX to hours
  # doepair$dom <- as.POSIXlt(doepair$date)$mday
  # hist(doepair$dom,breaks=0:31)
  # hist(doepair$hod,breaks=0:24) #Figure 2b
  # hist(doetemp$hod,breaks=0:24) #Figure 2c
  # hist(as.numeric(doetemp$duration), breaks = 10) #figure 2d
  # 
  # con_sf <- wildlifeDI::conSpatial(doephas,type='point')             # Get points of all contacts
  # 
  # sf_pt <- wildlifeDI::ltraj2sf(boar.traj.pair)  # Turn all fixes into sf points
  # plot(st_geometry(sf_pt),col='grey',pch=20)
  # plot(st_geometry(con_sf),col='black',pch=20,add=T)
  
  
  
  #######end test##############
  
  # Create fuzzy table to assess association value between all potential pairs of wild boar
  id1 <- unique(gps_i$animals_original_id)
  id2 <- id1
  data <- list(
    id1 = id1,
    id2 = id2
  )
  data %>%
    cross_df() -> data_df
  
  #Kenward's Coefficient of Sociality
  to_fill <- data.frame(ind1=NA, ind2=NA,si=NA,prox=NA, cs=NA, ca=NA, hai=NA, cr=NA,
                        nContacts = NA, propContacts = NA, longPhase_hr = NA, meanPhase_hr = NA)
  df_all <- NULL
  con_sf_i <- NULL
  dc <- 50
  tc <- 10*60
  for (j in 1:nrow(data_df) ) {
    #j <- 2
    print(j)
    pair <- data_df[j,]
    
    if (pair$id1==pair$id2) {next}
    pair1 <- boar.traj[adehabitatLT::id(boar.traj)==pair$id1]
    pair2 <- boar.traj[adehabitatLT::id(boar.traj)==pair$id2]
    #i
    check <- checkTO(pair1,pair2)
    if (check$TO==FALSE) {next}
    # Static interaction analysis (SI)
    pts_1 <- SpatialPoints(ld(pair1)[,1:2])
    pts_2 <- SpatialPoints(ld(pair2)[,1:2])
    #kde_1 <- kernelUD(pts_1)
    #kde_2 <- kernelUD(pts_2)
    #hr_1 <- getverticeshr(kde_1,95)
    #hr_2 <- getverticeshr(kde_2,95)
    
    hr_1 <- mcp(pts_1,percent=100)
    hr_2 <- mcp(pts_2,percent=100)
    #plot(hr_2); plot(hr_1, add=T)
    if (is.null(gIntersection(hr_1,hr_2))) {to_fill$ind1 <- pair$id1
    to_fill$ind2 <- pair$id2
    to_fill$si <- NA
    to_fill$prox <- NA
    to_fill$cs <- NA
    to_fill$hai <- NA
    to_fill$cr <- NA
    to_fill$nContacts <- NA
    to_fill$propContacts <- NA
    to_fill$longPhase_hr <- NA
    to_fill$meanPhase_hr <- NA
    #to_fill$iab_attract <- NA
    #to_fill$iab_avoid <- NA
    to_fill$study_area <- i
    df_all <- rbind(df_all, to_fill)
    next} else {SI <- gArea(gIntersection(hr_1,hr_2))/gArea(gUnion(hr_1,hr_2))}
    # Simultaneous Fixes in time! (SF)
    #SF <- try(GetSimultaneous(pair1,pair2,tc=tc), silent = TRUE)
    #if (class(SF)=="try-error") {SF <- NA}
    #SF_1 <- nrow(SF[[1]]) / nrow(pair1[[1]])
    #SF_2 <- nrow(SF[[2]]) / nrow(pair2[[1]])
    to_fill$si <- SI
    to_fill$ind1 <- pair$id1
    to_fill$ind2 <- pair$id2
    # Dynamic Interaction
    #Proximity Analysis (PX)
    
    PX <- try(Prox(pair1,pair2, tc=tc, dc=dc), silent = TRUE)
    if (class(PX)=="try-error") {PX <- NA}
    to_fill$prox <- PX
    # Ca - Coefficient of association
    CA <- try(Ca(pair1,pair2, tc=tc, dc=dc), silent = TRUE)
    if (class(CA)=="try-error") {CA <- NA}
    to_fill$ca <- CA
    # Cs - Coefficient of sociality
    CS <- try(Cs(pair1,  pair2, tc=tc), silent = TRUE)
    if (class(CS)=="try-error") {CS <- NA}
    if (!is.na(CS)) {to_fill$cs <- CS$Cs} else {to_fill$cs <- NA}
    # HAI - Half-weight association index
    oz <- gIntersection(hr_1, hr_2)
    HAI <- try(HAI(pair1, pair2, oz, tc=tc, dc=dc), silent = TRUE)
    if (class(HAI)=="try-error") {HAI <- NA}
    to_fill$hai <- HAI
    # IAB - Interaction Statistic
    #IAB <- try(IAB(pair1, pair2, tc=tc, dc=dc), silent = TRUE)
    #if (class(IAB)=="try-error") {IAB <- NA}
    # Cr - Correlation coefficient
    CR <- try(Cr(pair1, pair2, tc=tc), silent = TRUE)
    if (class(CR)=="try-error") {CR <- NA}
    if (!is.na(CR)) {to_fill$cr <- CR} else {to_fill$cr <- NA}
    # DI - Dynamic interaction index
    #DI <- try(DI(pair1, pair2, tc=tc))
    #if (class(DI)=="try-error") {DI <- NA}
    
    # Contact data
    # from https://cran.r-project.org/web/packages/wildlifeDI/vignettes/wildlifeDI-vignette-contact_analysis.html#processing-contacts
    gps_pair <- gps_i %>% filter(animals_original_id %in% c(pair$id1,pair$id2))
    pair_traj <- as.ltraj(st_coordinates(gps_pair),
                          date=gps_pair$acquisition_time,
                          id=gps_pair$animals_original_id, typeII=TRUE,
                          proj4string=CRS(SRS_string = "EPSG:32632"))
    doecons <- conProcess(pair_traj,dc=dc,tc=tc)
    #if (class(CR)=="doecons") {CR <- NA}
    doephas <- conPhase(doecons, pc=15*60)
    con_df <- try(conSummary(doephas), silent = TRUE)
    if (class(con_df )=="try-error") {to_fill$nContacts <- NA
    to_fill$propContacts <- NA
    to_fill$longPhase_hr <- NA
    to_fill$meanPhase_hr <- NA} else {
      con_sf <- wildlifeDI::conSpatial(doephas,type='point')      
      con_sf$study <- i
      # fill in the table
      to_fill$nContacts <- con_df$result[2]
      to_fill$propContacts <- con_df$result[2] / con_df$result[1]
      to_fill$longPhase_hr <- con_df$result[4] / 3600
      to_fill$meanPhase_hr <- con_df$result[5] / 3600
    }
    
    to_fill$study_area <- i
    df_all <- rbind(df_all, to_fill)
    con_sf_i <- bind_rows(con_sf_i, con_sf)
  }
  social_index <- bind_rows(social_index, df_all)
  con_sf_all <- bind_rows(con_sf_all, con_sf_i)
}
library(mapview)
mapview(con_sf_i)
View(df)
social_df <- social_index %>% distinct(ca, cr, cs, hai, prox, si,nContacts, study_area, .keep_all = T)
saveRDS(social_index, "tables/social_index_2022.rds" )
st_write(con_sf_all, "data/derived_data/contacts_sf.csv", layer_options = "GEOMETRY=AS_XY")

social_index %>% distinct(ca, cr, cs, hai, prox, si,nContacts, study_area, .keep_all = T) %>% group_by(study_area) %>% dplyr::summarise(maxSI = max(si,na.rm=T),
                                              meanSI = mean(si,na.rm=T),
                                                     minSI = min(si,na.rm=T))
social %>% filter(study_area =="HainichNP_15")
social_index %>% distinct(ca, cr, cs, hai, prox, si,nContacts, study_area, .keep_all = T) %>%
  pivot_longer(cols = c(si,prox,cs,ca,hai,cr), names_to="metric") %>% 
  #filter(value > 0) %>% 
  ggplot(aes(x = value)) + geom_histogram(bins = 20) +
  facet_wrap(~metric, scales = "free")
View(social_index)
# kmeans clustering ----
# from https://uc-r.github.io/kmeans_clustering
df <- USArrests
df <- social_index %>% distinct(ca, cr, cs, hai, prox, si,nContacts, study_area, .keep_all = T) %>%
  dplyr::select(ca, cr, cs, hai, prox, si) %>% na.omit() %>% scale()

k2 <- kmeans(df, centers = 4, nstart = 25)
str(k2)
fviz_cluster(k2, data = df)
# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(df, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

fviz_nbclust(df, kmeans, method = "wss")
fviz_nbclust(df, kmeans, method = "silhouette")
library(cluster)
gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
final <- kmeans(df, 3, nstart = 25)
fviz_cluster(final, data = df)
social_df %>% dplyr::select(study_area, ind1, ind2, ca, cr, cs, hai, prox, si) %>% na.omit() %>% 
  mutate(Cluster = final$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

cluster <- social_df %>% dplyr::select(study_area, ind1, ind2, ca, cr, cs, hai, prox, si) %>% na.omit() %>% 
  mutate(Cluster = final$cluster,
         cluster_class = case_when(Cluster == 1 ~ "between-group",
                                   Cluster == 2 ~ "intermediate",
                                   Cluster == 3 ~ "within-group")) %>% 
  dplyr::select(ind1, ind2, cluster_class)

social_final <- left_join(social_df, cluster, by= c("ind1", "ind2")) %>% 
  filter(!is.na(cluster_class))
View(social_final)
nrow(social_final)
social_final %>% 
  pivot_longer(cols = c(si,prox,cs,ca,hai,cr), names_to="metric") %>% 
  filter(value > 0) %>% 
  ggplot(aes(x = value, fill = cluster_class)) + geom_density( alpha = .5) +
  facet_wrap(~metric, scales = "free")

# associate the cluster to the sf object


# 4 Movescape ----
gps_traj <- readRDS(file = "D:/PROJECTS/02_Contacts/materials/Bastille_Rousseau/gps_ltraj.rds")
grid<-loop(gps_traj, 200)


# interpolation
grid2 <-interpolation(gps_traj, grid)
# layers = weight, degree, between, speed, dotp
plot(grid2[[3]])

# mosaic individual
mean_weight<-mosaic_network(grid, index=2, sc=T, fun=mean)
names(mean_weight) <- "Weight"
mean_degree <- mosaic_network(grid, index=4, sc=T, fun=mean)
names(mean_degree) <- "Degree"
max_between <- mosaic_network(grid, index=5, sc=T, fun=max)
names(max_between) <- "Between"
mean_speed <- mosaic_network(grid, index=11, sc=T, fun=mean)
names(mean_speed) <- "Speed"
mean_ta <- mosaic_network(grid, index=13, sc=T, fun=mean)
names(mean_ta) <- "Directionality"
par(mfrow=c(1,3))
plot(mean_weight, main="Weight");plot(mean_degree, main="Degree");plot(max_between, main="Betweeness")
# mosaic on interpolated value
mean_weight2 <-mosaic_network(grid2, index=1, sc=T, fun=mean)
names(mean_weight2) <- "Weight"
mean_degree2 <- mosaic_network(grid2, index=2, sc=T, fun=mean)
names(mean_degree2) <- "Degree"
max_between2 <- mosaic_network(grid2, index=3, sc=T, fun=max)
names(max_between2) <- "Between"
mean_speed2 <- mosaic_network(grid2, index=4, sc=T, fun=mean)
names(mean_speed2) <- "Speed"
mean_ta2 <- mosaic_network(grid2, index=5, sc=T, fun=mean)
names(mean_ta2) <- "Directionality"

par(mfrow=c(2,2))
plot(mean_weight2, main="Weight");plot(mean_degree2, main="Degree");
plot(max_between2, main="Betweeness");plot(mean_speed2, main="Speed")
par(mfrow=c(1,1))
plot(mean_ta, main="Directionality")

table_grid<-table_cluster(gps_traj, grid)
#Showing the first few rows of the table created.
head(table_grid)
#B- Individual-level clustering
#The first step of the analysis is to apply the clustering to each individual. *ind_clust* apply a mixture model to each individual. It is possible to specify the maximum number of clusters (here 8) and also the covariates to use for the clustering, but the function automatically selects the optimal number of clusters (based on BIC). In our case, 2 individuals had 6 clusters, 2 had seven, and 2 had eight clusters. *ls_ind* simply return a list object with each element representing a single individual.
ls_ind<-ind_clust(table_grid, max.n.clust=8)
#Showing the number of individuals with 6, 7, and 8 clusters
#   (i.e. no individual had less than 6 clusters)
table(unlist(lapply(ls_ind, function(x) x$G)))

#C- Population-level clustering
#After performing the individual clustering, a second clustering is applied via *pop_clust*. This second clustering takes the ouptut of *ind_clust* and will identify which individual clusters could be considered as one population-level clusters. The function automatically selects the optimal number of clusters (based on BIC). It is possible for two clusters from the same individual to be in the same population-level cluster. Likewise, it is possible that a population level cluster does not have all individuals. Here, 3 different population clusters were calculated. The second line extract the center (mean) of each cluster which is helpful in interpreting their meaning. The first cluster was heavily used (weight), well connected (degree), and important for connectivity (betweenness), but albatross were moving slowly and not linearly in them. The second cluster was a cluster with intermediate use, not important for connectivity and still with meandering movement. The third cluster was important for connectivty and albatross were  moving fast and linearly in it. We also extract the proportion of each cluster.
pop<-pop_clust(gps_traj, ls_ind)
pop
#Paremeters associated to each cluster (i.e. center)
pop[[1]]$parameters$mean
#Proportion of each cluster (how frequent it is spatially)
pop[[1]]$parameters$pro

#D- Mapping and results export
#After performing the population level cluster, the function  *clust_stack* recombines the individual and population level clustering and produce a *stack* object for each individual albatross showing the most likely cluster, and also the probability of observing each cluster (uncertainty) in any given pixel. This individual level data (but which contains the population level clustering) can be used in a regression based analysis as presented in the manuscript or simply mapped. We developed two functions to produce these maps. *pop_stack* generate for each population cluster, rasters showing if at least one individual is using this pixel for this specific cluster. *pop_overl* display for each pixel all potential use observed, for example a pixel having the value *123* will have at least one individual using this pixel as cluster 1, another individual using it as 2, and another individual as 3. We show how frequent each combination are using the *table* function. These object can be exported to be used in other software using the *writeRaster* function.

clust_stack<-clust_stack(grid, pop, ls_ind, table_grid)
#Plotting the stack for the first individual
plot(clust_stack[[1]])
pop_stack<-pop_stack(clust_stack)
#Plotting the stack object. Each raster shows if at least one
# individual is using the pixel as a specific cluster.
plot(pop_stack)
#writeRaster(pop_stack, "Network.tif", format="GTiff",  bylayer=T, suffix="names", overwrite=T)
pop_overl<-pop_overl(clust_stack)
#Table showing how frequent each pixels is. A value of 123 indicates
# that at least one individual is using the pixel as 1, one individual
# is using it as 2, and one as 3
table(values(pop_overl))
#Plotting of the overlap raster. The legend makes it hard to see,
# but values range from 0 (no animal present) to 123
# (each type of used is observed in the pixel).
plot(pop_overl)

## see others... Social grid ----
socgrid <- socialGrid(data=gps_traj, res=200, crs=32632, contact_threshold='5 minutes')
soc_tbl <- socgrid[[1]] %>% na.omit()

soc_ras <- socgrid[[2]]
crs(soc_ras) <- "+init=epsg:32632"
plot(soc_ras, box=FALSE, axes=F)
# join move and social scape variables
#grid_all2 <- stack(soc_ras, mean_weight2, mean_degree2, max_between2, mean_speed2, mean_ta2)
grid_all <- stack(soc_ras, mean_weight, mean_degree, max_between, mean_speed, mean_ta)
plot(grid_all)

# 3 Collinearity ----
library(virtualspecies)
removeCollinearity(soc_ras, plot = T, multicollinearity.cutoff = 0.5)
png("D:/PROJECTS/04_SocialScape/figures/collinearity.png",
    width=800, height=500, bg="transparent")
removeCollinearity(grid_all, plot = T, multicollinearity.cutoff = 0.5)
dev.off()

# or....
library(usdm)
vif(grid_all)
vif2 <- vifstep(raster_scale, th=10) 
sel3 <- dropLayer(grid_all, c("dist_sd","dist_median","Degree","freq_mean","n_ind_abs","n_ind_max_abs"))
png("D:/PROJECTS/04_SocialScape/figures/cluster_entry.png",
    width=800, height=500, bg="transparent")
plot(sel3)
dev.off()
vif(sel3)

# 4 Unsupervised classification ----
table_test <- as.data.frame(sel3)   %>% na.omit()
library(mclust)
dataBIC <- mclustBIC(table_test[,c(1:8)])
print(summary(dataBIC))
plot(dataBIC)
ls <- Mclust(table_test[,c(1:8)], G=6, modelNames="VVI")
ls_boot <- MclustBootstrap(ls, nboot = 100, type = "bs",verbose = TRUE)
summary(ls_boot, what = "se")
summary(ls_boot, what = "ci")

ls[["parameters"]][["mean"]]
ls[["parameters"]][["variance"]]
ls[["z"]]
class_uncertainty <- data.frame(uncertainty=round(ls[["uncertainty"]],2),
                          classification=ls[["classification"]]) %>% arrange(uncertainty) 
class_uncertainty %>% group_by(classification) %>% dplyr::summarise(mean_unc=mean(uncertainty))
write.csv(ls[["parameters"]][["mean"]], "D:/PROJECTS/04_SocialScape/tables/classification2.csv")

par(mfrow=c(1,1))
?predict.Mclust
ModPred <- predict.Mclust(ls, table_test[,c(1:8)]) # prediction
Pred_ras <- socgrid[[3]] # establishing a rediction raster
values(Pred_ras) <- NA # set everything to NA
# set values of prediction raster to corresponding classification according to rowname
values(Pred_ras)[as.numeric(rownames(table_test))] <- as.vector(ModPred$classification)
# uncertainty raster
unc_ras <- socgrid[[3]] # establishing a rediction raster
values(unc_ras) <- NA # set everything to NA
# set values of prediction raster to corresponding classification according to rowname
values(unc_ras)[as.numeric(rownames(class_uncertainty))] <- as.vector(class_uncertainty$uncertainty)
names(unc_ras) <- "uncertainty"
names(Pred_ras) <- "class"
plot(unc_ras)
library(rasterVis)
library(viridis)
library(httr)
png("D:/PROJECTS/04_SocialScape/figures/cluster_prediction.png",
    width=800, height=500, bg="transparent")
plot(Pred_ras, # what to plot
     col=viridis(6), # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
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


# Vizualization ----
gps_sub <- gps_df %>% filter(day==18 )
gps_sub %>% group_by(layer) %>% dplyr::summarise(nind=n_distinct(animals_original_id)) %>% View()
ggplot() +
  geom_point(data = subset(gps_df, day %in% c( 10:20) & layer %in% c(1862)), aes(x = x_, y = y_, color=animals_original_id), size = 5,shape = 20)  + # guides(color=FALSE) + theme_void()
  geom_sf(data = subset(sf_data, layer %in% c(1862)), fill = NA, size=.1) +
  guides(color=FALSE, size=FALSE) + 
  ggtitle("Grid 1862 - day 10 : 20") +
  theme_void()
ggsave("D:/PROJECTS/04_SocialScape/figures/pres_day10_20.png", width = 20, height = 12, units = "cm")

# OLD: function code preparation ----
### Steps 1: definition of a grid
traj <- gps_traj
res <- 200
tt<-SpatialPoints(ld(traj)[,1:2])
tt1<-apply(coordinates(tt), 2, min)
tt2<-apply(coordinates(tt), 2, max)
ras<-raster(xmn=floor(tt1[1]), ymn=floor(tt1[2]),xmx=ceiling(tt2[1]), ymx=ceiling(tt2[2]), res=res)
ras[] <- 1:ncell(ras)
plot(ras)
#text(ras, halo=TRUE, hc='white', size=0.2)
spdf_2 <- as(ras,'SpatialPolygonsDataFrame')
sf_data <- st_as_sf(spdf_2)
st_crs(sf_data) = 32632
# convert gps data to sf object
gps_sf <- st_as_sf(gps_utm, coords = c("x", "y"), crs = 32632, agr = "constant")
gps_sf <- st_join(gps_sf, sf_data, join = st_intersects)
library(lubridate)
gps_df <- as.data.frame(gps_sf) %>% mutate(yr_month=format(as.Date(acquisition_time), "%Y-%m"),
                                           day=yday(acquisition_time))


# metrics calculation ----
# tot_ind in the specific study area
tot_ind <- gps_df %>% distinct(animals_original_id) %>% nrow()
# n_ind_ -----
# = total number of distinct ind having visited a particular grid
n_ind <- gps_df %>% dplyr::select(grid_id=layer,animals_original_id, acquisition_time, age=age_class_code_capture, sex) %>%
  group_by(grid_id) %>%
  dplyr::summarise(n_ind=(n_distinct(animals_original_id)/tot_ind)*100,.groups='keep') %>% arrange(grid_id)
n_ind_empty <- data.frame(grid_id=getValues(ras), n_ind=0)
n_ind_df <- left_join(n_ind_empty, n_ind, by="grid_id") %>% #str()
  mutate(n_ind=replace_na(n_ind.y,0))

r_n_ind <- setValues(ras, n_ind_df$n_ind)
plot(r_n_ind, main="Total individuals having visited a grid cell")

# n_ind_max ----
# = max number of individual observed in a particular grid
# time_max = time at which it took place
# add timegroup variable by grid
DT <- data.table(gps_df)
DT[, datetime := as.POSIXct(date, tz = 'UTC')]
DT2 <- NULL

for (i in unique(DT$grid_id)) {
  #i = 2
  subDT <- DT[grid_id == i]
  #if (length(unique(subDT$animals_original_id))==1) {next}
  group_times(DT = subDT, datetime = 'datetime', threshold = '5 minutes')
  #nrow(subDT)
  DT2 <- rbind(DT2, subDT)
}
n_ind_max <- DT2 %>% dplyr::select(grid_id=layer,animals_original_id, acquisition_time, age=age_class_code_capture, sex, timegroup) %>%
  group_by(grid_id, timegroup) %>%
  dplyr::summarise(time_max=mean(acquisition_time),
                   n_ind=n_distinct(animals_original_id),.groups='keep') %>%
  group_by(grid_id) %>% dplyr::summarise(n_ind_max=(max(n_ind)/tot_ind)*100)#%>% ggplot(.,aes(x=n_ind))  + geom_histogram()

n_ind_max_empty <- data.frame(grid_id=getValues(ras), n_ind_max=0)
d <- left_join(n_ind_max_empty, n_ind_max, by="grid_id") %>% #str()
  mutate(n_ind_max=replace_na(n_ind_max.y,0))

r_n_ind_max <- setValues(ras, n_ind_max_df$n_ind_max)
plot(r_n_ind_max, main="Maximum number of individuals observed simultaneously")
cc# duration ----


freqdur_grp_df <- NULL

for (i in unique(DT2$layer)) {
  # subset study areas
  #i <- 10702
   DT3 <-  DT2 %>% filter(grid_id==i) #%>%
   group_pts(DT = DT3, threshold = 100, id = 'id',
            coords = c('x', 'y'), timegroup = 'timegroup')

   DT4 <- DT3 %>% dplyr::select(grid_id, timegroup, date, id) %>% group_by(timegroup) %>% # arrange(timegroup)
    dplyr::summarise(date=mean(date),
                     grid_id=max(grid_id),
                     n=n_distinct(id),
                     IDs=list(id)) %>% #View() # as_tibble()  %>%
    mutate(IDs2 = map2(IDs, lag(IDs), intersect),
           n2=lengths(IDs2)) %>% data.table() %>% View()
     filter(n>1 | n2 > 1) %>%
     #mutate(visit=case_when(as.numeric(difftime(acquisition_time,lag(acquisition_time), unit="hours"))< 1.1~1,
     #                       as.numeric(difftime(acquisition_time,lag(acquisition_time), unit="hours"))>= 1.1~0),
     #       timegroup1=timegroup) %>%
  data.table()
   if (nrow(DT4)<3) {next}
   DT4[, datetime := as.POSIXct(acquisition_time, tz = 'UTC')]
   group_times(DT = DT4, datetime = 'acquisition_time', threshold = '1 hour')
DT4 %>% arrange(timegroup) %>% View()
   # duration: select where 2 or > ind stay on a grid
   duration <- DT4 %>% filter(n2>1) %>% group_by(timegroup) %>%
      dplyr::summarise(start=min(datetime),
                       end=max(datetime),
                       duration=difftime(end,start, units="hours")) %>%
     dplyr::summarise(meandur=mean(duration),
                      maxdur=max(duration)) %>% mutate(grid_id=i) %>% data.frame()
   # frequency
   frequency <- DT4 %>% dplyr::summarise(visit=n_distinct(timegroup),
                            start=min(datetime),
                            end=max(datetime),
                            duration=difftime(end,start, units="days"),
                            freq_day=visit/as.numeric(duration)) %>% mutate(grid_id=i) %>%
     dplyr::select(grid_id, freq_day)


   # Group composition change
   DT5 <- DT3 %>% dplyr::select(layer, timegroup, acquisition_time, animals_original_id) %>% group_by(timegroup) %>% # arrange(timegroup)
     dplyr::summarise(acquisition_time=mean(acquisition_time),
                      layer=max(layer),
                      n=n_distinct(animals_original_id),
                      IDs=list(animals_original_id)) %>% #View() # as_tibble()  %>%
     mutate(IDs2 = map2(IDs, lag(IDs), intersect),
            n2=lengths(IDs2)) %>%
     filter(n>1 | n2 > 1) %>%
     mutate(IDs2 = map2(IDs, lag(IDs), intersect),
            n2=lengths(IDs2)) %>%
     #mutate(visit=case_when(as.numeric(difftime(acquisition_time,lag(acquisition_time), unit="hours"))< 1.1~1,
     #                       as.numeric(difftime(acquisition_time,lag(acquisition_time), unit="hours"))>= 1.1~0),
     #       timegroup1=timegroup) %>%
     data.table() #%>% group_by(n2) %>% tally()

   grp_comp <- DT5 %>% dplyr::summarise(sum_n2=sum(n2)/n(),
                            mean_n2=mean(n2)/n(),
                            count=n()) %>% mutate(grid_id=i)


   freqdur <- left_join(duration, frequency, by="grid_id")
   freqdur_grp <- left_join(freqdur, grp_comp, by="grid_id")

   freqdur_grp_df <- rbind(freqdur_grp_df, freqdur_grp)
}

empty_df <- data.frame(grid_id=getValues(ras), meandur=0, maxdur=0, freq_day=0, sum_n2=0,  mean_n2=0, count=0)
freqdur_grp_df2 <- left_join(empty_df, freqdur_grp_df, by="grid_id") %>% #str()
  mutate(meandur=replace_na(meandur.y,0),
         maxdur=replace_na(maxdur.y,0),
         freq_day=replace_na(freq_day.y,0),
         sum_n2=replace_na(sum_n2.y,0),
         mean_n2=replace_na(mean_n2.y,0),
         count=replace_na(count.y,0)) %>%
  dplyr::select(grid_id, meandur, maxdur, freq_day, sum_n2, mean_n2, count)

duration_mean <- setValues(ras, as.numeric(freqdur_grp_df2$meandur))
plot(duration_mean, main="Mean duration 2 or > individuals observed simultaneously (hours)")

duration_max <- setValues(ras, as.numeric(freqdur_grp_df2$maxdur))
plot(duration_max, main="Max. duration 2 or > individuals observed simultaneously (hours)")

freq_day <- setValues(ras, as.numeric(freqdur_grp_df2$freq_day))
plot(freq_day, main="Frequency of visit by (any) 2 or > individuals")

sum_n2 <- setValues(ras, as.numeric(freqdur_grp_df2$sum_n2))
plot(sum_n2, main="Average number of similar group member between consecutive visit")

mean_n2 <- setValues(ras, as.numeric(freqdur_grp_df2$mean_n2))
plot(mean_n2, main="")

count_n2 <- setValues(ras, as.numeric(freqdur_grp_df2$count))
plot(count_n2, main="")

rstack <-stack(r_n_ind,r_n_ind_max,duration_mean,duration_max,freq_day, sum_n2, count_n2)
names(rstack) <- c("r_n_ind","r_n_ind_max","duration_mean","duration_max","freq_day", "sum_n2", "count_n2")
plot(rstack, main=c("Total individuals having visited a grid cell",
                    "Maximum number of individuals observed simultaneously",
                    "Mean duration 2 or > individuals observed simultaneously (hours)",
                    "Max. duration 2 or > individuals observed simultaneously (hours)",
                    "Frequency of visit by (any) 2 or > individuals (/day)",
                    "Average number of similar group member between consecutive visit",
                    "Total count of visit by group (i.e. 2 or > ind.)"))



# frequency ----
i <- 3355
DT3 <-  DT2 %>% filter(layer==i) #%>%
group_pts(DT = DT3, threshold = 1000, id = 'animals_original_id',
          coords = c('x_', 'y_'), timegroup = 'timegroup')

DT4 <- DT3 %>% dplyr::select(layer, timegroup, acquisition_time, animals_original_id) %>% group_by(timegroup) %>% # arrange(timegroup)
  dplyr::summarise(acquisition_time=mean(acquisition_time),
                   layer=max(layer),
                   n=n_distinct(animals_original_id),
                   IDs=list(animals_original_id)) %>% #View() # as_tibble()  %>%
  mutate(IDs2 = map2(IDs, lag(IDs), intersect),
         n2=lengths(IDs2)) %>%
  filter(n>1 | n2 > 1) %>%

  #mutate(visit=case_when(as.numeric(difftime(acquisition_time,lag(acquisition_time), unit="hours"))< 1.1~1,
  #                       as.numeric(difftime(acquisition_time,lag(acquisition_time), unit="hours"))>= 1.1~0),
  #       timegroup1=timegroup) %>%
  data.table()
DT4[, datetime := as.POSIXct(acquisition_time, tz = 'UTC')]
group_times(DT = DT4, datetime = 'acquisition_time', threshold = '1 hour')

DT4 %>% dplyr::summarise(visit=n_distinct(timegroup),
                         start=min(datetime),
                         end=max(datetime),
                         duration=difftime(end,start, units="days"),
                         freq=visit/as.numeric(duration),
                         count_0 = )
# Group composition change ----
DT5 <- DT3 %>% dplyr::select(layer, timegroup, acquisition_time, animals_original_id) %>% group_by(timegroup) %>% # arrange(timegroup)
  dplyr::summarise(acquisition_time=mean(acquisition_time),
                   layer=max(layer),
                   n=n_distinct(animals_original_id),
                   IDs=list(animals_original_id)) %>% #View() # as_tibble()  %>%
  mutate(IDs2 = map2(IDs, lag(IDs), intersect),
         n2=lengths(IDs2)) %>%
  filter(n>1 | n2 > 1) %>%
  mutate(IDs2 = map2(IDs, lag(IDs), intersect),
         n2=lengths(IDs2)) %>%
  #mutate(visit=case_when(as.numeric(difftime(acquisition_time,lag(acquisition_time), unit="hours"))< 1.1~1,
  #                       as.numeric(difftime(acquisition_time,lag(acquisition_time), unit="hours"))>= 1.1~0),
  #       timegroup1=timegroup) %>%
  data.table() #%>% group_by(n2) %>% tally()

DT5 %>% dplyr::summarise(sum_n2=sum(n2)/n(),
                         mean_n2=mean(n2)/n(),
                         count=n())



  DT5<-data.table(DT4)
   frequency=nrow(DT5)
   DT5[, visit_id := rleid(visit)]
   DT5 %>% group_by(visit, visit_id) %>% tally()
  #  move to next study areas if no close interaction within dyads
  #if (min(dyads$distance) > 50) {next}
  # add ancillary information
  dyads_s1 <- left_join(dyads, timeloc, by="timegroup") %>% View()
    mutate(together=case_when(distance <= 50~1,
                              distance > 50~0),
           JDate=yday(acquisition_time),
           day=mday(acquisition_time),
           year=year(acquisition_time), month=month(acquisition_time),
           yr_month=format(as.Date(acquisition_time), "%Y-%m"),
           season=time2season(acquisition_time, out.fmt="seasons"),
           dist_class=case_when(distance <= 50~"0-50m",
                                distance > 50 & distance <= 200~"50-200m",
                                distance > 200 & distance <= 500~"200-500m",
                                distance > 500~"> 500m")) %>%
    arrange(dyadID, acquisition_time)

  #remove those dyads never encountering
  dyads_s2 <- dyads_s1 %>%  group_by(dyadID) %>% mutate(tog=sum(together)) %>% filter(tog>0)
  # consider duplicates, i.e. if > 2 animals, when one is leaving (fission), shall only one fission events be count, even if the animal leaving leaves multiple dyads
  # normally all these shall belong to the same group
  # plus the same should apply for fusion events (appear only once in case one ind. join 2 or more)
  #dyads_s3 <- dyads_s2 %>% #filter(grepl("alb_7091-alb_7995",dyadID)) %>% arrange(group) %>% # View()
  #  group_by(group) %>% add_tally(name="nInd") %>% #mutate(nInd= (round(n/2,0)-1)) %>%
  #  arrange(group) %>%
  #  distinct(group, .keep_all = T) %>% arrange(dyadID, acquisition_time)

  # assign a unique id along the fusion/fission event
  dyadsDT <-data.table(dyads_s2)
  dyadsDT[, together_id := rleid(dyadID,together)]

  # count N fusion/fission events per month
  freq_df <- dyadsDT %>%
    group_by(dyadID, yr_month) %>% mutate(join = ifelse(lag(together) == 0 & together != 0, 1, 0),
                                          leave = ifelse(lag(together) == 1 & together != 1, 1, 0)) %>%# View()
    dplyr::summarise(#Sum = sum(together),
      Nfission = sum(leave, na.rm = T),
      Nfusion = sum(join, na.rm = T))
  # duration
  stats_df <- dyadsDT %>%
    group_by(together_id) %>%
    dplyr::summarise(start=min(acquisition_time),
                     end=max(acquisition_time),
                     duration=difftime(end, start, units = "hours"),
                     dist_max=max(distance),
                     dist_min=min(distance),
                     dist_mean=mean(distance),
                     dist_med=median(distance),
                     nloc=n(),.groups='keep') %>% arrange(together_id)
  ID_join <- dyadsDT  %>% distinct(together_id, .keep_all = T) %>% arrange(together_id)
  duration <- left_join(ID_join, stats_df,by="together_id") %>% arrange(dyadID, start)
  # check if start/end of fusion/fission event takes place in the same month
  dur_df <- duration %>% group_by(dyadID, yr_month, together) %>% dplyr::summarise(sumDur=sum(duration),
                                                                                   mean_dist=mean(dist_mean),
                                                                                   max_dist=max(dist_max)) %>%
    pivot_wider(names_from=together,names_prefix="together", values_from=c(sumDur,mean_dist, max_dist)) %>%
    dplyr::select(dyadID, yr_month, Dfission=sumDur_together0, Dfusion=sumDur_together1,
                  MDfission=mean_dist_together0,MDfusion=mean_dist_together1,MxDfission=max_dist_together0,
                  MxDfusion=max_dist_together1)

  freqdur_df <- left_join(freq_df, dur_df, by=c("dyadID", "yr_month"))# %>% View()
  freqdur_df$study_areas <- i
  agesex <- dyadsDT %>% dplyr::select(dyadID, ID1, ID2,sex_ID1, age_ID1, sex_ID2, age_ID2) %>% distinct(dyadID, .keep_all = T)
  freqdur_df <- left_join(freqdur_df, agesex, by="dyadID")


  # HR overlap
  dyad_loop <- dyadsDT %>% dplyr::select(dyadID, ID1, ID2, yr_month) %>% distinct(dyadID, yr_month, .keep_all = T)
  ov2_df <- NULL
  for (j in 1:nrow(dyad_loop)) {

    dat <- DT %>% filter(animals_original_id %in% c(dyad_loop[j]$ID1, dyad_loop[j]$ID2) & yr_month==dyad_loop[j]$yr_month) %>%
      make_track(x, y, id = animals_original_id)
    trast <- make_trast(dat, res = 50)
    dat1 <- dat %>% nest(data = -id) %>%
      mutate(kde = map(data, ~ hr_kde(., trast = trast, level = c(0.5, 0.9))))
    ov2 <- hr_overlap(dat1$kde, type = "ba", labels = dat1$id, which = "one_to_all",
                      conditional = TRUE)
    ov2$yr_month <- dyad_loop[j]$yr_month
    ov2$dyadID <- dyad_loop[j]$dyadID
    ov2_df <- rbind(ov2_df, ov2)

  }

left_join(n_ind, n_ind_max) %>% View()

# Step 3: Clustering -----
# Unsupervised classification ----
table_test <- test %>% dplyr::select(-c(grid_id,dist_sd,dist_median)) %>% na.omit()
str(table_test)
library(mclust)
dataBIC <- mclustBIC(table_test[,c(1:7)])
print(summary(dataBIC))
plot(dataBIC)
ls <- Mclust(table_test[,c(1:7)], G=5, modelNames="VVE")
str(ls)
ls[["parameters"]][["mean"]]
#plot(ls, what = "uncertainty")
ModPred <- predict.Mclust(ls, table_test[,c(1:7)]) # prediction
Pred_ras <- ras # establishing a rediction raster
values(Pred_ras) <- NA # set everything to NA
# set values of prediction raster to corresponding classification according to rowname
values(Pred_ras)[as.numeric(rownames(table_test))] <- as.vector(ModPred$classification)
plot(Pred_ras)
colours <- rainbow(ls$G) # define 7 colours
plot(Pred_ras, # what to plot
     col = colours, # colours for groups
     colNA = "white", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2 # horizontal size of legend
)
table(values(Pred_ras))


mod1dr <- MclustDR(ls)
summary(mod1dr)
plot(mod1dr, what = "pairs")
plot(mod1dr, what = "boundaries", ngrid = 100)

# Supervised classification
library(rpart)
# Train the model
cart <- rpart(as.factor(table_test$n_ind_abs)~., data=table_test, method = 'class', minsplit = 5)
# print(model.class)
# Plot the trained classification tree
plot(cart, uniform=TRUE, main="Classification Tree")
text(cart, cex = 0.8)

#Paremeters associated to each cluster (i.e. center)
ls$parameters$mean
#Proportion of each cluster (how frequent it is spatially)
ls$parameters$pro


