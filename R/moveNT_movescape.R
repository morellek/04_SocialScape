library(moveNT)
source('D:/PROJECTS/04_SocialScape/R/functions_MoveNT.R')
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
# # data import ----
load("D:/PROJECTS/02_Contacts/data/raw_data/gps_animals.Rda")
# # subset to one study area for preliminary test
# head(gps_animals)
gps <- gps_animals %>% dplyr::filter(short_name=="Alb")
gps_utm <- mk_track(gps, .x = longitude, .y = latitude, .t =acquisition_time,id=animals_original_id,
                    crs = sp::CRS("+proj=longlat +ellps=WGS84"), all_cols = T) %>%
  transform_coords(sp::CRS("+init=epsg:32632")) %>%
  mutate(x=as.numeric(x_), y=as.numeric(y_), acquisition_time=as.POSIXct(t_),
         animals_original_id=case_when(is.na(animals_original_id)~animals_original_name,
                                       TRUE ~ as.character(animals_original_id)))
# gps_traj <- as.ltraj(xy = gps_utm[,c("x","y")], date = gps_utm$acquisition_time,
#                id = gps_utm$animals_original_id, typeII = TRUE)
# saveRDS(gps_traj, file = "D:/PROJECTS/02_Contacts/materials/Bastille_Rousseau/gps_ltraj.rds")
# 1 Movescape ----
gps_traj <- readRDS(file = "D:/PROJECTS/02_Contacts/materials/Bastille_Rousseau/gps_ltraj.rds")
gps_traj <- gps_traj[id=c("alb_7091","alb_7092","alb_7095") ]
grid<-loop(gps_traj, 200)
plot(grid[[1]])
#plot(gps_traj[],add=T)

# interpolation
grid2 <-interpolation(gps_traj, grid)
plot(grid2[[2]])
mean_weight<-mosaic_network(grid, index=2, sc=T, fun=mean) #Perform mean weight (not-interpolated)
plot(mean_weight)
#data(albatross)
#grid<-loop(albatross, 35000)

table_grid<-table_cluster(gps_traj, grid)
#Showing the first few rows of the table created.
head(table_grid)
#B- Individual-level clustering
#The first step of the analysis is to apply the clustering to each individual. *ind_clust* apply a mixture model to each individual. It is possible to specify the maximum number of clusters (here 8) and also the covariates to use for the clustering, but the function automatically selects the optimal number of clusters (based on BIC). In our case, 2 individuals had 6 clusters, 2 had seven, and 2 had eight clusters. *ls_ind* simply return a list object with each element representing a single individual.
ls_ind<-ind_clust(table_grid, max.n.clust=4)
#Showing the number of individuals with 6, 7, and 8 clusters
#   (i.e. no individual had less than 6 clusters)
table(unlist(lapply(ls_ind, function(x) x$G)))

#C- Population-level clustering
#After performing the individual clustering, a second clustering is applied via *pop_clust*. 
# This second clustering takes the ouptut of *ind_clust* and will identify which individual clusters could be considered as one population-level clusters. 
# The function automatically selects the optimal number of clusters (based on BIC). It is possible for two clusters from the same individual to be in the same population-level cluster. 
# Likewise, it is possible that a population level cluster does not have all individuals. Here, 3 different population clusters were calculated. 
# The second line extract the center (mean) of each cluster which is helpful in interpreting their meaning. 
# The first cluster was heavily used (weight), well connected (degree), and important for connectivity (betweenness), but albatross were moving slowly and not linearly in them. The second cluster was a cluster with intermediate use, not important for connectivity and still with meandering movement. The third cluster was important for connectivty and albatross were  moving fast and linearly in it. We also extract the proportion of each cluster.
pop<-pop_clust(gps_traj, ls_ind)
pop
#Paremeters associated to each cluster (i.e. center)
pop[[1]]$parameters$mean
#Proportion of each cluster (how frequent it is spatially)
pop[[1]]$parameters$pro
pop[[1]]$classification
pop[[1]]$uncertainty
pop[[2]]
#D- Mapping and results export
#After performing the population level cluster, the function  
# *clust_stack* recombines the individual and population level clustering and 
# produce a *stack* object for each individual albatross showing the most likely cluster, 
#and also the probability of observing each cluster (uncertainty) in any given pixel. 
#This individual level data (but which contains the population level clustering) 
# can be used in a regression based analysis as presented in the manuscript or simply mapped. 

# We developed two functions to produce these maps. 
# *pop_stack* generate for each population cluster, rasters showing if at least one individual is using this pixel for this specific cluster. 
# *pop_overl* display for each pixel all potential use observed, 
# for example a pixel having the value *123* will have at least 
# one individual using this pixel as cluster 1,
# another individual using it as 2, 
# and another individual as 3. 
# We show how frequent each combination are using the *table* function. These object can be exported to be used in other software using the *writeRaster* function.

clust_stack<-clust_stack(grid, pop, ls_ind, table_grid)
#Plotting the stack for the first individual
plot(clust_stack[[6]][["Clust"]])
pop_stack<-pop_stack(clust_stack)
#Plotting the stack object. Each raster shows if at least one
# individual is using the pixel as a specific cluster.
plot(pop_stack)
#writeRaster(pop_stack, "Network.tif", format="GTiff",  bylayer=T, suffix="names", overwrite=T)
pop_overl<-pop_overl(clust_stack)
plot(pop_overl)
#Table showing how frequent each pixels is. A value of 123 indicates
# that at least one individual is using the pixel as 1, one individual
# is using it as 2, and one as 3
data.frame(table(values(pop_overl)))
#Plotting of the overlap raster. The legend makes it hard to see,
# but values range from 0 (no animal present) to 123
# (each type of used is observed in the pixel).
plot(pop_overl)
pop_categ <- pop_overl
pop_categ[pop_categ == 3 | pop_categ == 13 | pop_categ == 23 | pop_categ == 34 |
            pop_categ == 123 | pop_categ == 134| pop_categ == 234| pop_categ == 1234] <- 3
pop_categ[pop_categ == 4 | pop_categ == 14 | pop_categ == 24 | pop_categ == 124] <- 2
pop_categ[pop_categ == 1 | pop_categ == 2 | pop_categ == 12] <- 1
plot(pop_categ)
library(usdm)
vif(stack(pop_categ, mean_weight))
## 2 Social grid ----
gps_traj <- readRDS(file = "D:/PROJECTS/02_Contacts/materials/Bastille_Rousseau/gps_ltraj.rds")
### Steps 1: definition of a grid ----
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

# Vizualization
gps_sub <- gps_df %>% filter(day==18 )
gps_sub %>% group_by(layer) %>% dplyr::summarise(nind=n_distinct(animals_original_id)) %>% View()
ggplot() +
  geom_point(data = subset(gps_df, day %in% c( 10:20) & layer %in% c(1862)), aes(x = x_, y = y_, color=animals_original_id), size = 5,shape = 20)  + # guides(color=FALSE) + theme_void()
  geom_sf(data = subset(sf_data, layer %in% c(1862)), fill = NA, size=.1) +
  guides(color=FALSE, size=FALSE) + 
  ggtitle("Grid 1862 - day 10 : 20") +
  theme_void()
ggsave("D:/PROJECTS/04_SocialScape/figures/pres_day10_20.png", width = 20, height = 12, units = "cm")
# Steps 2: metrics calculation ----
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
install.packages("mclust")
library(mclust)
data(diabetes)
class <- diabetes$class
table(class)
## class
## Chemical   Normal    Overt 
##       36       76       33
X <- diabetes[,-1]
class(X)
##   glucose insulin sspg
## 1      80     356  124
## 2      97     289  117
## 3     105     319  143
## 4      90     356  199
## 5      90     323  240
## 6      86     381  157
clPairs(X, class)
?clPairs
library(tidyverse)
stack_df <- as.data.frame(rstack) %>% filter_all(any_vars(. != 0))
saveRDS(stack_df, file = "D:/PROJECTS/04_SocialScape/data/stack_df.rds")
stack_df <- readRDS("D:/PROJECTS/04_SocialScape/data/stack_df.rds")
library(tidyverse)
stack_df <- stack_df %>% filter_all(any_vars(. != 0)) %>% dplyr::select(r_n_ind, r_n_ind_max)
tt <- scale(stack_df)
ind_clust
tt <- scale(table[table_grid$ID == id[i], vars])
ls <- Mclust(tt, G = 2:4)
ls$parameters$mean
#Proportion of each cluster (how frequent it is spatially)
ls$parameters$pro


BIC <- mclustBIC(stack_df)
plot(BIC)
mod1 <- Mclust(stack_df)
summary(mod1, parameters = TRUE)
library(moveNT)
loop
ind_clust


# 3 Extract environmental variables ----
# forest cover, forest type, road density, human footprint

# 4 Application ----
# 4.1 Bialowieza ----
# 4.2 Belgium ----