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
library(mapview)
library(virtualspecies)
library(spdplyr)
library(ENMeval)
library(lubridate)
library(amt)
library(osmextract)
library(rasterVis)
library(viridis)
# data import ---
intensity <- raster("D:/PROJECTS/04_SocialScape/data/derived_data/raster_use_intensity.tif")
corridor <- raster("D:/PROJECTS/04_SocialScape/data/derived_data/raster_corridor.tif")
socdirect <- raster("D:/PROJECTS/04_SocialScape/data/derived_data/raster_social_direct.tif")
socindir <- raster("D:/PROJECTS/04_SocialScape/data/derived_data/raster_social_indirect.tif")


plot(intensity)
arg <- list(at=c(1,2,3,4,5,6), labels=c("Low","Medium","High","High_Slow","High_Medium","High_Fast"))

plot(intensity, col=viridis(3), colNA = "white", legend.shrink=.5, legend.width=1,axis.args=list(at=c(1,2,3), labels=c("Low","Medium","High")),
     box=F, axes=F, main="Intensity of use")
plot(corridor, col=viridis(4), colNA = "white", legend.shrink=.5, legend.width=1,axis.args=list(at=c(2,4,5,6), labels=c("Medium","High_Slow","High_Medium","High_Fast")),
     box=F, axes=F, main="Corridor")
plot(socdirect, col=viridis(3), colNA = "white", legend.shrink=.5, legend.width=1,axis.args=list(at=c(1,2,3), labels=c("Low","Medium","High")),
     box=F, axes=F, main="Direct social interactions")
plot(socindir, col=viridis(3), colNA = "white", legend.shrink=.5, legend.width=1,axis.args=list(at=c(1,2,3), labels=c("Low","Medium","High")),
     box=F, axes=F, main="Indirect social interactions")




# landscape metrics -----
# https://cran.r-project.org/web/packages/landscapemetrics/index.html
# https://www.youtube.com/watch?v=NQZNpyEgVss
library(landscapemetrics)
library(landscapetools)

show_landscape(landcov)

# STACK covariates ----
stack <- stack("D:/cov_stack.grd")
plot(stack, axes=F, box=F)
names(cov_stack)
stack2 <- stack(stack, dForest, dEdges)
names(cov_stack)
writeRaster(stack2, "D:/PROJECTS/04_SocialScape/data/derived_data/cov_stack_bialowieza.grd")
