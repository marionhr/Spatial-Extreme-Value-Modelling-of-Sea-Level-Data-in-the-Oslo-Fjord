#!/usr/bin/env Rscript

u10_brick <- brick(paste0(path,"/../../Data/1950_2020.nc"), varname="u10")
u10_brick
v10_brick <- brick(paste0(path,"/../../Data/1950_2020.nc"), varname="v10")
v10_brick
msl_brick <- brick(paste0(path,"/../../Data/1950_2020.nc"), varname="msl")
msl_brick
tide_raster <- raster(paste0(path,"/../../Data/LAT_Oslofjorden.nc"), varname="height")
tide_raster

#find mean_annual_max/min for entire fields for u10/v10/msl
#extract the data
data_matrix <- rasterToPoints(u10_brick$X2020.12.31.23.00.00)
coordinates <- data_matrix[,1:2]
coordinates <- as.data.frame(coordinates)

data_u10 <- raster::extract(u10_brick, cbind(coordinates$x,
                                             coordinates$y))
data_v10 <- raster::extract(v10_brick, cbind(coordinates$x,
                                             coordinates$y))
data_msl <- raster::extract(msl_brick, cbind(coordinates$x,
                                             coordinates$y))

#data frame
years <- seq(1950,2020)
days_in_year <- as.numeric(strftime(paste0(years,"-",12,"-",31),format="%j"))

data_u10 <- as.data.frame(t(data_u10))
data_u10$year <- rep(years, days_in_year)
data_u10

data_v10 <- as.data.frame(t(data_v10))
data_v10$year <- rep(years, days_in_year)
data_v10

data_msl <- as.data.frame(t(data_msl))
data_msl$year <- rep(years, days_in_year)
data_msl

#study the data
years <- seq(1950,2020)

annual_maximum_u10 <- matrix(NA, nrow=length(years),ncol=66+1)
annual_maximum_v10 <- matrix(NA, nrow=length(years),ncol=66+1)
annual_minimum_msl <- matrix(NA, nrow=length(years),ncol=66+1)

for (i in 1:length(years)){
  annual_maximum_u10[i,] <- colMaxs(as.matrix(data_u10[data_u10$year==years[i],]),value=TRUE)
  annual_maximum_v10[i,] <- colMaxs(as.matrix(data_v10[data_v10$year==years[i],]),value=TRUE)
  annual_minimum_msl[i,] <- colMins(as.matrix(data_msl[data_msl$year==years[i],]),value=TRUE)
}

mean_annual_maximum_u10 <- colMeans(annual_maximum_u10[,1:66])
mean_annual_maximum_u10
u10_data_frame <- data.frame(x=coordinates$x,
                             y=coordinates$y,
                             z=mean_annual_maximum_u10)
u10_raster <- rasterFromXYZ(u10_data_frame, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

mean_annual_maximum_v10 <- colMeans(annual_maximum_v10[,1:66])
mean_annual_maximum_v10
v10_data_frame <- data.frame(x=coordinates$x,
                             y=coordinates$y,
                             z=mean_annual_maximum_v10)
v10_raster <- rasterFromXYZ(v10_data_frame, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
writeRaster(v10_raster, paste0(path,"/../Data/v10.nc"), overwrite=TRUE)

mean_annual_minimum_msl <- colMeans(annual_minimum_msl[,1:66])
mean_annual_minimum_msl
msl_data_frame <- data.frame(x=coordinates$x,
                             y=coordinates$y,
                             z=mean_annual_minimum_msl)
msl_raster <- rasterFromXYZ(msl_data_frame, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
writeRaster(msl_raster, paste0(path,"/../Data/msl.nc"), overwrite=TRUE)

#plot
pdf(paste0(path, "/../Plots/Data/u10.pdf"),width = 5, height = 4)
plot(u10_raster)
points(locations$longitude,locations$latitude)
dev.off()
pdf(paste0(path, "/../Plots/Data/v10.pdf"),width = 5, height = 4)
plot(v10_raster)
points(locations$longitude,locations$latitude)
dev.off()
pdf(paste0(path, "/../Plots/Data/msl.pdf"),width = 5, height = 4)
plot(msl_raster)
points(locations$longitude,locations$latitude)
dev.off()
pdf(paste0(path, "/../Plots/Data/tide.pdf"),width = 5, height = 4)
plot(tide_raster)
points(locations$longitude,locations$latitude)
dev.off()

#normalize and save rasters

u10_data_frame$z <- normalize_data(u10_data_frame$z)
u10_raster <- rasterFromXYZ(u10_data_frame, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
writeRaster(u10_raster, paste0(path,"/../Data/u10.nc"), overwrite=TRUE)

v10_data_frame$z <- normalize_data(v10_data_frame$z)
v10_raster <- rasterFromXYZ(v10_data_frame, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
writeRaster(v10_raster, paste0(path,"/../Data/v10.nc"), overwrite=TRUE)

msl_data_frame$z <- normalize_data(msl_data_frame$z)
msl_raster <- rasterFromXYZ(msl_data_frame, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
writeRaster(msl_raster, paste0(path,"/../Data/msl.nc"), overwrite=TRUE)

tide_data_frame <- as.data.frame(rasterToPoints(tide_raster$Height.of.Mean.Sea.Level..1996.2014..above.Lowest.Astronomical.Tide..LAT.))
names(tide_data_frame) <- c("x","y","z") 
tide_data_frame$z <- normalize_data(tide_data_frame$z)
tide_raster <- rasterFromXYZ(tide_data_frame, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
writeRaster(tide_raster, paste0(path,"/../Data/tide.nc"), overwrite=TRUE)
