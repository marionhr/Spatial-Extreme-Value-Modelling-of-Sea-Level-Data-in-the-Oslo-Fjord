#!/usr/bin/env Rscript


#coast data

new_map_data <- fhidata::norway_map_counties[fhidata::norway_map_counties$lat<60 & fhidata::norway_map_counties$long<11.249034 & fhidata::norway_map_counties$long>7.3,]
new_map_data_kyst <- plyr::count(new_map_data, c("lat","long"))
norwegian_coast_data <- new_map_data_kyst[new_map_data_kyst$freq==1,]
ggplot()+
  geom_point(data=norwegian_coast_data, mapping=aes(x=long, y=lat), color="black")+
  scale_fill_manual(values = c("white"))+
  theme_void()+ 
  coord_quickmap()

write.table(norwegian_coast_data, file=paste0(path,"/../Data/norwegian_south_coast_data.txt"))

norwegian_coast_data <- read.table(file=paste0(path,"/../Data/norwegian_south_coast_data.txt"))

N_samples <- 30
which_samples <- sample(1:nrow(norwegian_coast_data), N_samples)
samples <- norwegian_coast_data[which_samples,]

ggplot()+
  geom_point(data=norwegian_coast_data, mapping=aes(x=long, y=lat), color="black")+
  geom_point(data=samples, mapping=aes(x=long, y=lat), color="green")+
  scale_fill_manual(values = c("white"))+
  theme_void()+ 
  coord_quickmap()

N_train_samples <- 20
which_train_samples <- sample(1:nrow(samples), N_train_samples)
train_samples <- samples[which_train_samples,]
test_samples <- samples[-which_train_samples,]

ggplot()+
  geom_point(data=norwegian_coast_data, mapping=aes(x=long, y=lat), color="black")+
  geom_point(data=train_samples, mapping=aes(x=long, y=lat), color="green")+
  geom_point(data=test_samples, mapping=aes(x=long, y=lat), color="red")+
  scale_fill_manual(values = c("white"))+
  theme_void()+ 
  coord_quickmap()


#sorted coast data, for maximized distance between points
norwegian_coast_data <- read.table(file=paste0(path,"/../Data/norwegian_south_coast_data.txt"))
norwegian_coast_data <- norwegian_coast_data[norwegian_coast_data$long>9.6,]
norwegian_coast_data$location <- seq(1,nrow(norwegian_coast_data))
all_locations <- as.data.frame(rbind(as.matrix(locations[,c(3,4,2,1)]),as.matrix(norwegian_coast_data)))
names(all_locations) <- names(norwegian_coast_data)
all_locations

all_locations.sf = st_as_sf(all_locations,coords = c("long","lat"),
                            crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0")
transformed.all_locations.sf <- st_transform(all_locations.sf, crs = CRS("+proj=utm +zone=32 +datum=WGS84"))
all_distances <- as.matrix(st_distance(transformed.all_locations.sf))
all_distances  # unit m
rownames(all_distances) <- seq(1:nrow(all_distances))
colnames(all_distances) <- seq(1:nrow(all_distances))
all_distances <- all_distances/(10^3)  # unit km
units(all_distances) <- NULL

which.max(rowMins(all_distances[-c(list_of_indexes_added),list_of_indexes_added],TRUE))
list_of_indexes_added <- c(1:4)
while(length(list_of_indexes_added)<nrow(norwegian_coast_data)+4){
  new_index <- which.max(rowMins(all_distances[-c(list_of_indexes_added),list_of_indexes_added],TRUE))
  new_index <- as.numeric(rownames(all_distances[-c(list_of_indexes_added),list_of_indexes_added])[new_index])
  list_of_indexes_added <- c(list_of_indexes_added,new_index)
  if(any(duplicated(list_of_indexes_added))){
    print(which(duplicated(list_of_indexes_added)))
  }
}
list_of_indexes_added
length(list_of_indexes_added)
nrow(norwegian_coast_data)
nrow(sorted_norwegian_coast_data)

sorted_norwegian_coast_data <- norwegian_coast_data[(list_of_indexes_added[-c(1:4)]-4),]
write.table(sorted_norwegian_coast_data, file=paste0(path,"/../Data/sorted_norwegian_south_coast_data.txt"))
