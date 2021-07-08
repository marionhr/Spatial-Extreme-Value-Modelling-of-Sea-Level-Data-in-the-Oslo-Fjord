#!/usr/bin/env Rscript


#save coast data
new_map_data <- fhidata::norway_map_counties[fhidata::norway_map_counties$lat<60 & fhidata::norway_map_counties$long<11.249034 & fhidata::norway_map_counties$long>7.3,]
new_map_data_kyst <- plyr::count(new_map_data, c("lat","long"))
norwegian_coast_data <- new_map_data_kyst[new_map_data_kyst$freq==1,]
ggplot()+
  geom_point(data=norwegian_coast_data, mapping=aes(x=long, y=lat), color="black")+
  scale_fill_manual(values = c("white"))+
  theme_void()+ 
  coord_quickmap()

write.table(norwegian_coast_data, file=paste0(path,"/../Data/norwegian_south_coast_data.txt"))

##############################################################################################

#sorted coast data, for maximized distance between points

#load data
norwegian_coast_data <- read.table(file=paste0(path,"/../Data/norwegian_south_coast_data.txt"))
norwegian_coast_data <- norwegian_coast_data[norwegian_coast_data$long>9.6,]
norwegian_coast_data$location <- seq(1,nrow(norwegian_coast_data))
all_locations <- as.data.frame(rbind(as.matrix(locations[,c(3,4,2,1)]),as.matrix(norwegian_coast_data)))
names(all_locations) <- names(norwegian_coast_data)
all_locations

#find distances:
all_distances <- find_distances(all_locations)
rownames(all_distances) <- seq(1:nrow(all_distances))
colnames(all_distances) <- seq(1:nrow(all_distances))

#sort data by distance
list_of_indexes_added <- c(1:4)
while(length(list_of_indexes_added)<nrow(norwegian_coast_data)+4-1){
  new_index <- which.max(rowMins(all_distances[-c(list_of_indexes_added),list_of_indexes_added],TRUE))
  new_index <- as.numeric(rownames(all_distances[-c(list_of_indexes_added),list_of_indexes_added])[new_index])
  list_of_indexes_added <- c(list_of_indexes_added,new_index)
  if(any(duplicated(list_of_indexes_added))){
    print(which(duplicated(list_of_indexes_added)))
  }
}

#save results
sorted_norwegian_coast_data <- norwegian_coast_data[(list_of_indexes_added[-c(1:4)]-4),]
write.table(sorted_norwegian_coast_data, file=paste0(path,"/../Data/sorted_norwegian_south_coast_data.txt"))
