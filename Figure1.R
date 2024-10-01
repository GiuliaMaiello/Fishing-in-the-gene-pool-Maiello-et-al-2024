### Load libraries 
library("ggmap")
library("ggplot2")
library("ggrepel")
library("ggsn")
library("maps")
library("maptools")
library("marmap")
library("openxlsx")
library("RColorBrewer")

### Load data
mp_Staz =  read.xlsx("Coordinates_Medits_2022.xlsx", sep = ";")
is.num <- sapply(mp_Staz, is.numeric)
mp_Staz[is.num] <- lapply(mp_Staz[is.num], round, 2)

register_stadiamaps(key= "7271dcb5-8253-4f4d-a3c7-6e1ce32244ce")

### GSA 17 - Northern Adriatic Sea
mp_Staz_17 <- mp_Staz[which(mp_Staz$GSA %in% 'Adr'),]
map17 <- get_stadiamap(bbox = c(left = 11.7,
                                bottom = 43.7,
                                right =  14.2,
                                top = 45.8),
                      zoom = 10,
                      maptype = "stamen_terrain_background")
minlon = 11.7
minlat = 43.7
maxlon = 14.2
maxlat = 45.8

pt.baty = getNOAA.bathy(lon1 = minlon, lon2 =  maxlon,
                        lat1 = minlat, lat2 = maxlat, resolution = 1)
pt.baty = -1*pt.baty

ggmap17 = ggmap(map17) +
  geom_contour_filled(data = pt.baty,
                      aes(x=x, y=y, z=z),
                      breaks=c(0, 25, 50, 100),
                      size=c(0.3), alpha = 0.7,
                      colour="grey",
                      show.legend = TRUE) +
  scale_fill_manual("Depth strata", values = brewer.pal(8, "Blues")) +
  geom_point(data = mp_Staz_17, aes(x = Longitude, y = Latitude), size = 2) + 
  geom_text_repel(data = mp_Staz_17, aes(x = Longitude, y = Latitude,label=Staz),size = 5)+
  xlab("") +
  ylab("") + 
  theme(legend.position = "right", legend.key.size = unit(1,"cm"))
ggmap17 <- ggmap17 + theme(legend.position = "none")

### GSA 9 -  Northern Tyrrhenian Sea 
mp_Staz_9 <- mp_Staz[which(mp_Staz$GSA %in% 'Tyrr'),]
map9 <- get_stadiamap(bbox = c(left = min(mp_Staz_9$Longitude)-0.1,
                               bottom = min(mp_Staz_9$Latitude)-0.2,
                               right =  max(mp_Staz_9$Longitude)+0.1,
                               top = max(mp_Staz_9$Latitude)+0.2),
                      zoom = 10,
                      maptype = "stamen_terrain_background")
minlon = 9.85
minlat = 40.75
maxlon = 13.53
maxlat = 43.05

pt.baty = getNOAA.bathy(lon1 = minlon, lon2 =  maxlon,
                        lat1 = minlat, lat2 = maxlat, resolution = 1)
pt.baty = -1*pt.baty

ggmap09 = ggmap(map9) +
  geom_contour_filled(data = pt.baty,
                      aes(x=x, y=y, z=z),
                      breaks=c(0,25, 50, 100, 200, 500, 1000,2000,4000),
                      size=c(0.3), alpha = 0.7,
                      colour="grey",
                      show.legend = TRUE) +
  scale_fill_manual("Depth strata", values = brewer.pal(8, "Blues")) +
  geom_point(data = mp_Staz_9, aes(x = Longitude, y = Latitude), size = 2) + 
  geom_text_repel(data = mp_Staz_9, aes(x = Longitude, y = Latitude,label=Staz), size = 5)+
  xlab("") +
  ylab("") + 
  theme(legend.position = "right", legend.key.size = unit(1,"cm"))
ggmap09 <- ggmap09 + theme(legend.position = "none")

### GSA 11 - Sardinian Sea
mp_Staz_11 <- mp_Staz[which(mp_Staz$GSA %in% 'Sard'),]
map11 <- get_stadiamap(bbox = c(left = min(mp_Staz_11$Longitude)-0.4,
                               bottom = min(mp_Staz_11$Latitude)-0.3,
                               right =  max(mp_Staz_11$Longitude)+0.8,
                               top = max(mp_Staz_11$Latitude)+0.1),
                      zoom = 10,
                      maptype = "stamen_terrain_background")
minlon = 7.57
minlat = 38.35
maxlon = 9.21
maxlat = 39.64

pt.baty = getNOAA.bathy(lon1 = minlon, lon2 =  maxlon,
                        lat1 = minlat, lat2 = maxlat, resolution = 1)
pt.baty = -1*pt.baty

ggmap11 = ggmap(map11) +
  geom_contour_filled(data = pt.baty,
                      aes(x=x, y=y, z=z),
                      breaks=c(0,25, 50, 100, 200, 500, 1000,2000,4000),
                      size=c(0.3), alpha = 0.7,
                      colour="grey",
                      show.legend = TRUE) +
  scale_fill_manual("Depth strata", values = brewer.pal(8, "Blues")) +
  geom_point(data = mp_Staz_11, aes(x = Longitude, y = Latitude), size = 2) + 
  geom_text_repel(data = mp_Staz_11, aes(x = Longitude, y = Latitude,label=Staz), size = 5)+
  xlab("") +
  ylab("") + 
  theme(legend.position = "right", legend.key.size = unit(1,"cm"))
ggmap11 <- ggmap11 + theme(legend.position = "none")
