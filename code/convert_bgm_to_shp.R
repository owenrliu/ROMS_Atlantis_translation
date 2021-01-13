# Convert CalCurrentV3_utm.bgm to a shapefile for viewing later
# Using R package 'rbgm' from Michael Sumner
#https://cran.r-project.org/web/packages/rbgm/
# install.packages('rbgm')
library(rbgm)
library(dplyr)
library(sf)

# file name here of downloaded .bgm file from current Atlantis model
fl <- "C:/Users/Owen.Liu/Documents/NWFSC Research/Packard/Atlantis/CalCurrentV3_utm.bgm"

# load the file
bgm <- read_bgm(fl)
names(bgm)

# can look at some of the info about the bgm's boxes and faces
print(bgm$boxes)
print(bgm$faces)

# the spatial projection info (proj4 string) is in "extra"
bgm$extra

# Let's make a shapefile to save
# convert to a spatial object
emocc_shp <- box_sp(bgm)

# or an sf file
emocc_sf <- box_sf(bgm)
st_crs(emocc_sf) <- st_crs(attr(emocc_sf$geometry, "crs")$proj)

emocc_sf

# save the shapefile
st_write(emocc_sf,"C:/Users/Owen.Liu/Documents/NWFSC Research/Packard/Atlantis/emocc_112720.shp")

# could also look at the "inside" points for each box (and save them if we want)
inside <- bgm$boxes %>% dplyr::select(insideX,insideY)
glimpse(inside)
# make into a spatial file
inside_sf <- inside %>% st_as_sf(coords=c('insideX','insideY'),crs=st_crs(emocc_sf)) %>% 
  # add box number
  mutate(boxNum=bgm$boxes$.bx0)

plot(inside_sf)

# is inside? matrix of inside points (rows) and which polygons they are inside (columns)
# if correct, every row should have exactly 1 "TRUE" (one polygon it is inside)
check_inside <- st_within(inside_sf,emocc_sf,sparse=F)
# do all rows sum to 1 (i.e., one TRUE value)
all(rowSums(check_inside)==1)
# yes

# Plot
library(ggplot2)
emocc_plot <- ggplot()+
  geom_sf(data=emocc_sf,aes(fill=.bx0))+
  geom_sf(data=inside_sf)+
  labs(fill="Box Number",title="EMOCC Boxes and Inside Points")
ggsave("C:/Users/Owen.Liu/Documents/NWFSC Research/Packard/Atlantis/emocc_polys.png",w=6,h=8)
