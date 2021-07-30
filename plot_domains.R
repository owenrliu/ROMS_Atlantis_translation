# Make a map of the ROMS and Atlantis domains

library(sf)
library(here)
library(tidyverse)
library(rbgm)
library(tidync)
library(rnaturalearth)

# Boundary geometry file from Atlantis
atlantis_bgm <- read_bgm(here('data','atlantis','CalCurrentV3_utm.bgm'))
#Atlantis geometry as an sf shapefile
atlantis_sf <- atlantis_bgm %>% box_sf()


# cell lat/lons from ROMS
romsfile <- here::here('data','roms','wc12_ccsra31_his_month_avg_1981_2010_Jan.nc')
roms <- tidync(romsfile)
roms_vars <- hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables asssociated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

latlon_rhogrd <- roms_vars %>% filter(name=="lat_rho") %>% pluck('grd')
roms_rho <- roms %>% activate(latlon_rhogrd) %>% hyper_tibble() %>%
  dplyr::select(lon_rho,lat_rho,xi_rho,eta_rho) %>% 
  st_as_sf(coords=c('lon_rho','lat_rho'),crs=4326) %>%  #make into spatial
  # bounding box
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_transform(st_crs(atlantis_sf))

# continent
wc <- ne_countries(continent="North America",returnclass = 'sf') %>% 
  filter(name %in% c('Canada','United States','Mexico')) %>% 
  st_transform(st_crs(atlantis_sf))

# plot
bbox=st_bbox(atlantis_sf)
p <- ggplot()+
  geom_sf(data=atlantis_sf,fill='blue',alpha=0.6)+
  geom_sf(data=roms_rho,fill='light blue',col=NA)+
  geom_sf(data=wc,fill='gray70')+
  xlim(bbox[1],bbox[2])+ylim(bbox[2],bbox[4])+
  annotate('text',x=0,y=2e6,angle=320,label="atlantis")+
  annotate('text',x=-4e5,y=1.5e6,angle=320,label="ROMS")+
  theme_minimal()
p

ggsave(here('atlantis_ROMS_overlap_map.png'),p,w=6,h=8)
