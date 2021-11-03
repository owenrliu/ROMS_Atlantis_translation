# Alberto Rovellini 11/2/2021
#################################
# This is step 3B of the ROMS workflow
#################################
# This is basically Hem's code

library(tidyverse)
library(rbgm)
library(tidync)
library(ncdf4)
library(RNetCDF)

source("code/make_statevars_functions.R")

atlantis_bgm <- read_bgm('../../data/atlantis/GOA_WGS84_V4_final.bgm')
atlantis_box <- atlantis_bgm$boxes

bgm_polygons <- atlantis_box %>% select(.bx0) %>% distinct() %>% pull() # GOA

#folder with all dat files

avg_files <- list.files(path="../../outputs/short/", pattern="*avg.dat$", recursive = TRUE, full.names = TRUE, include.dirs = TRUE)

# general
fill_value <- 0 # is this appropriate?
this_geometry <- "GOA_WGS84_V4_final.bgm"
this_title <- "ROMS"

# the below is entered manually but it may be automated from the BGM file instead, especially using rbgm
#options depth dimension
depth_bins <- 1:7 # 0 30 70 100 300 500 2969 in GOA. It seems to include sediment in PS
d_units <- "depth layers"

#options time dimension
timestep <- 12 # 12 hour timesteps
t_units <- "seconds since 2017-01-01 00:00:00" # be careful with your start date here
# time.unit.length <- 2 # years
time_length <- 31*2 # how do we handle files of different length here?
seconds_timestep <- 60*60*12
time_vector <- 1:time_length

#options polygon dimension
pol_units <- "spatial polygons"
#pol.bins <- 1:109

time_array <- make_time_array(time_vector)

# pol.array <- pol.bins %>% 
#   as.array

depth_array <- depth_bins %>% 
  as.array()

# this stitches all the dat files into one table with all time steps from ROMS
avg_data <- lapply(avg_files,clean_datafiles) %>% 
  bind_rows() %>% 
  arrange(time_step,polygon,depth_layer)

# this makes a unique index for the time step based on y-m-ts
# avg_data_index <- avg_data %>% 
#   mutate(index=paste(year,month,day, sep="_"))

time_steps <- avg_data %>% 
  distinct(time_step) %>% 
  pull(time_step)

# these_years <- avg_data %>% 
#   distinct(year) %>% 
#   pull(year)

index_list <- avg_data %>% 
  distinct(time_step) %>% 
  pull(time_step)

temp_list_array <- get_array("temperature") # this lists variables from the bottom-up

temp_result_array <- array(dim=c(length(depth_bins),length(bgm_polygons),length(temp_list_array)))

for(eachindex in 1:length(temp_list_array)){
  
  this_vector <- do.call("cbind", temp_list_array[eachindex])
  temp_result_array[,,eachindex] <- this_vector
  
}

salt_list_array <- get_array("salinity")

salt_result_array <- array(dim=c(length(depth_bins),length(bgm_polygons),length(salt_list_array)))

for(eachindex in 1:length(salt_list_array)){
  
  this_vector <- do.call("cbind", salt_list_array[eachindex])
  salt_result_array[,,eachindex] <- this_vector
  
}

make_statevars("temperature", 
               nc_name="goa_roms_temp.nc", 
               t_units, 
               seconds_timestep, 
               this_title, 
               this_geometry, 
               time_array, 
               this_array=temp_result_array, 
               nbox=length(bgm_polygons), 
               nlayer=length(depth_bins))

make_statevars("salinity", 
               nc_name="goa_roms_salt.nc", 
               t_units, seconds_timestep, 
               this_title, 
               this_geometry, 
               time_array, 
               this_array=salt_result_array, 
               nbox=length(bgm_polygons), 
               nlayer=length(depth_bins))
