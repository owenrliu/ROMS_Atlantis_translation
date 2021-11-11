# Alberto Rovellini 
# 11/2/2021
#################################
# Step 4: from temp and salt DAT files to temp.nc and salt.nc
#################################
# Based on Hem's code

library(tidyverse)
library(rbgm)
library(tidync)
library(ncdf4)
library(RNetCDF)

select <- dplyr::select

# source("code/make_statevars_functions.R") # figure out scoping

atlantis_bgm <- read_bgm('../../data/atlantis/GOA_WGS84_V4_final.bgm')
atlantis_box <- atlantis_bgm$boxes

roms_files <- list.files('../../outputs/2017/monthly/post-interp/',full.names = TRUE)
avg_files <- roms_files[grepl('avg',roms_files)]

make_statevars_nc <- function(avg_file, atlantis_sf=atlantis_box){
  
  #########################################################################################################
  # TODO: move this!!!
  # functions called in make_statevars.
  # this is Hem's code.
  # consider having a similar setup for the other steps - with one or several function files (pros and cons) 
  
  isEmpty <- function(x) {
    return(identical(x, numeric(0)))
  }
  
  #create array for time
  time_list <- function(eachtimestep){
    
    this_step <- eachtimestep * seconds_timestep
    return(this_step)
    
  }
  
  make_time_array <- function(time_vector) {
    
    time_array <- lapply(time_vector,time_list) %>% 
      unlist() %>% 
      as.array()
    
    return(time_array)
    
  }
  
  #read input file
  clean_datafiles <- function(this_file){
    
    this_data <- read.csv(this_file,sep="\t",header=TRUE) %>% 
      select(-layer) %>%
      setNames(c("time_step","polygon","depth_layer","vertical_velocity","temperature","salinity")) #%>% 
    
    return(this_data)
  }
  
  
  get_array <- function(eachvariable) {
    
    print(eachvariable)
    
    variable_list <- list()
    
    for(eachindex in time_steps)  {
      
      print(eachindex)
      this_index_data <- avg_data %>% 
        filter(time_step==eachindex)
      
      pol_list <- list()
      
      for(eachpolygon in bgm_polygons) {
        
        this_value <- this_index_data %>% 
          filter(polygon == eachpolygon) %>% 
          arrange(desc(depth_layer)) %>% 
          select(all_of(eachvariable)) %>% 
          pull(1)
        
        if(isEmpty(this_value)){
          
          length(this_value) <- 7
          
        } else if(!isEmpty(this_value)){
          
          length(this_value) <- 7
          
          this_value[[7]] <- this_value[[1]]
          
        }
        
        pol_array <- this_value %>% 
          as.array()
        
        pol_list[[eachpolygon+1]] <- pol_array
        
      }
      
      this_index <- do.call("cbind", pol_list) %>% as.array()
      
      variable_list[[eachindex]] <- this_index
      
    }
    
    return(variable_list)
    
  }
  
  # create and write the netCDF file -- ncdf4 version
  # define dimensions
  
  make_statevars <- function(eachvariable, nc_name, t_units, seconds_timestep, this_title, this_geometry, time_array, this_array, nbox, nlayer) {
    
    nc_file <- create.nc(nc_name)
    
    dim.def.nc(nc_file, "t", unlim=TRUE)
    dim.def.nc(nc_file, "b", nbox) # manual 
    dim.def.nc(nc_file, "z", nlayer) # manual 
    
    var.def.nc(nc_file, "t", "NC_DOUBLE", "t")
    var.def.nc(nc_file, eachvariable, "NC_DOUBLE", c("z","b","t"))
    
    att.put.nc(nc_file, eachvariable, "_FillValue", "NC_DOUBLE", 0)
    att.put.nc(nc_file, "t", "units", "NC_CHAR", t_units)
    att.put.nc(nc_file, "t", "dt", "NC_DOUBLE", seconds_timestep)
    att.put.nc(nc_file, "NC_GLOBAL", "title", "NC_CHAR", this_title)
    att.put.nc(nc_file, "NC_GLOBAL", "geometry", "NC_CHAR", this_geometry)
    att.put.nc(nc_file, "NC_GLOBAL", "parameters", "NC_CHAR", "")
    
    var.put.nc(nc_file, "t", time_array)
    var.put.nc(nc_file, eachvariable, this_array)
    close.nc(nc_file)
  }
  #########################################################################################################
  
  # read in data
  # this stitches all the dat files into one table with all time steps from ROMS
  avg_data <- clean_datafiles(avg_file) %>% 
    arrange(time_step,polygon,depth_layer)
  
  monthyear <- sub('.*post-interp/ *(.*?) *_avg.*','\\1',avg_file) # get the month-year index
  
  bgm_polygons <- atlantis_box %>% select(.bx0) %>% distinct() %>% pull() # GOA
  
  # general
  # TODO: set these as args of make_statevars_nc? Or global vars?
  fill_value <- 0 # is this appropriate?
  this_geometry <- "GOA_WGS84_V4_final.bgm"
  this_title <- "ROMS"
  
  # the below is entered manually but it may be automated from the BGM file instead, especially using rbgm
  #options depth dimension
  depth_bins <- 1:7 # 0 30 70 100 300 500 2969 in GOA. It seems to include sediment in PS
  d_units <- "depth layers"
  
  #options time dimension
  time_steps <- avg_data %>% select(time_step) %>% distinct() %>% pull()
  time_length <- length(time_steps)
  time_vector <- 1:time_length
  t_units <- paste("seconds since",time_steps[1],sep=" ")
  timestep <- 12 # 12 hour timesteps
  seconds_timestep <- 60*60*12
  
  #options polygon dimension
  pol_units <- "spatial polygons"
  #pol.bins <- 1:109
  
  time_array <- make_time_array(time_vector)
  
  # pol.array <- pol.bins %>% 
  #   as.array
  
  depth_array <- depth_bins %>% 
    as.array()
  
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
                 nc_name=paste0("../../outputs/2017/monthly/forcings/goa_roms_temp_",monthyear,".nc"), 
                 t_units, 
                 seconds_timestep, 
                 this_title, 
                 this_geometry, 
                 time_array, 
                 this_array=temp_result_array, 
                 nbox=length(bgm_polygons), 
                 nlayer=length(depth_bins))
  
  make_statevars("salinity", 
                 nc_name=paste0("../../outputs/2017/monthly/forcings/goa_roms_salt_",monthyear,".nc"), 
                 t_units, 
                 seconds_timestep, 
                 this_title, 
                 this_geometry, 
                 time_array, 
                 this_array=salt_result_array, 
                 nbox=length(bgm_polygons), 
                 nlayer=length(depth_bins))
  
  
}

lapply(avg_files, make_statevars_nc)