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
    # mutate(day = day(time_step),
    #        month = month(time_step),
    #        year = year(time_step))
  
  return(this_data)
}


get_array <- function(eachvariable) {
  
  print(eachvariable)
  
  variable_list <- list()
  
  for(eachindex in index_list)  {
    
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
