# Alberto Rovellini 
# 3/31/2022
#################################
# Step 3: from DAT exchange files to hydro.nc
#################################
# # 
# This code reworks the exchange file from DAT format prepared for HC to NetCDF files to force Atlantis with. 
# It also needs the average file for the state variables because that one contains $w$. We will need to somehow account for $w$.
# 
# Steps:
#   
# 1. Read data from the transformation code, interpolated at 12h time steps if needed
# 2. Fix the sign of the horizontal exchanges. The sign of the exchanges in the transformation code was evaluated by face with the LR convention (negative LR, positive RL). Here we move to source and destination boxes, so changing the sign depending on whether the exchange leaves (negative) or enters (positive) the source box. 
# 3. List each horizontal flux only once. The transformation code listed all fluxes twice, once LR and once RL for a face. Atlantis needs us to list fluxes from a source to a destination cell only once per time step.
# 4. List the vertical fluxes between layers of the same box. Signs: ROMS convention is that negative fluxes go upwards. In the translation code, $w$ is calculated through the bottom of a cell. See below for how this was considered in setting the sign of the vertical exchange from a source cell to a destination cell above it.
# 5. All bottom layers (layer 0) have a vertical flux through the bottom. To simulate advection from outside the model domain, we use $C_{0,0}$ as source of these fluxes. 
# 6. Sum up the three components (horizontal fluxes, vertical fluxes, and vertical fluxes from outside the model domain) to obtain net flux from each source cell to each destination cell. 
# 7. Create three arrays: one for exchanges, one for destination boxes (dest_b), one for destination layers (dest_k).
# 8. Pack this information in a NetCDF file. 
# 

# NOTE 11/04/2021. Fluxes in Atlantis seem to follow the convention that they are positive when they leave a box.
# Flip the sign of the horizontal fluxes, and map the vertical exchanges accordingly. 
# Mark every change related to this with 11/04/2021, in case we need to find this again and change

# NOTE 11/08/2022
# Added new option to correct for hyperdiffusion
# Now setting hyperdiff=2 in the function calls divides the flux across each face by the length across the box that receives the flux
# This accounts for box shape when scaling the fluxes
# in addition, vertical fluxes are divided by the dz of the layer that receives the flux, following the same logic
# This will cause fluxes of larger magnitude across the board
# proportions across fluxes will also be different now, especially in long and skinnny boxes
# this may result in overall different patterns of circulation and is not equivalent to scaling all fluxes the same way across the board
# check that the dominant fluxes are still E-W. This may no longer be the case for long boxes oriented E-W

library(tidyverse)
library(rbgm)
library(tidync)
library(sf)
library(ncdf4)
library(RNetCDF)

select <- dplyr::select

# Read data
#make lists of the avg and hydro files, post-interpolation

roms_files <- list.files('../../outputs/complete_from_Emily/Script1/2017/monthly/post-interp/',full.names = TRUE)
hydro_files <- roms_files[grepl('flux',roms_files)]
avg_files <- roms_files[grepl('avg',roms_files)]

# read in depth lookup key
depth_key <- read.csv('../../outputs/complete_from_Emily/Script1/2017/depth_layer.csv') # TODO: bypass this

# read in bgm
atlantis_bgm <- read_bgm('C:/Users/Alberto Rovellini/Documents/GOA/ROMS/data/atlantis/GOA_WGS84_V4_final.bgm')
atlantis_box <- atlantis_bgm %>% box_sf()

# read in hyperdiffusion correction table
hd_by_length <- read.csv('hd_by_length.csv')

# # make folder
# dir.create('../../outputs/complete_from_Emily/Script1/2017/monthly/forcings')
# dir.create('../../outputs/complete_from_Emily/Script1/2017/monthly/forcings/hydro')


#Extract faces of each polygon from Atlantis and organize them
# information about each face, including its angular coords, and which boxes are to its left and right
faces <- atlantis_bgm$faces %>% select(-label)

# function to apply the whole routine
# TODO: move functions to a different file

make_exchange_nc <- function(hydro_file, atlantis_sf=atlantis_box, faces_tmp=faces, depth=depth_key, fluxperts=F, hyperdiff=0, hd=hd_by_length){
  
  monthyear <- sub('.*post-interp/ *(.*?) *_flux.*','\\1',hydro_file) # get the month-year index
  
  avg_file <- avg_files[grepl(paste0('/',monthyear),avg_files)]
  
  hydro <- read.table(hydro_file,header=TRUE,sep='\t')
  avg <- read.table(avg_file,header=TRUE,sep='\t') # need this for the vertical fluxes
  
  # Set up some dimensions for later.
  layers <- hydro %>% select(Depth.Layer) %>% distinct() %>% pull()
  boxes <- atlantis_sf %>% select(.bx0) %>% distinct()  %>% pull() #hydro %>% select(Polygon.number) %>% distinct() %>% pull() # NOTE! You want to get the islands in the data!!!
  time <- hydro %>% select(Time.Step..12.hr) %>% distinct() %>% pull()
  
  nlayer <- length(layers)
  nbox <- length(boxes)
  ntime <- length(time)
  
  #Construct a table with maximum and minimum depth for each depth layer of each box.
  depth[depth=='_']<- NA
  
  depth <- depth %>% select(-layer6) %>% # ditch sediment - do not think it is needed in the fluxes
    mutate_if(is.character, as.numeric) %>% 
    pivot_longer(cols = starts_with('layer'), names_to = 'layer', values_to = 'depth') %>%
    mutate(layer = as.numeric(str_replace(layer,'layer',''))) %>%
    drop_na() %>%
    arrange(box_id,desc(layer)) %>%
    group_by(box_id) %>%
    mutate(lower = cumsum(depth),
           upper = lag(lower,default = 0),
           dz = lower-upper)
  
  depth_all <- expand.grid(box_id=unique(depth$box_id),layer=unique(depth$layer)) %>%
    left_join(depth,by=c('box_id','layer'))
  
  #####################################################################################
  # Make dest_b, dest_k, and exchange arrays for NetCDF
  
  # This must include the correction for hyperdiffusion, based on box area and shape (e.g. EW length vs NS length).
  # 
  
  hydro <- hydro %>% left_join(faces_tmp %>% select(.fx0,left,right), by = c('FaceID'='.fx0')) # sometimes the source box is on the right, sometimes on the left
  
  hydro <- hydro %>% mutate(adjacent.box = case_when(Polygon.number==left ~ right, # determine which one is the destination box
                                                     Polygon.number==right ~ left)) %>%
    select(-left,-right)
  
  # 11/04/2021: flip sign: a flux going out of a source box is interpreted as positive by Atlantis
  hydro <- hydro %>% mutate(Flux..m3.s.=-1*Flux..m3.s.)
  
  # This step should leave us with half as many exchanges.
  hydro_with_depth <- hydro %>% left_join(depth_all %>% select(box_id,upper,layer),
                                          by=c('Polygon.number'='box_id','Depth.Layer'='layer')) # layer of the source box
  
  # keep both fluxes, but set the second exchange to 0. This way we keep the source-destination interaction, but there is no water getting exchanged
  # if the upper depth is NA, keep the (empty) fluxes as they are because they are between empty layers but need to be there
  # if there is one flux only across a face, it means that it's a layer that is not shared by both boxes (e.g. one box is deeper) and the flux out of it is 0.
  handle_duplicate_flows <- function(upper,data){
    if(!is.na(upper) & nrow(data)>1) data[2,]$Flux..m3.s.<-0 else data
    return(data)
  }
  
  hydro_one_flux <- hydro_with_depth %>% 
    group_by(Time.Step..12.hr,FaceID,upper) %>% 
    nest() %>%
    mutate(OneFlux = purrr::map2(upper,data,handle_duplicate_flows)) %>%
    select(-data) %>%
    unnest(cols = OneFlux)
  
  # update 11/08/2022
  # adding option to divide by section of the receiving box for horizontal fluxes, or by dz for vertical fluxes
  # remember that layer 0 is the bottom
  # this will result in far larger fluxes, but directions should be the same
  # if you want to do it by face you need to do it here, before integrating fluxes across faces into source/destination boxes
  if(hyperdiff==2){
    hydro_one_flux <- hydro_one_flux %>%
      left_join(hd, by = c('FaceID'='.fx0')) %>%
      rowwise() %>%
      mutate(Flux..m3.s. = ifelse(adjacent.box==rbox, Flux..m3.s./distances_r, Flux..m3.s./distances_l)) %>%
      ungroup() %>%
      select(-(rbox:distances_l))
  }
  
  ## Horizontal fluxes
  # sum up fluxes bewteen the same source and destination but across different faces
  dest_b_horizontal <- hydro_one_flux %>% 
    group_by(Time.Step..12.hr,Polygon.number,Depth.Layer,upper,adjacent.box) %>%
    summarise(Exchange = sum(Flux..m3.s.,na.rm=TRUE)) %>% 
    ungroup()
  
  # now need to introduce the depth layer of the destination box that matches the depth layer of the source box.
  # FIXME: can we preserve this information from the R code instead?
  # If we just kept consistent layering from the surface down we could probably avoid all this
  map_dest_k <- function(dest_box,source_upper){
    dest_k <- depth_all %>% filter(box_id==dest_box,upper==source_upper) %>% pull(layer)
    dest_k <- ifelse(identical(dest_k, numeric(0)),NA,dest_k) # we may want to remove this for speed and change after
    return(dest_k)
  }
  
  dest_horizontal_tmp <- dest_b_horizontal %>% # this step takes a long time, ~15-20 minutes per month
    mutate(dest_k=purrr::map2(adjacent.box,upper,map_dest_k)) %>%
    unnest(cols = c(dest_k))
  
  # rename the columns to understand better
  dest_horizontal <- dest_horizontal_tmp %>% 
    select(Time.Step..12.hr,Polygon.number,Depth.Layer,adjacent.box,dest_k,Exchange) %>%
    set_names('ts','source_b','source_k','dest_b','dest_k','exchange') 
  
  if(dest_horizontal %>% filter(is.na(dest_k)) %>% pull(exchange) %>% sum() != 0) stop("There are non-zero exchanges to a non-existing dest_k")
  
  # For now: set the NA layers to 0. This will not matter when we add up the fluxes to the destinations
  dest_horizontal[is.na(dest_horizontal)]<-0
  
  ## Vertical fluxes
  vert <- avg %>% select(Time.Step:Vertical.velocity..m3.s.,layer)
  
  dest_vertical <- vert %>% 
    select(Time.Step,Polygon.number,layer,Vertical.velocity..m3.s.) %>% # all source box
    arrange(Time.Step,Polygon.number,layer) %>%
    set_names(c('ts','source_b','source_k','exchange_vert')) %>%
    mutate(dest_b=source_b) %>% # these are vertical fluxes within a box, so source_b and dest_b are the same
    group_by(ts,source_b) %>%
    mutate(dest_k=lead(source_k,default=0), # note: default=0 here denotes a flux from layer 5 of a box (the surface) to layer 0 of the same box (bottom). This does not exist and it serves the purpose of handling NA's for now - it works because w out of layer 5 of all boxes is NA too so no impact on actual flux calcs
           #exchange = lead(exchange_vert,default=0)) %>% # these are all flows from focal box to box above. A negative value means an upward flux, which means a negative flux for the focal cell ("give"); a positive value indicates a downward flux (because of ROMS), so a flux into the focal cell from the cell above. There should be no flow out of the top layer.
           exchange = -1*lead(exchange_vert,default=0)) %>% # 11/04/2021 these are all flows from focal box to box above. A negative value means an upward flux in ROMS, which means a positive flux from the focal cell ("give"); 
    select(-exchange_vert)
  
  # 11/08/2022 correcting for hyperdiffusion
  # for vertical fluxes we divide by dz. Note that this will cause the system to have proportionally stronger vertical fluxes than horizontal ones
  if(hyperdiff==2){
    dest_vertical <- dest_vertical %>%
      left_join(depth_all %>% select(box_id, layer, dz), by = c('source_b'='box_id', 'source_k'='layer')) %>% # source_b and dest_b are the same box for vertical flows
      left_join(depth_all %>% select(box_id, layer, dz), by = c('source_b'='box_id', 'dest_k'='layer')) %>%
      rename(dz_source = dz.x, dz_dest = dz.y) %>%
      mutate(exchange = case_when(is.na(exchange) ~ exchange,
                                  exchange == 0 ~ 0,
                                  exchange > 0 ~ exchange / dz_dest,
                                  exchange < 0 ~ exchange / dz_source)) %>%
      select(-dz_source, -dz_dest)
  }
  
  ## Flows from outside the model domain through the bottom of Atlantis
  bottom_flows <- vert %>% filter(layer==0) # fluxes through the bottom of all 0 layers (deepest in each box)
  
  # treat Polygon.number as source, 0,0 as dest, and make sure that the sign makes sense:
  # from ROMS, flows up (into 0 through its bottom) are negative, and flows down (out of 0 through its bottom) are positive
  # this should align with how it is for hydro.nc (positive flux leaves the box)
  
  dest_bottom <- bottom_flows %>% 
    mutate(dest_b = 0, dest_k = 0) %>%
    select(Time.Step,Polygon.number,layer,dest_b,dest_k,Vertical.velocity..m3.s.) %>%
    set_names(c('ts','source_b','source_k','dest_b','dest_k','exchange'))
  
  # 11/08/2022 correction for hyperdiffusion
  # dividing by 46 m if the flux enters box 0 layer 0, else by the receiving dz
  if(hyperdiff == 2){
    dest_bottom <- dest_bottom %>%
      left_join(depth_all %>% select(box_id, layer, dz), by = c('source_b'='box_id', 'source_k'='layer')) %>%
      rename(dz_source = dz) %>%
      mutate(dz_dest = 46) %>% # because layer 0 of box 0 is 46 m
      mutate(exchange = case_when(is.na(exchange) ~ exchange,
                                  exchange == 0 ~ 0,
                                  exchange > 0 ~ exchange / dz_source,
                                  exchange < 0 ~ exchange / dz_dest)) %>%
      select(-dz_source, -dz_dest)
  }
  
  #####################################################################################
  ## Bring it all together
  dest_all <- rbind(dest_horizontal,dest_vertical,dest_bottom) %>%
    arrange(ts,source_b,source_k,dest_b,dest_k)
  
  # correct for hyperdiffusion
  # 4/8/2022 using the area of the box of arrival of the flux:
  # if the flux is negative (incoming) divide by area of source box,
  # if the flux is positive (outgoing) divide by area of dest box.
  if(hyperdiff>0){
    if(hyperdiff==1){
      
      dest_all <- dest_all %>%
        left_join((atlantis_sf %>% st_set_geometry(NULL) %>% select(box_id, area)), by = c('source_b'='box_id')) %>%
        left_join((atlantis_sf %>% st_set_geometry(NULL) %>% select(box_id, area)), by = c('dest_b'='box_id')) %>%
        rename(area_source = area.x, area_dest = area.y) %>%
        rowwise() %>%
        mutate(exchange = ifelse(exchange<0, exchange/area_source, exchange/area_dest)) %>%
        ungroup() %>%
        select(-area_source,-area_dest)
        
    } else if(hyperdiff==2) {
      
      print('Hyperdiffusion correction done on individual faces')
      
    } else {
      
      stop('This method of hyperdiffusion correction has not been implemented yet')
      
    }
  }
  
  # set NA exchanges to 0
  dest_all[is.na(dest_all)]<-0
  
  #sum fluxes from the same source to the same destination
  #for example 0,0 talks to 4,0 and 6,0 through horizontal and vertical exchange
  all_fluxes <- dest_all %>% 
    group_by(ts,source_b,source_k,dest_b,dest_k) %>%
    summarise(exchange_tot = sum(exchange,na.rm = TRUE)) %>%
    ungroup()
  
  # Now reshape this horizontally 
  # 1.
  all_fluxes1 <- all_fluxes %>% 
    group_by(ts,source_b,source_k) %>%
    nest() %>% 
    mutate(dest_b = purrr::map(data, ~t(.x$dest_b)),
           dest_k = purrr::map(data, ~t(.x$dest_k)),
           exchange = purrr::map(data, ~t(.x$exchange_tot))) %>%
    select(-data)
  
  if(nrow(all_fluxes1)-(ntime*nbox*nlayer)!=0) stop('The rows of the CDF file are not time*box*layer')
  
  # 2. get the longest vector
  ndest <- all_fluxes1 %>% mutate(ndest=purrr::map_dbl(dest_b,length)) %>% pull(ndest) %>% max() # 10 with new method 
  
  all_fluxes2 <- all_fluxes1 %>%
    mutate(dest_b_long = purrr::map(dest_b,function(x) c(x,rep(NA,(ndest-length(x))))),# is NA the best choice here?
           dest_k_long = purrr::map(dest_k,function(x) c(x,rep(NA,(ndest-length(x))))),# is NA the best choice here?
           exchange_long = purrr::map(exchange,function(x) c(x,rep(0,(ndest-length(x)))))) %>% # is 0 the best choice here?
    select(-c(dest_b,dest_k,exchange)) %>%
    ungroup()
  
  #3. Set up 3 arrays, one for each variable
  dest_b <- all_fluxes2 %>% # destination boxes
    select(dest_b_long) %>% 
    unlist() %>% 
    matrix(ncol=ndest,byrow=T) %>% 
    t() %>%
    as.array(dim = c(ndest,nlayer*nbox*ntime)) # turn to array to pack to NetCDF
  
  dest_k <- all_fluxes2 %>% # destination layers
    select(dest_k_long) %>% 
    unlist() %>% 
    matrix(ncol=ndest,byrow=T) %>% 
    t() %>%
    as.array(dim = c(ndest,nlayer*nbox*ntime)) # turn to array to pack to NetCDF
  
  exchange <- all_fluxes2 %>% #exchanges
    select(exchange_long) %>% 
    unlist() %>% 
    matrix(ncol=ndest,byrow=T) %>% 
    t() %>%
    as.array(dim = c(ndest,nlayer*nbox*ntime)) # turn to array to pack to NetCDF
  
  # the flux is still in m3/s here. Convert to m3 per time step optionally
  if(isTRUE(fluxperts)) {
    exchange <- exchange * 60*60*12
  }
  
  #####################################################################################
  # Pack to NetCDF
  
  # Here we need to write out the variables and the flows to a NetCDF file.
  # TODO: a lot of these should be args of the main function
  
  # change the dimension of the arrays for dest and exchange so that instead of being a large matrix they become a 4-dimensional array. The below seems to work
  dim(dest_b) <- c(ndest,nlayer,nbox,ntime)
  dim(dest_k) <- c(ndest,nlayer,nbox,ntime)
  dim(exchange) <- c(ndest,nlayer,nbox,ntime)
  
  # info on the run
  this_geometry <- "GOA_WGS84_V4_final.bgm"
  this_title <- "Advection between boxes"
  
  # the below is entered manually but it may be automated from the BGM file instead, especially using rbgm
  #options depth dimension
  # depth_bins <- 1:7 # 0 30 70 100 300 500 2969 in GOA. It seems to include sediment in PS
  d_units <- "depth layers"
  
  #options exchanges
  ex_units <- "m3"
  
  #options time dimension
  timestep <- 12 # 12 hour timesteps
  t_units <- paste("seconds since",time[1],sep=" ") # should it be absolute or relative to the file?
  # time.unit.length <- 2 # years
  seconds_timestep <- 60*60*12
  time_array <- array((1:ntime)*seconds_timestep)
  
  make_hydro <- function(eachvariable, nc_name, t_units, seconds_timestep, this_title, this_geometry, time_array, exchange_array, dest_b_array, dest_k_array, nbox, nlayer, ndest) {
    
    nc_file <- create.nc(nc_name)
    
    dim.def.nc(nc_file, "t", unlim=TRUE)
    dim.def.nc(nc_file, "b", nbox) # manual 
    dim.def.nc(nc_file, "z", nlayer) # manual I am not sure about this, do we need the sediment layer in the fluxes??
    dim.def.nc(nc_file, "dest", ndest) # manual 
    
    var.def.nc(nc_file, "t", "NC_DOUBLE", "t")
    var.def.nc(nc_file, eachvariable, "NC_DOUBLE", c("dest","z","b","t"))
    var.def.nc(nc_file, "dest_b", "NC_INT", c("dest","z","b","t"))
    var.def.nc(nc_file, "dest_k", "NC_INT", c("dest","z","b","t"))
    
    att.put.nc(nc_file, eachvariable, "_FillValue", "NC_DOUBLE", 0)
    att.put.nc(nc_file, "dest_b", "_FillValue", "NC_INT", -1)
    att.put.nc(nc_file, "dest_k", "_FillValue", "NC_INT", -1)
    att.put.nc(nc_file, "exchange", "units", "NC_CHAR", ex_units)
    att.put.nc(nc_file, "t", "units", "NC_CHAR", t_units)
    att.put.nc(nc_file, "t", "dt", "NC_DOUBLE", seconds_timestep)
    att.put.nc(nc_file, "NC_GLOBAL", "title", "NC_CHAR", this_title)
    att.put.nc(nc_file, "NC_GLOBAL", "geometry", "NC_CHAR", this_geometry)
    att.put.nc(nc_file, "NC_GLOBAL", "parameters", "NC_CHAR", "")
    
    var.put.nc(nc_file, "t", time_array)
    var.put.nc(nc_file, "dest_b", dest_b_array)
    var.put.nc(nc_file, "dest_k", dest_k_array)
    var.put.nc(nc_file, eachvariable, exchange_array)
    
    close.nc(nc_file)
  }
  
  make_hydro("exchange", 
             nc_name=paste0("../../outputs/complete_from_Emily/Script1/2017/monthly/forcings/hydro_hdNov2022/goa_hydro_",monthyear,".nc"), 
             t_units, 
             seconds_timestep, 
             this_title, 
             this_geometry, 
             time_array, 
             exchange_array=exchange, 
             dest_b_array=dest_b, 
             dest_k_array=dest_k, 
             nbox=nbox, 
             nlayer=nlayer, 
             ndest=ndest)
  
}

lapply(hydro_files, make_exchange_nc, fluxperts=T, hyperdiff=2)

# NOTE: for some reason the line

# :parameters = "" ;

# does not stick in the global attributes. Also, after this script we get monthly nc files. They can be concatenated with:
# cdo mergetime *.nc goa_hydro_2017.nc
# However, doing so causes the loss of line

# t:dt = 43200. ;

# For now the only way is to turn the merged file to cdf, add the missing lines, and turn back to .nc
# Using cat instead of mergetime jumbles the time steps due to how we are naming the files/folders (TODO: fix this, although mergetime is safer in general)

