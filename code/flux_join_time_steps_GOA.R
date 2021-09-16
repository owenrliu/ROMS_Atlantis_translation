# Translation of GOA ROMS to input for HydroConstruct

start_time <- Sys.time()

# Packages and settings
library(raster)
library(tidyverse)
library(sf)
library(here)
library(tidync)
library(rbgm)
library(viridis)
library(angstroms)
library(tabularaster)

select <- dplyr::select
here <- here::here
map <- purrr::map
options(dplyr.summarise.inform=FALSE)

# read ROMS data
romsfile <- here::here('data','roms','nep5_avg_0804.nc')
roms <- tidync(romsfile)
# read GOA ROMS grid
romsfile2 <- here::here('data','roms','NEP_grid_5a.nc')
roms2 <- tidync(romsfile2)
# read Atlantis BGM
atlantis_bgm <- read_bgm(here('data','atlantis','GOA_WGS84_V4_final.bgm'))
#Atlantis geometry as an sf shapefile
atlantis_sf <- atlantis_bgm %>% box_sf()

# get list of ROMS variables
roms_vars <- hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })
# ... and grid variables
roms2_vars <- hyper_grids(roms2) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables asssociated with that grid and make a reference table
  purrr::map_df(function(x){
    roms2 %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

# extract time
time_grd <- roms_vars %>% filter(name=="ocean_time") %>% pluck('grd')
roms_time <- roms %>% activate(time_grd) %>% hyper_tibble() %>% pull()

# ROMS grids
# Rho grid
# Horizontal coordinates: xi, eta, lat and lon of rho points
latlon_rhogrd <- roms2_vars %>% filter(name=="lat_rho") %>% pluck('grd')
roms_rho <- roms2 %>% activate(latlon_rhogrd) %>% hyper_tibble() %>%
  select(lon_rho,lat_rho,xi_rho,eta_rho) %>% 
  mutate(rhoidx=row_number()) # add index

# Vertical coordinates: h, s_rho, Cs_r, Cs_w
h <- roms2 %>% activate(latlon_rhogrd) %>% hyper_tibble() %>% select(xi_rho,eta_rho,h)
s_rho_grd <- roms_vars %>% filter(name=="s_rho") %>% pluck('grd')
s_rho <- roms %>% activate(s_rho_grd) %>% hyper_tibble()
## Cs_r is the S-coord stretching (length is num layers from -1 to 0 describing what portion of the w.c. each layer spans)
## We pull the Cs_r values from the roms ncdf
## one Cs_r value per s-coordinate
Cs_r <- s_rho %>% pluck('Cs_r')
s_w_grd <- roms_vars %>% filter(name=="s_w") %>% pluck('grd')
s_w <- roms %>% activate(s_w_grd) %>% hyper_tibble()
Cs_w <- s_w %>% pluck('Cs_w')

# u and v grids
# u and v grids from ROMS data
latlon_ugrd <-roms2_vars %>% filter(name=="lat_u") %>% pluck('grd')
latlon_vgrd <-roms2_vars %>% filter(name=="lat_v") %>% pluck('grd')

# pull the lon/lats
roms_u <- roms2 %>% activate(latlon_ugrd) %>% hyper_tibble() %>% select(lon_u,lat_u,xi_u,eta_u) %>% mutate(uidx=row_number())
roms_v <- roms2 %>% activate(latlon_vgrd) %>% hyper_tibble() %>% select(lon_v,lat_v,xi_v,eta_v) %>% mutate(vidx=row_number())

# append coordinates to the native lat lon from ROMS
append_xy_coords <- function(lonlatdat,xyproj=atlantis_bgm$extra$projection,lon_col="lon_rho",lat_col="lat_rho"){
  lonlatdat %>% 
    st_as_sf(coords=c(lon_col,lat_col),crs=4326,remove=F) %>%  # convert to spatial object
    st_transform(xyproj) %>%  # convert to Atlantis coords
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2]) # grab x and y coordinates and add them as attributes
}

rhoxy<- append_xy_coords(roms_rho,lon_col="lon_rho",lat_col="lat_rho") %>% mutate(rhoidx=row_number())
uxy <- append_xy_coords(roms_u,lon_col="lon_u",lat_col="lat_u")%>% mutate(uidx=row_number())
vxy <- append_xy_coords(roms_v,lon_col="lon_v",lat_col="lat_v")%>% mutate(vidx=row_number())

# join ROMS grids with Atlantis boxes (rho grid) and faces (u and v grids)
# Boxes
boxes_rho_join <- atlantis_sf %>% st_join(rhoxy)

# Faces
faces <- atlantis_bgm$faces %>% select(-label)

faces_sf <- atlantis_bgm %>% face_sf() %>% 
  mutate(label = 0:(length(label)-1)) %>% # creates a new index 'face_id' starting from 0 and increasing, as the 'label' column produced by rbgm::face_sf() is incorrect (tested in R 4.0.4)
  # join attribute data
  left_join(faces,by=c('label'='.fx0')) %>% 
  rename(.fx0=label)

# construct a buffer around each face
faces_buffer <- st_buffer(faces_sf,dist=10000) #TODO: make this an argument
# join u points
faces_u_join <- faces_buffer %>% st_join(uxy)
# join v points
faces_v_join <- faces_buffer %>% st_join(vxy)

# just a couple of checks
# insert a check for which boxes have no overlapping rho points
empty_boxes<- boxes_rho_join %>% 
  st_set_geometry(NULL) %>% 
  filter(is.na(rhoidx)) %>% 
  select(.bx0) %>% 
  distinct() %>% pull()
paste0("Atlantis boxes with no ROMS points are boxes ",paste(empty_boxes,collapse = ","))

# ... and one for which faces do not intercept u and v points
empty_faces<- faces_u_join %>% 
  st_set_geometry(NULL) %>% 
  filter(is.na(uidx)) %>% 
  select(.fx0) %>% 
  distinct() %>% pull()
paste0("Atlantis faces with no u points are faces ",paste(empty_faces,collapse = ","))

# get indeces of rho, u, and v points that overlap with Atlantis geometry, to subset large ROMS files and reduce memory chokes
min_xi_rho <- min(boxes_rho_join$xi_rho, na.rm = TRUE)
max_xi_rho <- max(boxes_rho_join$xi_rho, na.rm = TRUE)
min_eta_rho <- min(boxes_rho_join$eta_rho, na.rm = TRUE)
max_eta_rho <- max(boxes_rho_join$eta_rho, na.rm = TRUE)

min_xi_u <- min(faces_u_join$xi_u, na.rm = TRUE)
max_xi_u <- max(faces_u_join$xi_u, na.rm = TRUE)
min_eta_u <- min(faces_u_join$eta_u, na.rm = TRUE)
max_eta_u <- max(faces_u_join$eta_u, na.rm = TRUE)

min_xi_v <- min(faces_v_join$xi_v, na.rm = TRUE)
max_xi_v <- max(faces_v_join$xi_v, na.rm = TRUE)
min_eta_v <- min(faces_v_join$eta_v, na.rm = TRUE)
max_eta_v <- max(faces_v_join$eta_v, na.rm = TRUE)

# Atlantis depth
# enter the depth breaks in the model
atlantis_z <- c(-30,-100,-200,-500,-1000,-4000) #TODO: make this an argument

# small function to build enough depth layers for each box, starting with the defined layers above
# in this function, botz is the bottom depth given for each box, which is available in the .bgm file 
build_Atlantis_depths <- function(botz,lyrs){
  # bottom of each layer, starting from shallowest
  lyrbot<-lyrs
  # layers to use are all those that are shallower than the given botz
  lyr_vec <- lyrbot[lyrbot>botz]
  # the depth of the deepest layer is equal to botz
  lyr_vec <- c(lyr_vec,botz)
  # in Atlantis, each box has the same number of depth layers, but some layers have zero thickness
  # so we have to pad with zeroes to make all boxes have the same number of layers
  nzeroes <- length(lyrs)-length(lyr_vec)
  lyr_vec <- c(lyr_vec,rep(0,nzeroes))
  return(lyr_vec)
}
# construct the depth profiles of each Atlantis box

atlantis_depths <- atlantis_bgm$boxes %>% select(.bx0,botz) %>% 
  # apply the function above to create the layers
  mutate(maxz=purrr::map(botz,~build_Atlantis_depths(.,lyrs=atlantis_z))) %>% 
  unnest(cols=c(maxz)) %>% 
  # add a minimum depth for each layer
  group_by(.bx0) %>% 
  mutate(minz=lag(maxz,1,default = 0),atlantis_layer=1:length(atlantis_z)) %>% 
  # add a layer thickness calculation
  mutate(dz=minz-maxz) %>% 
  # "dummy" layers (layers too deep for a given Atlantis box) should have minz and dz=0
  mutate(minz=ifelse(maxz==0,0,minz),dz=ifelse(maxz==0,0,dz)) %>% 
  ungroup() %>% 
  select(.bx0,atlantis_layer,minz,maxz,dz)

# face depths (for later flux calcs)
face_depths <- faces %>% 
  left_join(atlantis_depths %>% select(.bx0,atlantis_layer,dz),by=c("left"=".bx0")) %>% 
  mutate(left_area=length*dz) %>% 
  rename(dz_left=dz) %>% 
  left_join(atlantis_depths %>% select(.bx0,atlantis_layer,dz),by=c("right"=".bx0","atlantis_layer")) %>% 
  mutate(right_area=length*dz) %>% 
  rename(dz_right=dz) %>% 
  # area of the face is the smallest area, maintaining NAs (if one box is deeper than its neighbor, no flux)
  rowwise() %>%
  mutate(dz_max = max(dz_left,dz_right),
         face_area=ifelse((left_area>0 & right_area >0), min(left_area, right_area), NA)) 

# For GOA: custom version of romshcoords()
# TODO: put this in a separate file when packaging

ncget <- function(x, varname) {
  nc <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_get(nc, varname)
}

set_indextent <- function(x) {
  setExtent(x, extent(0, ncol(x), 0, nrow(x)))
}

romshcoords_goa <- function(x, y, grid_type = "rho", slice, ..., S = "Cs_r", depth = "h", simple = FALSE){
  h <- romsdata(x, varname = depth)
  Cs_r <- ncget(y, S)
  v <- values(h)
  if (simple) {
    ## simplistic, early version - probably should be defunct
    out <- set_indextent(brick(array(rep(rev(Cs_r), each = length(v)) * v, 
                                     c(ncol(h), nrow(h), length(Cs_r))), transpose = TRUE))
  } else {
    grid_type <- match.arg(tolower(grid_type),c("rho","psi","u","v","w"))
    
    Vtransform <- as.integer(ncget(y,"Vtransform"))
    if (!Vtransform %in% c(1,2)) stop("Vtransform must be 1 or 2")
    
    hc <- ncget(y,"hc")
    
    depth_grid <- if (grid_type=="w") "w" else "rho"
    
    zeta <- if (missing(slice)) 0 else stop("slice not supported yet")##angstroms::romsdata2d(x,"zeta",slice=slice,transpose=FALSE)
    N <- length(ncget(y,"Cs_r"))
    Np <- N+1
    
    h <- ncget(x,"h")
    hmin <- min(h)
    hmax <- max(h)
    
    Lp <- dim(h)[1]
    Mp <- dim(h)[2]
    L <- Lp-1
    M <- Mp-1
    
    z <- array(NA,dim=c(Lp,Mp,if (grid_type=="w") Np else N))
    
    ## Compute vertical stretching function, C(k):
    ##stretch <- stretching(x,depth_grid)
    if (depth_grid=="w") {
      stretch <- list(C=ncget(y,"Cs_w"),s=ncget(y,"s_w"))
    } else {
      stretch <- list(C=ncget(y,"Cs_r"),s=ncget(y,"s_rho"))
    }
    
    ## Average bathymetry and free-surface at requested C-grid type.
    if (grid_type=="rho") {
      hr <- h
      zetar <- zeta
    } else if (grid_type=="psi") {
      hp <- 0.25*(h[1:L,1:M]+h[2:Lp,1:M]+h[1:L,2:Mp]+h[2:Lp,2:Mp])
      zetap <- 0.25*(zeta[1:L,1:M]+zeta[2:Lp,1:M]+zeta[1:L,2:Mp]+zeta[2:Lp,2:Mp])
    } else if (grid_type=="u") {
      hu <- 0.5*(h[1:L,1:Mp]+h[2:Lp,1:Mp])
      zetau <- 0.5*(zeta[1:L,1:Mp]+zeta[2:Lp,1:Mp])
    } else if (grid_type=="v") {
      hv <- 0.5*(h[1:Lp,1:M]+h[1:Lp,2:Mp])
      zetav <- 0.5*(zeta[1:Lp,1:M]+zeta[1:Lp,2:Mp])
    } else if (grid_type=="w") {
      hr <- h
      zetar <- zeta
    } else {
      stop("unsupported grid_type: ",grid_type)
    }
    
    ## Compute depths (m) at requested C-grid location.
    
    if (Vtransform == 1) {
      if (grid_type=="rho") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hr
          z[,,k] <- z0 + zetar*(1.0 + z0/hr)
        }
      } else if (grid_type=="psi") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hp
          z[,,k] <- z0 + zetap*(1.0 + z0/hp)
        }
      } else if (grid_type=="u") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hu
          z[,,k] <- z0 + zetau*(1.0 + z0/hu)
        }
      } else if (grid_type=="v") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hv;
          z[,,k] <- z0 + zetav*(1.0 + z0/hv)
        }
      } else if (grid_type=="w") {
        z[,,1] <- -hr
        for (k in seq(from=2,to=Np,by=1)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hr
          z[,,k] <- z0 + zetar*(1.0 + z0/hr)
        }
      } else {
        stop("unsupported grid_type: ",grid_type)
      }
    } else if (Vtransform == 2) {
      if (grid_type=="rho") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hr)/(hc+hr)
          z[,,k] <- zetar+(zeta+hr)*z0
        }
      } else if (grid_type=="psi") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hp)/(hc+hp)
          z[,,k] <- zetap+(zetap+hp)*z0
        }
      } else if (grid_type=="u") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hu)/(hc+hu)
          z[,,k] <- zetau+(zetau+hu)*z0
        }
      } else if (grid_type=="v") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hv)/(hc+hv)
          z[,,k] <- zetav+(zetav+hv)*z0
        }
      } else if (grid_type=="w") {
        for (k in seq_len(Np)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hr)/(hc+hr)
          z[,,k] <- zetar+(zetar+hr)*z0
        }
      } else {
        stop("unsupported grid_type: ",grid_type)
      }
    } else {
      stop("Vtransform must be 1 or 2")
    }
    ## FIXME all these flips and twirls can be applied more efficiently (or avoided)
    ## though should layers start at the surface and go down or ...
    
    out <- raster::flip(set_indextent(raster::brick(z, transpose = TRUE)), "y")
    ## NO - we want to start at the bottom, so we match romsdata3d
    #out <- raster::subset(out, rev(seq_len(raster::nlayers(out))))
    
  } 
  
  out
}

# convert ROMS s-coordinates to depth with Mike Sumner's angstroms package

romsdepths <- romshcoords_goa(x = romsfile2, y = romsfile, S = "Cs_r", depth = "h")

# using tabularaster to convert to tibble
# and a indexing template with "by_column" filling
romsi <- crossing(xi_rho=1:dim(romsdepths)[2],eta_rho=1:dim(romsdepths)[1]) %>% arrange(-eta_rho) %>% mutate(cellindex=row_number()) # making sure that the join by cellindex below is correct - doing this for consistency with the way tabularaster::as_tibble() unpacks the raster cells 
romsdepthsdf <- tabularaster::as_tibble(romsdepths,dim=F) %>% 
  arrange(cellindex) %>% 
  left_join(romsi,by='cellindex') %>% 
  set_names(c("romsdepth","cellindex","xi_rho","eta_rho")) %>% 
  group_by(cellindex,xi_rho,eta_rho) %>% 
  nest(romsdepth=c(romsdepth)) %>% ungroup() %>% 
  mutate(romsdepth=purrr::map(romsdepth,function(x)x[['romsdepth']])) %>%
  filter(between(xi_rho, min_xi_rho, max_xi_rho) & between(eta_rho, min_eta_rho, max_eta_rho))

# do the same for depth at w
romsdepths_w <- romshcoords_goa(x = romsfile2, y = romsfile, grid_type = "w", S = "Cs_w", depth = "h")
#romsdepths_w <- subset(romsdepths_w, dim(romsdepths_w)[3]:1) # angstroms returns depths from shallowest to deepest, while below we extract ROMS variables from deepest to shallowest. flipping here or else it maps ROMS variables upside-down

# using tabularaster to convert to tibble
romsdepthsdf_w <- tabularaster::as_tibble(romsdepths_w,dim=F) %>% 
  arrange(cellindex) %>% 
  left_join(romsi,by='cellindex') %>% 
  set_names(c("romsdepth","cellindex","xi_rho","eta_rho")) %>% 
  group_by(cellindex,xi_rho,eta_rho) %>% 
  nest(romsdepth=c(romsdepth)) %>% ungroup() %>% 
  mutate(romsdepth=purrr::map(romsdepth,function(x)x[['romsdepth']])) %>%
  filter(between(xi_rho, min_xi_rho, max_xi_rho) & between(eta_rho, min_eta_rho, max_eta_rho))

# Matching keys for boxes and faces
# TODO: shift this above - this will not need repeating at each time step
# build a matching key for Atlantis boxes...
boxes_rho_thin <- boxes_rho_join %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(.bx0,xi_rho,eta_rho,rhoidx,area) %>% 
  drop_na()
boxes_rho_join_with_depth <- atlantis_depths %>% 
  left_join(boxes_rho_thin,by=c(".bx0")) %>% 
  ungroup() %>% 
  drop_na()

# make vectors of rho point indices
boxes_depths_rhopts <- boxes_rho_join_with_depth %>% 
  select(-dz) %>% 
  group_by(.bx0,atlantis_layer,minz,maxz) %>% 
  nest(rhopts=c(xi_rho,eta_rho,rhoidx)) %>% 
  mutate(rhovec=purrr::map(rhopts,~pluck(.,'rhoidx')))

# ... and for Atlantis faces
faces_u_thin <- faces_u_join %>%
  st_set_geometry(NULL) %>%
  select(.fx0,xi_u,eta_u,uidx) %>%
  drop_na()
faces_u_join_with_depth <- face_depths %>%
  left_join(faces_u_thin, by = c(".fx0")) %>%
  ungroup() %>%
  drop_na() %>% # turn this on or off depending on whether we want to have NA fluxes in non-existing layers or not, cannot recall what HydroCOnstruct wants
  select(.fx0,atlantis_layer,dz_max,xi_u,eta_u,uidx,face_area) %>%
  group_by(uidx,.fx0) %>%
  mutate(maxz=-cumsum(dz_max),minz=-lag(-maxz,default=0)) %>%
  ungroup()%>%
  select(-dz_max)

# again for v
faces_v_thin <- faces_v_join %>%
  st_set_geometry(NULL) %>%
  select(.fx0,xi_v,eta_v,vidx) %>%
  drop_na()
faces_v_join_with_depth <- face_depths %>%
  left_join(faces_v_thin, by = c(".fx0")) %>%
  ungroup() %>%
  drop_na() %>%
  select(.fx0,atlantis_layer,dz_max,xi_v,eta_v,vidx,face_area) %>%
  group_by(vidx,.fx0) %>%
  mutate(maxz=-cumsum(dz_max),minz=-lag(-maxz,default=0)) %>%
  ungroup()%>%
  select(-dz_max)

# make vectors of u point indices
faces_depths_upts <- faces_u_join_with_depth %>% 
  group_by(.fx0,atlantis_layer,minz,maxz) %>% 
  nest(upts=c(xi_u,eta_u,uidx)) %>% 
  mutate(uvec=map(upts,~pluck(.,'uidx'))) %>%
  ungroup()
# and v
faces_depths_vpts <- faces_v_join_with_depth %>% 
  group_by(.fx0,atlantis_layer,minz,maxz) %>% 
  nest(vpts=c(xi_v,eta_v,vidx)) %>% 
  mutate(vvec=map(vpts,~pluck(.,'vidx'))) %>%
  ungroup()
# now join
faces_depths_uvpts <- faces_depths_upts %>% full_join(faces_depths_vpts, by = c('.fx0','atlantis_layer','maxz','minz','face_area'))

# pull ROMS static variables (salt, temp, w, etc.) from Vertical interpolation
# function to wrap 'spline' and apply it to ROMS data and depth vector (1 m intervals)
interp_foo <- function(romsdepths,romsvar) {
  depths_out <- seq(round(min(romsdepths)),0,by=1) # 1m interpolation, starting from deepest
  interp <- spline(romsdepths,romsvar,xout=depths_out) %>% pluck('y')
  return(tibble(depth=depths_out,val=interp))
}

# function that pulls variables from ROMS at rho points, and interpolates the vertical values at each 1 m
interpolate_var <- function(variable, time_step){
  grd <- roms_vars %>% filter(name==variable) %>% pluck('grd')
  # pull the env data
  # interpolate the env data
  # do this step conditional to join with the appropriate depth data frame depending on the variable
  # if variable is horizontal velocity
  if(variable == "u") {
    dat <- roms %>% activate(grd) %>%
      hyper_tibble(select_var=variable, 
                   xi_u = between(xi_u, min_xi_u, max_xi_u), 
                   eta_u = between(eta_u, min_eta_u, max_eta_u),
                   ocean_time = ocean_time == time_step)
    
    interp_dat <- dat %>% 
      dplyr::select(xi_u,eta_u,!!variable) %>% 
      nest(data=c(!!variable))%>% 
      mutate(evar=purrr::map(data,~.[[1]])) %>% 
      inner_join(romsdepthsdf,by=c('xi_u'='xi_rho','eta_u'='eta_rho')) %>% 
      inner_join(roms_u,by=c('xi_u','eta_u')) 
  } else if (variable == "v") {
    dat <- roms %>% activate(grd) %>%
      hyper_tibble(select_var=variable, 
                   xi_v = between(xi_v, min_xi_v, max_xi_v), 
                   eta_v = between(eta_v, min_eta_v, max_eta_v),
                   ocean_time = ocean_time == time_step)
    
    interp_dat <- dat %>% 
      dplyr::select(xi_v,eta_v,!!variable) %>% 
      nest(data=c(!!variable))%>% 
      mutate(evar=purrr::map(data,~.[[1]])) %>% 
      inner_join(romsdepthsdf,by=c('xi_v'='xi_rho','eta_v'='eta_rho')) %>% 
      inner_join(roms_v,by=c('xi_v','eta_v')) 
  } else { # if variable is a state variable or vertical velocity
    dat <- roms %>% activate(grd) %>%
      hyper_tibble(select_var=variable, 
                   xi_rho = between(xi_rho, min_xi_rho, max_xi_rho), 
                   eta_rho = between(eta_rho, min_eta_rho, max_eta_rho),
                   ocean_time = ocean_time == time_step)
    
    interp_dat <- dat %>% 
      dplyr::select(xi_rho,eta_rho,!!variable) %>% 
      nest(data=c(!!variable))%>% 
      mutate(evar=purrr::map(data,~.[[1]]))
    if (variable == 'w') { # if the variable is w we need a different vertical mapping
      interp_dat <- interp_dat %>%
        inner_join(romsdepthsdf_w,by=c('xi_rho','eta_rho')) %>% 
        inner_join(roms_rho,by=c('xi_rho','eta_rho'))
      # drop NAs - there are a lot in w - might want to check why
      # idx <- unlist(lapply(interp_dat$evar, function(x)length(which(is.na(x)))/length(x)))
      # interp_dat <- interp_dat[-which(idx==1),]
    } else { # for all other state variables
      interp_dat <- interp_dat %>%
        inner_join(romsdepthsdf,by=c('xi_rho','eta_rho')) %>% 
        inner_join(roms_rho,by=c('xi_rho','eta_rho')) 
    }
  }
  interp_dat <- interp_dat %>% 
    mutate(interp = purrr::map2(romsdepth,evar,interp_foo)) %>% 
    dplyr::select(-data,-evar,-romsdepth)
  return(interp_dat)
}

# From here down: repeat at each time step
# TODO: put this in a separate file / code block. 
# Also the functions that we call in here will need to go elsewhere and only be called once

interpolate_ts <- function(time_step) {

  # apply interpolate_var
  # TODO: operationalize this so that it works with whatever state variables people need to pull
  salt_interp <- interpolate_var('salt', time_step)
  temperature_interp <- interpolate_var('temp', time_step) 
  w_interp <- interpolate_var('w', time_step) 
  u_interp <- interpolate_var('u', time_step)
  v_interp <- interpolate_var('v', time_step)
  
  # match interpolated data to Atlantis boxes, faces, and layers
  match_interpolated_data <- function(interp_dat,idx_vec,minz,maxz){
    #idx <- names(interp_dat[grep('idx',names(interp_dat))])
    interp_var<- interp_dat
    if('uidx'%in%names(interp_var)){
      interp_var<-interp_var%>%filter(uidx%in%idx_vec)
    }else if ('vidx'%in%names(interp_var)){
      interp_var<-interp_var%>%filter(vidx%in%idx_vec)
    }else{
      interp_var<-interp_var%>%filter(rhoidx%in%idx_vec)
    }
    interp_var<-interp_var %>%    # filter depths
      mutate(interp_within_depths=purrr::map(interp, ifelse(maxz==0,function(df) df %>% mutate(val = NA), function(df) df %>% filter(between(depth,maxz,minz))))) %>% # write out NAs if the thickness of the depth layer is 0 m, i.e. for non-existing Atlantis layers
      # unpack the interpolated data and calculate a mean across relevant rho points and depths
      pluck('interp_within_depths') %>% 
      bind_rows() %>% 
      pluck('val') %>% 
      mean(na.rm=T)
  }
  
  # state variables, but not w
  boxes_statevars_interp <- boxes_depths_rhopts %>% 
    # apply the function to the interpolated salinity and temperature data
    mutate(salt_mean=purrr::pmap_dbl(list(idx_vec=rhovec,minz=minz,maxz=maxz),match_interpolated_data,interp_dat=salt_interp)) %>% 
    mutate(temperature_mean=purrr::pmap_dbl(list(idx_vec=rhovec,minz=minz,maxz=maxz),match_interpolated_data,interp_dat=temperature_interp))
  
  # now extract w at the interface between Atlantis layers
  w_at_interface <- boxes_rho_join_with_depth %>% left_join(w_interp, by = c('xi_rho','eta_rho','rhoidx'))
  
  # function to match w at rho points within each box at the exact depth where Atlantis layers meet
  pull_w_at_depth <- function(maxz,interp){
    w_at_depth <- interp %>% filter(depth==round(maxz)) %>% select(val) %>% pull()
    return(w_at_depth)
  }
  
  # get mean w at the interface between layers
  # TODO: is there a way to extract info regarding the direction of w from the ROMS file? I.e., is positive w upward or downward?
  boxes_w_interp <- w_at_interface %>% 
    drop_na() %>% # get rid of NAs
    filter(!maxz==0) %>% # drop all rows where maxz = 0, because these are "empty" layers and we do not need w at the surface
    group_by(rhoidx) %>%
    mutate(w_at_depth = purrr::map2(maxz,interp,pull_w_at_depth)) %>%
    select(-interp) %>%
    unnest(cols = c(w_at_depth)) %>%
    ungroup() %>% 
    group_by(.bx0,atlantis_layer,minz,maxz,area) %>%
    summarise(mean_w = mean(w_at_depth,na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(net_w=mean_w*area)
  
  # pull u and v
  faces_uv_interp <- faces_depths_uvpts %>% 
    # apply the function to the interpolated salinity and temperature data
    mutate(u_mean=purrr::pmap_dbl(list(idx_vec=uvec,minz=minz,maxz=maxz),match_interpolated_data,interp_dat=u_interp),
           v_mean=purrr::pmap_dbl(list(idx_vec=vvec,minz=minz,maxz=maxz),match_interpolated_data,interp_dat=v_interp))
  #mutate(temperature_mean=purrr::pmap_dbl(list(rhoidx_vec=rhovec,minz=minz,maxz=maxz),match_interpolated_data,interp_dat=temperature_interp))
  
  # calclate horizontal fluxes across faces
  # from the faces_sf spatial object, let's retrieve geometry information of the points that define each face, and calculate the angle. Also calculate its sine and cosine 
  
  get_angle <- function(geomstring) {
    coords <- st_coordinates(geomstring)
    this_atan <- atan2(coords[2,2]-coords[1,2], coords[2,1]-coords[1,1]) # as per R help, atan2(y,x)
    return(this_atan)
  }
  
  faces_sf <- faces_sf %>% rowwise() %>% mutate(face_angle = get_angle(geometry),
                                                sine_new = sin(face_angle),
                                                cosine_new = cos(face_angle)) %>%
    select(.fx0, face_angle) %>% st_set_geometry(NULL) # only keep face id and angle
  
  # now use the same code as above, except the line that works out the angles from sine and cosine
  
  uv_interp <- faces_uv_interp %>% select(.fx0,atlantis_layer,face_area,u_mean,v_mean) %>%
    left_join(faces_sf,by=c('.fx0')) %>%
    mutate(uv_angle=atan2(v_mean, u_mean),uv_mag=sqrt(u_mean^2+v_mean^2)) %>%
    mutate(diff_angle = uv_angle - face_angle, # beta - alpha
           angle_new = ifelse((diff_angle>=0 & diff_angle<=pi) | (diff_angle<0 & diff_angle>=-pi), diff_angle, # add or subtract 2*pi as needed (9)see above)
                              ifelse(diff_angle < -pi, diff_angle+2*pi, diff_angle-2*pi)),
           orthogonal_flux = uv_mag * sin(angle_new)) # get the component of the flux that is perpendicular to the face. Because of the angle transformation (+/- 2*pi), multiplying by the sine of the angle between uv and the face returns the correct sign for the flux: negative for LR, positive for RL.
  
  fluxes_interp <- uv_interp %>% 
    mutate(gross = orthogonal_flux * face_area) %>%
    select(.fx0, atlantis_layer, gross)
  
  # write results to DAT file
  #TODO: look into what HydroConstruc needs. Tab-delimited? How are the faces numbered? What happens to w in HydroConstruct?
  
  variables_out <- boxes_statevars_interp %>%
    left_join(boxes_w_interp %>% select(.bx0,atlantis_layer,net_w),by=c('.bx0','atlantis_layer')) %>%
    drop_na() %>%
    mutate(time_step = time_step, maxz = -maxz) %>% # how to turn this into iteration over time stesp? Either ts from NetCDF, or just iteration # of the function
    ungroup() %>%
    select(time_step,.bx0,maxz,net_w,temperature_mean,salt_mean) %>%
    mutate_all(formatC,format='e',digits=8) %>%
    set_names(c('Time Step','Polygon number','Depth Layer [m]','Vertical velocity [m3/s]','Average Temperature [Celsius]','Average Salinity [PartPer1000]'))

  # think about the below when we get to do this for more than one time step
  if(time_step == roms_time[1]) {
    write.table(variables_out, 'outputs/state_vars_test.dat', quote=FALSE, row.names = FALSE, sep = '\t', append = TRUE)
  } else {
    write.table(variables_out, 'outputs/state_vars_test.dat', quote=FALSE, row.names = FALSE, sep = '\t', append = TRUE, col.names = FALSE)
  }
  
  # write out fluxes
  # We need to relabel the faces. Format 9 in HC takes faces numbered as 0,1,2,...,n for each box, instead of with their unique identifier .fx0. Start from faces_sf.
  # Workflow:
  # * Join flux data set with face data set for box to the left, by face. 
  # * Join flux data set with face data set for box to the right, by face. **Flip the sign to the flux**, because negative fluxes are LR, but here the box refers to the box to the right, and the LR flux is actually entering the box, thus requiring to be positive.
  # * Do a rbind()
  # * Pad empty depth layers with 0's
  # This way fluxes will be duplicated, which I believe is what we need.
  fluxes_out_left <- fluxes_interp %>%
    left_join(faces %>% select(.fx0,left), by='.fx0') %>%
    mutate(time_step=1) %>%
    select(left,.fx0,atlantis_layer,gross) %>%
    set_names(c('Polygon_number','Face_number','Depth_layer','Flux_m3s'))
  
  fluxes_out_rigth <- fluxes_interp %>%
    left_join(faces %>% select(.fx0,right), by='.fx0') %>%
    mutate(gross = -gross) %>%
    select(right,.fx0,atlantis_layer,gross) %>%
    set_names(c('Polygon_number','Face_number','Depth_layer','Flux_m3s'))
  
  fluxes_out <- rbind(fluxes_out_left,fluxes_out_rigth) %>% 
    distinct() %>% 
    arrange(Polygon_number,Face_number,Depth_layer) 
  
  # make an index for the new faces
  face_idx <- fluxes_out %>%
    select(Polygon_number,Face_number) %>% 
    distinct() %>%
    group_by(Polygon_number) %>%
    mutate(Face_new=row_number()) %>%
    ungroup()
  
  # join the index to the fluxes_out data set
  fluxes_out <- fluxes_out %>% left_join(face_idx,by=c('Polygon_number','Face_number'))
  
  # pad missing depth layers
  fluxes_out <- fluxes_out %>% complete(Depth_layer, nesting(Polygon_number,Face_new)) %>% 
    select(Polygon_number,Face_new,Depth_layer,Flux_m3s) %>%
    arrange(Polygon_number,Face_new,Depth_layer)
  
  # add time step column (placeholder value here)
  fluxes_out <- fluxes_out %>% mutate(Time_step=1)
  
  # change NAs to 0s
  fluxes_out[is.na(fluxes_out)] <- 0
  
  # set columns in the right order, sort, and rename them as HC needs
  fluxes_out <- fluxes_out %>% 
    select(Polygon_number,Face_new,Time_step,Depth_layer,Flux_m3s) %>%
    mutate_all(formatC,format='e',digits=8) %>%
    set_names('Polygon number','Face number','Time Step (12)hr','Depth Layer','Flux [m3/s]')

  # think about the below when we get to do this for more than one time step
  if(time_step == roms_time[1]) {
    write.table(fluxes_out, 'outputs/transport_test.dat', quote=FALSE, row.names = FALSE, sep = '\t', append = TRUE)
  } else {
    write.table(fluxes_out, 'outputs/transport_test.dat', quote=FALSE, row.names = FALSE, sep = '\t', append = TRUE, col.names = FALSE)
  }
}

# apply function to time steps

# purrr::map(roms_time, possibly(interpolate_ts, NA))

end_time <- Sys.time()
end_time-start_time
