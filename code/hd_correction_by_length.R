# Calculate correction factor for fluxes for hyperdiffusion correction
# Inspired by Owen Liu's approach, we get the length of the segment perpendicular to all faces that meets the opposite, nearest face
# we divide fluxes by that. It will have to be the segment in the dest box

# read in geometry
atlantis_bgm <- read_bgm('C:/Users/Alberto Rovellini/Documents/GOA/ROMS/data/atlantis/GOA_WGS84_V4_final.bgm')
atlantis_box <- atlantis_bgm %>% box_sf()
atlantis_crs <- atlantis_bgm$extra$projection

# Procedure:
# find midpoint
# find slope of face
# find slope of perpendicular line
# given slope and midpoint, get equation of perpendicular line
# create 2 points on this line: get a range of X +- 1000 km or something like that
# cast to LINESTRING object
# see where each linestring intersects with the opposing face
# if it intersects more than once, keep only the closest
# calculate the distance

# find midpoint
faces <- atlantis_bgm$faces %>% select(-label)

faces_geom <- atlantis_bgm %>% face_sf() %>%
  mutate(.fx0 = 0:(nrow(.)-1)) %>%
  select(.fx0)

orth_line <- faces_geom %>% st_cast('POINT') %>% 
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL) %>%
  group_by(.fx0) %>%
  mutate(X2 = lead(X),
         Y2 = lead(Y)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(midpoint_X = (X+X2)/2,
         midpoint_Y = (Y+Y2)/2,
         slope_face = (Y2-Y)/(X2-X), # slope of the face
         slope_orth = -1/slope_face, # slope of perpendicular line m2 = -1/m
         intercept_orth = midpoint_Y - slope_orth*midpoint_X, # solve for q = y - m2x
         dX = 1e6 * (1/(sqrt(1+slope_orth^2))), # calculate the increment in X necessary for a 1,000 km line departing in each direction along the intersect form the midpoint
         X3 = midpoint_X+dX,
         X4 = midpoint_X-dX,
         Y3 = slope_orth * X3 + intercept_orth,
         Y4 = slope_orth * X4 + intercept_orth)
  
  
orth_line_thin <- orth_line %>% select(.fx0, X3, Y3, X4, Y4)

tmp1 <- orth_line_thin %>% select(.fx0, X3, Y3) %>% mutate(End = 'End3') %>% rename(X = X3, Y = Y3) 
tmp2 <- orth_line_thin %>% select(.fx0, X4, Y4) %>% mutate(End = 'End4') %>% rename(X = X4, Y = Y4) 

orth <- rbind(tmp1,tmp2) %>% arrange(.fx0,End) %>% select(-End)

get_linestring <- function(xcol, ycol){
  z <- st_linestring(matrix(data = c(xcol, ycol), ncol = 2))
  z <- z %>% st_sfc() %>% st_set_crs(atlantis_crs)
  return(z)
}

# make list of linestring geometries
orth_list <- orth %>%
  group_by(.fx0) %>%
  nest() %>%
  mutate(geometry = purrr::map(data, ~get_linestring(.x$X, .x$Y))) %>%
  ungroup() %>%
  select(-data)

# now for each face, find the point intersecting on the opposing face

data_list <- vis_list <- list()

for(i in 1:nrow(faces)){
  
  this_face <- faces$.fx0[i]
  
  # where is the midpoint?
  this_midpoint <- orth_line %>% 
    filter(.fx0 == this_face) %>% 
    select(midpoint_X, midpoint_Y) %>% 
    st_as_sf(coords = c('midpoint_X', 'midpoint_Y'), crs = atlantis_crs)
  
  # the orthogonal line
  this_orth <- orth_list$geometry[[this_face+1]] 
  
  # boxes to the r and l
  rbox <- faces %>% filter(.fx0 == this_face) %>% pull(right)
  lbox <- faces %>% filter(.fx0 == this_face) %>% pull(left)
  
  rbox_geom <- atlantis_box %>% filter(.bx0 == rbox) %>% select(.bx0)
  lbox_geom <- atlantis_box %>% filter(.bx0 == lbox) %>% select(.bx0)
  
  intersect_r <- st_intersection(this_orth, rbox_geom) # this contains all intersections with the box including the midpoint
  intersect_l <- st_intersection(this_orth, lbox_geom) # this contains all intersections with the box including the midpoint
  
  intersect_pts_r <- intersect_r %>% st_cast('POINT') # turn the intersection from a line to points
  intersect_pts_l <- intersect_l %>% st_cast('POINT')
  
  distances_r <- as.numeric(st_distance(intersect_pts_r, this_midpoint, by_element = T)) # calculate distances
  distances_l <- as.numeric(st_distance(intersect_pts_l, this_midpoint, by_element = T))
  
  distances_r <- distances_r[distances_r > 10] # get rid of the actual midpoint - allow for 10 m imprecision
  distances_l <- distances_l[distances_l > 10] 
  
  distances_r <- min(distances_r) # keep the closest intersection only
  distances_l <- min(distances_l) 
  
  # as df?
  data_list[[i]] <- data.frame(".fx0" = this_face, rbox, lbox, distances_r, distances_l)
  
  # for visualisation purposes
  # this one is a little mind-bending
  pt_r <- intersect_pts_r[-st_nearest_feature(this_midpoint, intersect_pts_r)] # remove midpoint
  pt_r <- pt_r[st_nearest_feature(this_midpoint, pt_r)] # now get next closest 
  pt_l <- intersect_pts_l[-st_nearest_feature(this_midpoint, intersect_pts_l)] # remove midpoint
  pt_l <- pt_l[st_nearest_feature(this_midpoint, pt_l)] # now get next closest 
  # make 3 sets of points
  coord_midpoint <- st_coordinates(this_midpoint)
  coord_r <- st_coordinates(pt_r)
  coord_l <- st_coordinates(pt_l)
  # join mid-point to r and l closest points
  segment_r <- st_linestring(rbind(coord_midpoint, coord_r)) %>% st_sfc() %>% st_set_crs(atlantis_crs)
  segment_l <- st_linestring(rbind(coord_midpoint, coord_l)) %>% st_sfc() %>% st_set_crs(atlantis_crs)
  
  # bind into one
  vis_list[[i]] <- rbind(data.frame('.fx0'=this_face, to_box = 'r', 'segment' = segment_r),
                    data.frame('.fx0'=this_face, to_box = 'l', 'segment' = segment_l)) 
  
}

all_data <- do.call(rbind, data_list)

write.csv(all_data, 'hd_by_length.csv', row.names = F)

# all_vis <- do.call(rbind, vis_list) %>% 
#   st_as_sf() %>%
#   full_join(faces %>% select(.fx0, right, left), by = '.fx0')
# 
# # # view
# # ggplot()+
# #   geom_sf(data = (atlantis_box %>% filter(box_id == 50)), aes(fill= NULL))+
# #   geom_sf(data = (all_vis %>% filter(right == 50 | left == 50)), aes(color = to_box))
