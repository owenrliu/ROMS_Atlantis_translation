library(ncdf4)
library(pracma)
library(maptools)
if (!rgeosStatus()) gpclibPermit()
library(methods)
library(geosphere)

#using corner information
#you need a file containing information on the faces of your polygons, where the coordinates are given in lon/lat (the shp file could be used as well). 
corners <- read.table("/data/felles/NoBaAtlantis/Common_files/corners_neighbours_nordic.txt", header=T)
section = 1:nrow(corners)

h_get = c(25,100,200,300, 425, 750, 1200)*-1  # Depth ranges (m) to interpolate currents to
dz = 20 		# Depth interval (m) for vertical decimation, NB! define this as large as possible, it messes up the computational costs
Nt_get = c(1)		# Time ranges to process from NetCDF file - ant tidssteg
nlay=length(h_get) # nr of layers in the model
max_neigh=20 #max nr of neighbours for a polygon
max_Ni=6 # max nr of timesteps in you roms nc file
daymean=5 #NB!! Check this in your ROMS files, this is for the NoBa model - GISS_AOM
# calculate on subset instead of whole area - reduces the computational time 
#NB - needs to be updated to the specific model
x1=1
x2=390
xt=x2-x1 #
       #
y1=1
y2=380
yt=y2-y1 

year_all=c(2009:2015) # the years you need info from

#This is information that should be in your roms file - but check the lon/lat if you get issues
mask.file=("/data/NoBaAtlantis/Common_files/AA_10km_grid.nc") #NB! 
m.f=nc_open(mask.file)
h = ncvar_get(m.f, "h") # bathymetry at rho points dim( 507 329 )
lon_rho = ncvar_get(m.f, "lon_rho") # bathymetry at rho points dim( 507 329 )
lat_rho = ncvar_get(m.f, "lat_rho") # bathymetry at rho points dim( 507 329 )
angle = ncvar_get(m.f, "angle") # bathymetry at rho points dim( 507 329 )
nc_close(m.f)

x_grab=c(380,770)
y_grab=c(570,950)


#information on C_sr, hc, s_cr etc
load('/data/NoBaAtlantis/Common_files/missing_roms_vars.RData') 

for (y in year_all){
  print(y)
  year=y
  print(year)

# arrays for creating daily values for fluxes

  cnt=1
  f.name = paste("/data/NoBaAtlantis/Forcing/NorESM_2006_2070/NorESM_year_",year,'.nc',sep="")
  ex.cb = nc_open(f.name)   

  ocean_time = ncvar_get( ex.cb, "ocean_time" )  # 
  u = ncvar_get( ex.cb, "u" ) # u-momentum component   
  v = ncvar_get( ex.cb, "v" ) # v-momentum component 
  zeta = ncvar_get( ex.cb, "zeta") # 
  s_rho=ncvar_get( ex.cb, "s_rho" )

  nc_close(ex.cb)

  Nx = dim(lon_rho)[1]
  Ny = dim(lon_rho)[2]
  Nz=dim(u)[3]
  Nt=length(ocean_time)
#
  t_exc=array(0,dim=c(length(section),length(h_get),length(ocean_time)*max_Ni))
  t_db=array(NA,dim=c(length(section),length(h_get),length(ocean_time)*max_Ni))
  t_dk=array(NA,dim=c(length(section),length(h_get),length(ocean_time)*max_Ni))
  t_boxn=array(NA,dim=c(length(section),length(h_get),length(ocean_time)*max_Ni))
  w=seq(0,1,by=(1/(ceil(365/Nt))))

       # Get the horizontal, vertical, and time dimensions of the data

      # if the time record chosen is not in the file, then exit:
      if (max(Nt_get) > Nt) {
      	print('The time record exceeds the available number ....')	
       } else {
       	print('Ready to Go!')
      }

      # Calculate the number of depth intervals to average u, v, T and S over:
      Ni=length(h_get)
      # Calculate the number of time intervals to average u, v, T and S over:
      Ti=Nt # Several time steps in each file
      

      # define arrays at internal rho-points to average velocities to and get the
      # variables at these point locations
      x_sub=array(1,dim=c(xt-2, yt-2,Ni,Ti))      # 
      y_sub=array(1,dim=c(xt-2, yt-2,Ni,Ti))      # 
      u_sub=array(1,dim=c(xt-2, yt-2,Ni,Ti))      # 
      v_sub=array(1,dim=c(xt-2, yt-2,Ni,Ti))      # 

      u_tmp=array(1,dim=c(xt-2, yt - 2, Nz,Nt))  #
      v_tmp=array(1,dim=c(xt-2,yt-2,Nz,Nt))

      x_sub=lon_rho[(x_grab[1]+2):(x_grab[2]-2),(y_grab[1]+2):(y_grab[2]-2)]
      y_sub=lat_rho[(x_grab[1]+2):(x_grab[2]-2),(y_grab[1]+2):(y_grab[2]-2)]

      h_sub=drop(h[(x_grab[1]+2):(x_grab[2]-2),(y_grab[1]+2):(y_grab[2]-2)]) #
      zang_sub=drop(angle[(x_grab[1]+2):(x_grab[2]-2),(y_grab[1]+2):(y_grab[2]-2)])
                                           
      uu= u #
      vv= v #


      uu[which(is.na(uu))]=0.0
      vv[which(is.na(uu))]=0.0

      if(Ti==1){
         zeta_sub= drop(zeta[(x1+1):(x2-1),(y1+1):(y2-2)])     
      }else{
         zeta_sub= drop(zeta[(x1+1):(x2-1),(y1+1):(y2-2),])     
      }


      zetazeta<- zeta #drop(zeta[,,Nt_get])

      #Exclude land velocity values and set the .. to zero if necessary:

      length(which(uu<1000000))/length(uu) # 

      II= which(abs(uu)>1000000) # 
      if (length(II)!=0) uu[II]=0

      II= which(abs(vv)>1000000)
      if (length(II)!=0) vv[II] =0

      # Average velocities to internal rho-point locations:

      #for (i in 1:(xt-2)) {
      #    if(Ti==1){
      #       u_tmp[i,,,]=0.5*(drop(uu[i,(y1+2):(y2-1),])+drop(uu[(i+1),(y1+2):(y2-1),]))
      #    }else{
      #       u_tmp[i,,,]=0.5*(drop(uu[i,(y1+2):(y2-1),,])+drop(uu[(i+1),(y1+2):(y2-1),,]))
      #    }
      #}

      #for (j in 1:(yt-2)) {
      #    if(Ti==1){
      #       v_tmp[,j,,]=0.5*(drop(vv[(x1+2):(x2-1),j,])+drop(vv[(x1+2):(x2-1),j+1,]))
      #    }else{
      #       v_tmp[,j,,]=0.5*(drop(vv[(x1+2):(x2-1),j,,])+drop(vv[(x1+2):(x2-1),j+1,,]))
      #    }
      #}

      u_tmp=uu[(x1+1):(x2-2),(y1+1):(y2-2),,]
      v_tmp=vv[(x1+1):(x2-2),(y1+1):(y2-2),,]


#remove NAs
     #zeta_sub[is.na(zeta_sub)]=0
     u_tmp[is.na(u_tmp)]=0
     v_tmp[is.na(v_tmp)]=0


      if (Ti==1) TiLoop = 1
      if (Ti > 1 ) TiLoop = Ti 

      for (i in 1:(xt-2)) {  
         for (j in 1:(yt-2)) {  #
            if (!is.na(zeta_sub[i,j,1])) { 
	       for (tm in 1:TiLoop){  
	          tm1 = tm
	          tm2 = tm+1
	          z0=(s_rho-Cs_r)*(hc)+ Cs_r*h_sub[i,j] 
                  if (Ti > 1) {
		     zr=z0+mean(zeta_sub[i,j,tm])*(1.0 + z0/h_sub[i,j])  
                  }else{
		     zr=z0+mean(zeta_sub[i,j])*(1.0 + z0/h_sub[i,j])     
	          }

	           # Top depth interval:
	          if (Ni >1 )	{ # if you have > 1 depth layer
		     tailzr = length(zr)
		     if (tail(zr,1) < h_get[2])  zr[tailzr] <- tail(zr,1) + 1.5*dz 
	          }	

	          # Decimate the vertical u, v, T, S values

	          zrv=seq(zr[1],tail(zr,1),by=dz)

                  ui=interp1(as.vector(zr),drop(u_tmp[i,j,,tm]),zrv,method='linear')  
	          vi=interp1(as.vector(zr),drop(v_tmp[i,j,,tm]),zrv,method='linear')

	          uE=ui*cos(zang_sub[i,j])-vi*sin(zang_sub[i,j])
	          vN=ui*sin(zang_sub[i,j])+vi*cos(zang_sub[i,j])
                  #uE=ui
                  #vN=vi

	          if (Ni > 1) { 
	              for (km in 1:(Ni-1)) { # Ni = length(h_get), i.e. depth
		         II = which((zrv >= h_get[km+1]) & (zrv < h_get[km]))
	                 if ( length(II) != 0 ) {
		            u_sub[i,j,km,tm] = mean(uE[II])
			    v_sub[i,j,km,tm] = mean(vN[II])
		         }
	              }
                  }

	           II=which(zrv <= h_get[Ni])
	           if (length(II) != 0) {
		      u_sub[i,j,Ni,tm]=mean(uE[II])
		      v_sub[i,j,Ni,tm]=mean(vN[II])
	           }
                } 
	     }
	 } # end j
      } # end i


   ############# remove NA-verdier u_sub and v_sub

   u_sub2=u_sub
   u_sub2[which(u_sub2==1.00000,arr.ind=TRUE)]=NA

   v_sub2=v_sub
   v_sub2[which(v_sub2==1.00000,arr.ind=TRUE)]=NA

   vel_section=array(NA,dim=c(max(section), length(h_get),Ti))
   exchange=array(0,dim=c(max(section), length(h_get),Ti))
   angle_sec = array(NA,dim=c(max(section), length(h_get),Ti))
   angle_vel = array(NA,dim=c(max(section), length(h_get),Ti))
   dest_k = array(NA,dim=c(max(section),length(h_get),Ti))
   dest_b = array(NA,dim=c(max(section),length(h_get),Ti))
   boxnr = array(NA,dim=c(max(section),length(h_get),Ti))
   hyperdiff=array(NA,dim=c(max(section)))

   for (tm in 1:TiLoop) { # loop by time
   for (km in 1:Ni) { # loop by depth
     qq = u_sub2[,,km,tm] # get the depth and time levels looping through now
     if (length(which(is.finite(qq)))>0) {
        II=is.finite(qq) # true false
        x_fin=x_sub[II] #  lon's
        y_fin=y_sub[II] # lats
        qq=drop(u_sub2[,,km,tm]);u_fin=qq[II]
        qq=drop(v_sub2[,,km,tm]);v_fin=qq[II]

        for(i in 1:max(section)){
           pol.factor.x = 1.01  #factor to extend line in lon/lat directions 
           pol.factor.y = 1.01
           if(i==3|i==32|i==48|i==107|i==115|i==124|i==150|i==240|i==284|i==293|i==307|i==332|i==458|i==475|i==559){ #NB!!!!! this is ONLY for the NoBa model
              pol.factor.x = 1.07 #need larger factors, otherwise too slim/short to get any points  
              pol.factor.y = 1.02 
           }
           p_f = c(pol.factor.x, pol.factor.y)
           x = c(corners$Lon1[i], corners$Lon2[i])
           y = c(corners$Lat1[i], corners$Lat2[i])
           xy=cbind(x,y)
           xy.pol = rbind(xy, rbind(xy[2,]*p_f,xy[1,]*p_f), xy[1,])    #first and last coordinate have to be the same
           pol = SpatialPolygons(list(Polygons(list(Polygon(xy.pol)), ID="1")))
           coordssat=SpatialPoints(coords=cbind(x_fin,y_fin))
	   pts.pol<-over(coordssat,pol)
           mean.u = mean(u_fin[pts.pol==1], na.rm=T)
           mean.v = mean(v_fin[pts.pol==1], na.rm=T)
           np = sum(pts.pol==1, na.rm=T)

           angle_sec[i,km,tm]=(180*atan2(x[2]-x[1],y[2]-y[1]))/pi
           angle_vel[i,km,tm]=(180*atan2(mean.u,mean.v))/pi
           vel_section[i,km,tm]=sqrt(mean.u^2+mean.v^2)
           if(!is.na(angle_vel[i,km,tm])&angle_vel[i,km,tm]>angle_sec[i,km,tm]) vel_section[i,km,tm]=(-1)*vel_section[i,km,tm]

           #compute flux over sections
           #exchange[i,km,tm]=vel_section[i,km,tm]*(distMeeus(xy[1,],xy[2,])*corners$MeanDepth[i])
           exchange[i,km,tm]=vel_section[i,km,tm]*(distMeeus(x,y)*corners$MeanDepth[i])

           #add hyperdiffusion fix
	   #hyperdiff[i]=distMeeus(xy[1,],xy[2,])/(corners$Size.km[i]*10e6)
	   hyperdiff[i]=distMeeus(x,y)/(corners$Size.km[i]*10e6)

           #fix flux with hyperdiffuction
           exchange[i,km,tm]=hyperdiff[i]*exchange[i,km,tm]

           dest_k[i,km,tm]=km-1 #NB!!! start couting from 0, NOT from 1
           dest_b[i,km,tm]=corners$neighbour[i]
           boxnr[i,km,tm]=corners$Area[i]
        } # end loop over i
      }#end if over length
    }#end loop over depth
    }#end loop over time

   t_exc[,,cnt:(cnt+Ti-1)]=exchange
   t_db[,,cnt:(cnt+Ti-1)]=dest_b
   t_dk[,,cnt:(cnt+Ti-1)]=dest_k
   t_boxn[,,cnt:(cnt+Ti-1)]=boxnr
   cnt=cnt+Ti


#remove access entries
t_exc=t_exc[,,1:(cnt-1)]
t_db=t_db[,,1:(cnt-1)]
t_dk=t_dk[,,1:(cnt-1)]
t_boxn=t_boxn[,,1:(cnt-1)]

#interpolate between - to get daily values. Each field is a mean over three days
t_fin=array(0,dim=c(length(section),length(h_get),(cnt-1)*daymean))
db_fin=array(NA,dim=c(length(section),length(h_get),(cnt-1)*daymean))
dk_fin=array(NA,dim=c(length(section),length(h_get),(cnt-1)*daymean))
bn_fin=array(NA,dim=c(length(section),length(h_get),(cnt-1)*daymean))


jcnt=1
for (i in 1:(cnt-2)){
for(j in 1:(length(w)-1)){
       t_fin[,,jcnt+j-1]=t_exc[,,i]*w[(length(w)-j+1)]+t_exc[,,i+1]*w[j]  
       db_fin[,,jcnt:(jcnt+j-1)]=t_db[,,i]
       dk_fin[,,jcnt:(jcnt+j-1)]=t_dk[,,i]
       bn_fin[,,jcnt:(jcnt+j-1)]=t_boxn[,,i]
}
       jcnt=jcnt+daymean
}

#length <365 for years will create a wrong seasonal cycle, especially in long runs
if((jcnt-1)<365 & cnt>360){
   t_fin[,,(jcnt:365)]=t_fin[,,jcnt-1]
   db_fin[,,(jcnt:365)]=db_fin[,,jcnt-1]
   dk_fin[,,(jcnt:365)]=dk_fin[,,jcnt-1]
   bn_fin[,,(jcnt:365)]=bn_fin[,,jcnt-1]
}

#temporary arrays
exch_nc=array(0,dim=c(max_neigh,length(h_get),(max(corners$Area)+1),dim(t_fin)[3])) # max neighbours set to 20 for now
db_nc=array(NA,dim=c(max_neigh,length(h_get),(max(corners$Area)+1),dim(t_fin)[3])) # max neighbours set to 20 for now
dk_nc=array(NA,dim=c(max_neigh,length(h_get),(max(corners$Area)+1),dim(t_fin)[3])) # max neighbours set to 20 for now

start_b=0
it_prev=1
tot_bnr=max(corners$Area)+1 #total nr of boxes - start counting on zero, hence need 1 additional
for (i in 1:tot_bnr){
   it=findInterval((i-1),bn_fin[,1,1])
   tot_it=it_prev:it
   for (k in 1:dim(t_fin)[2]){
   for (j in 1:dim(t_fin)[3]){
      exch_nc[1:length(tot_it),k,i,j]=t_fin[it_prev:it,k,j]
      db_nc[1:length(tot_it),k,i,j]=db_fin[it_prev:it,k,j]
      dk_nc[1:length(tot_it),k,i,j]=dk_fin[it_prev:it,k,j]
   }
   }
   it_prev=it+1
}


exch_nc[which(is.na(exch_nc))]=0.0
#fday=days[1]
#f.path=paste("path_to_romsfiles",year,sep="/")
#f.time = paste(f.path,"filename.nc", sep="")
#nc_t=nc_open(f.time)
#time=ncvar_get(nc_t,"ocean_time")
t_start=ocean_time[1]
dt=86400
t_stop=ocean_time[1]+(dim(exch_nc)[4]-1)*86400
t_tot=seq(t_start,t_stop,dt)

filename=paste("NoBa_NorESM_hydro_",year,".nc",sep="")
print(filename)

#define dimensions
dimd=ncdim_def("dest","",1:max_neigh,create_dimvar=FALSE)
dimb=ncdim_def("b","",1:tot_bnr,create_dimvar=FALSE)
dimz=ncdim_def("z","",1:nlay,create_dimvar=FALSE)
dimt=ncdim_def("t1","",1:length(t_tot),unlim=TRUE)#,create_dimvar=FALSE)

#create variables
#NB!!!!!! Unlimited rec needs to be on the right - otherwise R complains!
#origMissVal_ex=0.0
var.t=ncvar_def("t","seconds since 2008-01-01 00:00:00 +10",dimt,0,prec="double")
var.e=ncvar_def("exchange","m/sˆ3",list(dimd,dimz,dimb,dimt),0,prec="double")
var.b=ncvar_def("dest_b","destb",list(dimd,dimz,dimb,dimt),0,prec="integer")
var.k=ncvar_def("dest_k","destk",list(dimd,dimz,dimb,dimt),0,prec="integer")

dt=86400;
nc_transp=nc_create(filename,list(var.t,var.e,var.b,var.k))

#assign global attributes to file
ncatt_put(nc_transp,0,"title","Transport file, Nordic,1998")
ncatt_put(nc_transp,0,"geometry","Nordic.bgm")
ncatt_put(nc_transp,0,"parameters","")

#assign attributes to variables
ncatt_put(nc_transp,var.t,"dt",86400,prec="double")

#assign variables to file
ncvar_put(nc_transp,var.e,exch_nc)
ncvar_put(nc_transp,var.t,t_tot)
ncvar_put(nc_transp,var.b,db_nc)
ncvar_put(nc_transp,var.k,dk_nc)

nc_close(nc_transp)

} # end loop years
