# Alberto Rovellini
# 9/27/2021
# code to interpolate the outputs from the ROMS-to-Atlantis transformation code and bring them on a 12 h time step
# It also changes the time step column from ocean time to numbered time steps
library(tidyverse)

# read data - consider that you may have a number of data input files - how does HC decide which came first? Merge them all together?

file1 <- '../../outputs/state_vars_test.dat'
file2 <- '../../outputs/transport_test.dat'

roms_state_vars <- read.csv(file1,sep='\t')
roms_hydro <- read.csv(file2,sep='\t')

# sometimes the ROMS time series may have gaps. If we ignore them, time steps may end up being shifted
# Write a function that: 
# 1. Takes as arguments the ROMS-to-Atlantis outputs and whether they are state variables
# 1. Builds a sequence of time steps from t0 (the first of the series whenever that would be) and going forward every 12 h
# 2. fills the time series at each 12 h time step by box and depth layer (and face for the fluxes)


fill_time_steps_12h <- function(roms_data,statevars){ 
  ts <- unique(roms_data[,grep('Time',colnames(roms_data))])
  t_0 <- as.POSIXct(ts[1],origin='1900-01-01',tz='UTC') # this has to be specific to your ROMS, so check your origin and tz
  t_end <- as.POSIXct(ts[length(ts)],origin='1900-01-01',tz='UTC')
  complete <- seq(from=t_0,to=t_end,by=60*60*12) # 12 hours is the target for HC, units from ts are in seconds
  
  # change time column to actual times
  roms_data[,grep('Time',colnames(roms_data))] <- as.POSIXct(roms_data[,grep('Time',colnames(roms_data))],origin='1900-01-01',tz='UTC')
  
  boxes <- unique(roms_data$Polygon.number)
  df_all <- as.data.frame(matrix(ncol=ncol(roms_data)))[-1,] # prepare an empty data frame to fill
  
  if(isTRUE(statevars)){ # for state variables (salt, temp, w, etc.)
    for(Ibox in 1:length(boxes)){
      layers <- unique(roms_data[roms_data$Polygon.number==boxes[Ibox],]$Depth.Layer..m.)
      for(Ilyr in 1:length(layers)){
        this_cell <- roms_data[roms_data$Polygon.number==boxes[Ibox] & roms_data$Depth.Layer..m.==layers[Ilyr],]
        # maintaining names as they are in the data
        Vertical.velocity..m3.s. <- approx(this_cell$Time.Step,this_cell$Vertical.velocity..m3.s.,xout=complete,rule=2)$y
        Average.Temperature..Celsius. <- approx(this_cell$Time.Step,this_cell$Average.Temperature..Celsius.,xout=complete,rule=2)$y
        Average.Salinity..PartPer1000. <- approx(this_cell$Time.Step,this_cell$Average.Salinity..PartPer1000.,xout=complete,rule=2)$y
        # put together into a data frame
        df_long <- data.frame('Time.Step'=complete,
                              'Polygon.number'=boxes[Ibox],
                              'Depth.Layer..m.'=layers[Ilyr],
                              Vertical.velocity..m3.s.,
                              Average.Temperature..Celsius.,
                              Average.Salinity..PartPer1000.)

        df_all <- rbind(df_all,df_long)
      }
    }
    df_all <- df_all %>% arrange(Time.Step,Polygon.number,Depth.Layer..m.) %>%
      mutate(Time.Step=(as.numeric(difftime(Time.Step,t_0,units = 'hours'))/12)+1) # replace the true time with a time step counter - seems to be what HC needs

  } else { # ...and for fluxes between boxes
    for(Ibox in 1:length(boxes)){
      faces <- unique(roms_data[roms_data$Polygon.number==boxes[Ibox],]$Face.number)
      for(Iface in 1:length(faces)){
        layers <- unique(roms_data[roms_data$Polygon.number==boxes[Ibox],]$Depth.Layer)
        for(Ilyr in 1:length(layers)){
          this_cell <- roms_data[roms_data$Polygon.number==boxes[Ibox] & 
                                   roms_data$Face.number==faces[Iface] &
                                   roms_data$Depth.Layer==layers[Ilyr],]
          # maintaining names as they are in the data
          Flux..m3.s. <- approx(this_cell$Time.Step..12.hr,this_cell$Flux..m3.s.,xout=complete,rule=2)$y
          # put together into a data frame
          df_long <- data.frame('Polygon.number'=boxes[Ibox],
                                'Face.number'=faces[Iface],
                                'Time.Step..12.hr'=complete,
                                'Depth.Layer'=layers[Ilyr],
                                Flux..m3.s.)
          df_all <- rbind(df_all,df_long)
        }
      }
    }
    df_all <- df_all %>% arrange(Time.Step..12.hr,Polygon.number,Face.number,Depth.Layer) %>%
      mutate(Time.Step..12.hr=(as.numeric(difftime(Time.Step..12.hr,t_0,units = 'hours'))/12)+1) # replace the true time with a time step counter - seems to be what HC needs

  }
  return(df_all)
}

# apply (it is slow with the transport files)

roms_state_vars_interp <- fill_time_steps_12h(roms_data = roms_state_vars,statevars = TRUE)
roms_hydro_interp <- fill_time_steps_12h(roms_data = roms_hydro,statevars = FALSE)

# rename columns to replace '.' with ' ' (probably we can change the ROMS transformation code to skip a couple of these flips back and forth)

colnames(roms_state_vars_interp) <- colnames(read.csv(file1,sep='\t',check.names = FALSE))
colnames(roms_hydro_interp) <- colnames(read.csv(file2,sep='\t',check.names = FALSE))

# export
write.table(roms_state_vars_interp, 'C:/Users/Alberto Rovellini/Documents/GOA/ROMS/HydroConstruct/Inputs/state.dat', 
            quote=FALSE, row.names = FALSE, sep = '\t')

write.table(roms_hydro_interp, 'C:/Users/Alberto Rovellini/Documents/GOA/ROMS/HydroConstruct/Inputs/transport.dat', 
            quote=FALSE, row.names = FALSE, sep = '\t')

