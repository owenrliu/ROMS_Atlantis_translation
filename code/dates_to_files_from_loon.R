# Alberto Rovellini
# 3/6/2023
# This code reads each file output from loon after Script 1 and it maps the dates for each year to specific files
# It will be then read in by Script 2 to read the correct files for each year, stitching different 

library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)

# define origin and tz
this_origin <- '1900-01-01'
this_tz <- 'UTC'

# list all state var files
statevar_files <- list.files('../Script1/hindcast_revised/', recursive = T, pattern = 'state_vars', full.names = T)

# list all hydro files
transport_files <- list.files('../Script1/hindcast_revised/', recursive = T, pattern = 'transport', full.names = T)

statevar_list <- list()
for(i in 1:length(statevar_files)){
  
  this_statevar_file <- statevar_files[i]
  
  this_dat <- read.csv(this_statevar_file, sep='\t', header = T)
  
  this_ts <- this_dat %>% 
    mutate(Time.Step=as.POSIXct(Time.Step, origin=this_origin,tz=this_tz)) %>% 
    select(Time.Step) %>% 
    distinct() %>%
    mutate(year = year(Time.Step), month = month(Time.Step), day = day(Time.Step))
  
  # pull dates that have Jan 1st and Dec 31st and add file name
  date_limits <- this_ts %>%
    filter(month == 1 & day == 1 | month == 12 & day == 31) %>%
    mutate(filename = this_statevar_file)
  
  statevar_list[[i]] <- date_limits
  
}

statevar_frame <- rbindlist(statevar_list)

# same for hydro files

transport_list <- list()
for(i in 1:length(transport_files)){
  
  this_transport_file <- transport_files[i]
  
  this_dat <- read.csv(this_transport_file, sep='\t', header = T)
  
  this_ts <- this_dat %>% 
    mutate(Time.Step..12.hr=as.POSIXct(Time.Step..12.hr, origin=this_origin,tz=this_tz)) %>% 
    select(Time.Step..12.hr) %>% 
    rename(Time.Step = Time.Step..12.hr) %>%
    distinct() %>%
    mutate(year = year(Time.Step), month = month(Time.Step), day = day(Time.Step))
  
  # pull dates that have Jan 1st and Dec 31st and add file name
  date_limits <- this_ts %>%
    filter(month == 1 & day == 1 | month == 12 & day == 31) %>%
    mutate(filename = this_transport_file)
  
  transport_list[[i]] <- date_limits
  
}

transport_frame <- rbindlist(transport_list)

# add missing dates
# need to go into the raw files and check the time steps...
incomplete_years <- statevar_frame %>% group_by(year) %>% summarise(n = n()) %>% ungroup() %>% filter(n < 2) %>% pull(year)
incomplete_years
# incomplete in hindcast_revised: 1990 1994 2010 2020

# 1990: leave it we did not translate
# 1994: 1994-12-28 12:00:00 files: ../Script1/hindcast_revised//out9/state_vars_test.dat, ../Script1/hindcast_revised//out9/transport_test.dat
# 2010: 2010-12-27 12:00:00 files: ../Script1/hindcast_revised//out11/state_vars_test.dat, ../Script1/hindcast_revised//out11/transport_test.dat
# 2020: 2020-12-24 12:00:00 files: ../Script1/hindcast_revised//out15/state_vars_test.dat, ../Script1/hindcast_revised//out15/transport_test.dat

statevar_fill <- data.frame("Time.Step" = as.POSIXct(c("1994-12-28 12:00:00", "2010-11-26 12:00:00", "2020-12-24 12:00:00")),
                            "year" = c(1994, 2010, 2020),
                            "month" = c(12, 12, 12),
                            "day" = c(28, 27, 24),
                            "filename" = c("../Script1/hindcast_revised//out9/state_vars_test.dat",
                                           "../Script1/hindcast_revised//out11/state_vars_test.dat",
                                           "../Script1/hindcast_revised//out15/state_vars_test.dat"))
transport_fill <- data.frame(statevar_fill[,-ncol(statevar_fill)], "filename" = c("../Script1/hindcast_revised//out9/transport_test.dat",
                                                                                  "../Script1/hindcast_revised//out11/transport_test.dat",
                                                                                  "../Script1/hindcast_revised//out15/transport_test.dat"))

# bind
statevar_frame <- rbind(statevar_frame, statevar_fill)
transport_frame <- rbind(transport_frame, transport_fill)

# write out
write.csv(statevar_frame, 'statevar_files_list.csv', row.names = F)
write.csv(transport_frame, 'transport_files_list.csv', row.names = F)