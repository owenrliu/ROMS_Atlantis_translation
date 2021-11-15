#################################################################################################
# Bias-correct GCM projections of SST and CHL from GFDL, IPSL and HADL GCMs
# for the whole North Pacific
# Procedure:
# SST:
#   1. Calculate historical monthly climatology from COBE2 SST (1976-2005 to match Mer's projections)
#   2. Calculate historical climatology from ESM (using same years as obs: 1976-2005)
#   3. Calculate ESM delta: subtract ESM historical climatology from ESM projection to get monthly change
#   4. Add ESM delta to historical climatology from observations
# CHL: (annual ESM output):
#   5. Calculate historical monthly climatology from observations (1997-2020 ESA 4.2)
#   6. Calculate historical mean from ESM (average across all years, using same years as obs: 1997-2020)
#   7. Calculate ESM fractional change: divide ESM projection by ESM historical mean to get annual fractional change
#   8. Multiply historical monthly climatology from observations by ESM fractional change 
#         (same fractional change will be applied to all months of each year)
# Per Mike, use the 2nd option (fractional change) for all biogeochemical variables
##################################################################################################

library(ncdf4)
library(lubridate)
library(reshape2)
library(plyr)
library(tidyverse)

# Define where I store all the large netcdfs
extDrive <- "E:/"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First SST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################################################################################################
# Step 1: load historical SST observations (here I used COBE SST)
# This particular nc is already trimmed to the North Pacific (0-60N, 120-250E),
# and includes years 1900 - 2018
cobeSST <- nc_open("./treaty/COBE2SST_060N_120250E_1900_2018.nc") 
# print(cobeSST) # check the details if you want
latsst <- ncvar_get(cobeSST, "lat")
lonsst <- ncvar_get(cobeSST, "lon")
timesst <- ncvar_get(cobeSST, "time")
datesst <- substr(as.POSIXlt(timesst * 86400, origin = '1891-01-01', tz = "GMT"), 1, 10) 
# Figure out which date range we want. Not subsetting lon or lat in this example
dtMin <- which(datesst == as.Date("1976-01-01"))
dtMax <- which(datesst == as.Date("2005-12-01"))
# Extract SST. Index order is lon, lat, time
sstObs <- ncvar_get(cobeSST, "sst", start = c(1, 1, dtMin), 
                                    count = c(length(lonsst), length(latsst), (dtMax - dtMin) + 1)) 
nc_close(cobeSST)
# Reshape 3D array to a 2d df, easier to work with
dimnames(sstObs) <- list(lon = lonsst, lat = latsst, date = datesst[dtMin:dtMax])
sstlist <- reshape2::melt(sstObs, value.name = "obsHistSST")
sstlist$date <- as.Date(sstlist$date)
sstlist$month <- month(sstlist$date)
sstlist$year <- year(sstlist$date)
sstlist$date <- NULL
# As we're working in the CCS for this project, let's work with longitude in degrees west
sstlist$lon360 <- sstlist$lon
sstlist$lon <- ifelse(sstlist$lon360 > 180, sstlist$lon360 - 360, sstlist$lon360)
# Now calculate COBE2 SST monthly climatology (1976-2005)
sstClim <- aggregate(obsHistSST ~ lon + lat + month, sstlist, FUN = mean, na.rm = TRUE, na.action = na.omit)

###############################################################################################
# Step 2: Next calculate 1976-2005 monthly SST climatology for 3 GCMs
# These are "flooded" so where pixels were missing nearshore, values from the next pixel offshore
# were carried shorewards to fill in any gaps
# First the GFDL ESM (sorry, hard path as on lge external drive)
gfdlSST <- nc_open(
            paste0(extDrive, 
            "projections/1x1_NorthPacificFlooded/NEP_sst_Omon_GFDL-ESM2M_historical_r1i1p1_1976-2005_1x1_ext.nc"))
latgcm <- ncvar_get(gfdlSST, "lat")
longcm <- ncvar_get(gfdlSST, "lon")
# Fun! Each of the ESMs deals with time/date in a different way. GFDL is a 365 (no leap) calendar
timegcm <- ncvar_get(gfdlSST, "time")
# I added 28 days to counteract the leap effect (in a kludgy way). Works on a monthly level
dategcm <- substr(as.POSIXlt((timegcm + 28) * 86400, origin = '1861-01-01', tz = "GMT"), 1, 10) 
# Now extract SST from the GFDL ESM. Index order is lon, lat, depth level, time
# Note that when Mer initially sent me an earlier version of these files, the variable name was still "thetao" 
# (not "sst"), can check using "print"
sstGFDL <- ncvar_get(gfdlSST, "sst", start = c(1, 1, 1, 1), 
                     count = c(length(longcm), length(latgcm), 1, length(timegcm)))
nc_close(gfdlSST)
dimnames(sstGFDL) <- list(lon = longcm, lat = latgcm, date = dategcm)
# Now IPSL. Can just use the same dates as for GFDL
ipslSST <- nc_open(
            paste0(extDrive, 
            "projections/1x1_NorthPacificFlooded/NEP_sst_Omon_IPSL-CM5A-MR_historical_r1i1p1_1976-2005_1x1_ext.nc"))
sstIPSL <- ncvar_get(ipslSST, "sst", start = c(1, 1, 1, 1), 
                     count = c(length(longcm), length(latgcm), 1, length(timegcm)))
nc_close(ipslSST)
dimnames(sstIPSL) <- list(lon = longcm, lat = latgcm, date = dategcm)
# Now Hadley
hadlSST <- nc_open(
            paste0(extDrive, 
            "projections/1x1_NorthPacificFlooded/NEP_sst_Omon_HadGEM2-ES_historical_r2i1p1_1976-2005_1x1_ext.nc"))
sstHADL <- ncvar_get(hadlSST, "sst", start = c(1, 1, 1, 1), 
                     count = c(length(longcm), length(latgcm), 1, length(timegcm)))
nc_close(hadlSST)
dimnames(sstHADL) <- list(lon = longcm, lat = latgcm, date = dategcm)

# Reshape all ESM outputs, add month/year fields, create climatologies
gfdlSSTlist <- reshape2::melt(sstGFDL, value.name = "gfdlHistSST")
ipslSSTlist <- reshape2::melt(sstIPSL, value.name = "ipslHistSST")
hadlSSTlist <- reshape2::melt(sstHADL, value.name = "hadlHistSST")
gfdlSSTlist$date <- ipslSSTlist$date <- hadlSSTlist$date <- as.Date(gfdlSSTlist$date)
gfdlSSTlist$month <- ipslSSTlist$month <- hadlSSTlist$month <- month(gfdlSSTlist$date)
gfdlSSTlist$year <- ipslSSTlist$year <- hadlSSTlist$year <- year(gfdlSSTlist$date)
gfdlSSTlist$date <- ipslSSTlist$date <- hadlSSTlist$date <- NULL
gcmsSSTlist <- data.frame(cbind(gfdlSSTlist, "ipslHistSST" = ipslSSTlist$ipslHistSST, 
                                "hadlHistSST" = hadlSSTlist$hadlHistSST))
sstClimGCMs <- aggregate(cbind(gfdlHistSST, ipslHistSST, hadlHistSST) ~ lon + lat + month, 
                         gcmsSSTlist, FUN = mean, na.rm = TRUE, na.action = na.omit)

# # Check units for longitude (can be deg East or deg West), and SST (can be deg Celcius or deg Kelvin)
# # Has varied depending on versions that Mer has sent me!
# hist(sstClimGCMs$lon)
# hist(sstClimGCMs$gfdlHistSST)

# Join SST from ESMs with SST observations. SST observations only include ocean pixels, so this will have the effect
# of trimming the ESMs (flooded eastwards over land) to the correct extent
sstObsGCMs <- inner_join(sstClim, sstClimGCMs, by = c("lon", "lat", "month"))
# Use a map to check the extent of the data, if you want
pac.coast <- borders("world", colour="gray50", fill="gray50", xlim = c(-134, -110), ylim = c(20, 60))
ggplot(sstObsGCMs) + pac.coast + geom_point(aes(x = lon, y = lat))

################################################################################################
# Step 3: Calculate historical/future delta from ESMs
# First load GFDL future projections
gfdlSSTpr <- nc_open(
              paste0(extDrive, 
              "projections/1x1_NorthPacificFlooded/NEP_sst_Omon_GFDL-ESM2M_rcp85_r1i1p1_2006-2100_1x1_ext.nc"))
timegcmpr <- ncvar_get(gfdlSSTpr, "time")
# I added 9 days to counteract the leap effect again 
dategcmpr <- substr(as.POSIXlt((timegcmpr + 9) * 86400, origin = '2006-01-01', tz = "GMT"), 1, 10)
# Get the GFDL SST projections
sstGFDLpr <- ncvar_get(gfdlSSTpr, "sst", start = c(1, 1, 1, 1), 
                     count = c(length(longcm), length(latgcm), 1, length(timegcmpr)))
nc_close(gfdlSSTpr)
dimnames(sstGFDLpr) <- list(lon = longcm, lat = latgcm, date = dategcmpr)
# Now IPSL
ipslSSTpr <- nc_open(
              paste0(extDrive,
              "projections/1x1_NorthPacificFlooded/NEP_sst_Omon_IPSL-CM5A-MR_rcp85_r1i1p1_2006-2100_1x1_ext.nc"))
sstIPSLpr <- ncvar_get(ipslSSTpr, "sst", start = c(1, 1, 1, 1), 
                       count = c(length(longcm), length(latgcm), 1, length(timegcmpr)))
nc_close(ipslSSTpr)
dimnames(sstIPSLpr) <- list(lon = longcm, lat = latgcm, date = dategcmpr)
# Now HADL
hadlSSTpr <- nc_open(
              paste0(extDrive,
              "projections/1x1_NorthPacificFlooded/NEP_sst_Omon_HadGEM2-ES_rcp85_r2i1p1_2006-2100_1x1_ext.nc"))
sstHADLpr <- ncvar_get(hadlSSTpr, "sst", start = c(1, 1, 1, 1), # was thetao before
                       count = c(length(longcm), length(latgcm), 1, length(timegcmpr)))
nc_close(hadlSSTpr)
dimnames(sstHADLpr) <- list(lon = longcm, lat = latgcm, date = dategcmpr)

# Melt all to 2D
gfdlSSTlistpr <- reshape2::melt(sstGFDLpr, value.name = "gfdlFutSST")
ipslSSTlistpr <- reshape2::melt(sstIPSLpr, value.name = "ipslFutSST")
hadlSSTlistpr <- reshape2::melt(sstHADLpr, value.name = "hadlFutSST")
gfdlSSTlistpr$date <- ipslSSTlistpr$date <- hadlSSTlistpr$date <- as.Date(gfdlSSTlistpr$date)
gfdlSSTlistpr$month <- ipslSSTlistpr$month <- hadlSSTlistpr$month <- month(gfdlSSTlistpr$date)
gfdlSSTlistpr$year <- ipslSSTlistpr$year <- hadlSSTlistpr$year <- year(gfdlSSTlistpr$date)
gfdlSSTlistpr$date <- ipslSSTlistpr$date <- hadlSSTlistpr$date <- NULL
gcmsSSTlistpr <- data.frame(cbind(gfdlSSTlistpr, "ipslFutSST" = ipslSSTlistpr$ipslFutSST, 
                                "hadlFutSST" = hadlSSTlistpr$hadlFutSST))
# Check SST units again. Whoops, GFDL is still in degrees Kelvin
hist(gcmsSSTlistpr$gfdlFutSST)
gcmsSSTlistpr$gfdlFutSST <- gcmsSSTlistpr$gfdlFutSST - 273.15

# Add in historical and observed SSTs from sstObsGCMs created in steps 1-2
gcmsSSTlistpr <- left_join(gcmsSSTlistpr, sstObsGCMs, by = c("lon", "lat", "month"))

###############################################################################################
# Step 4: calculate deltas and add to observations. ("BC" denotes bias-corrected)
gcmsSSTlistpr$gfdlBCsst <- (gcmsSSTlistpr$gfdlFutSST - gcmsSSTlistpr$gfdlHistSST) + gcmsSSTlistpr$obsHistSST
gcmsSSTlistpr$ipslBCsst <- (gcmsSSTlistpr$ipslFutSST - gcmsSSTlistpr$ipslHistSST) + gcmsSSTlistpr$obsHistSST
gcmsSSTlistpr$hadlBCsst <- (gcmsSSTlistpr$hadlFutSST - gcmsSSTlistpr$hadlHistSST) + gcmsSSTlistpr$obsHistSST

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Now chlorophyll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################################################################################################
# Step 5: Load chlorophyll: using ESA satellite reanalysis which was coarsened to 1x1 degrees
# Satellite chl only available 1997-present (2020), so that's what I used for my "historical" period
chlClim <- readRDS(paste0(extDrive, "/esa/monthly/esaChlMonthly42_060N_140250E_clim_19972020_1x1.rds"))
# As above: we're working in degW this time, so change lon chlClim
chlClim$lon360 <- chlClim$lon
chlClim$lon <- ifelse(chlClim$lon360 > 180, chlClim$lon360 - 360, chlClim$lon360)
###############################################################################################
# Step 6: Calculate historical mean CHL from ESMs 
gfdlCHL <- nc_open(
            paste0(extDrive,
            "projections/1x1_NorthPacificFlooded/NEP_surfchl_Oyr_GFDL-ESM2M_historical_r1i1p1_1976-2005_1x1_ext.nc"))
latgcm <- ncvar_get(gfdlCHL, "lat") # Same as SST netcdfs above
longcm <- ncvar_get(gfdlCHL, "lon")
timegcmCHL <- ncvar_get(gfdlCHL, "time") # Different to SST: BGC variables are output annually, not monthly
dategcmCHL <- substr(as.POSIXlt((timegcmCHL) * 86400, origin = '1861-01-01', tz = "GMT"), 1, 10) 
chlGFDL <- ncvar_get(gfdlCHL, "chl", start = c(1, 1, 1, 1), 
                     count = c(length(longcm), length(latgcm), 1, length(dategcmCHL)))
nc_close(gfdlCHL)
dimnames(chlGFDL) <- list(lon = longcm, lat = latgcm, date = dategcmCHL)
# Now IPSL
ipslCHL <- nc_open(
          paste0(extDrive,
          "projections/1x1_NorthPacificFlooded/NEP_surfchl_Oyr_IPSL-CM5A-MR_historical_r1i1p1_1976-2005_1x1_ext.nc"))
chlIPSL <- ncvar_get(ipslCHL, "chl", start = c(1, 1, 1, 1), 
                     count = c(length(longcm), length(latgcm), 1, length(dategcmCHL)))
nc_close(ipslCHL)
dimnames(chlIPSL) <- list(lon = longcm, lat = latgcm, date = dategcmCHL)
# Now HADL
hadlCHL <- nc_open(
           paste0(extDrive,
           "projections/1x1_NorthPacificFlooded/NEP_surfchl_Oyr_HadGEM2-ES_historical_r2i1p1_1976-2005_1x1_ext.nc"))
chlHADL <- ncvar_get(hadlCHL, "chl", start = c(1, 1, 1, 1), 
                     count = c(length(longcm), length(latgcm), 1, length(dategcmCHL)))
nc_close(hadlCHL)
dimnames(chlHADL) <- list(lon = longcm, lat = latgcm, date = dategcmCHL)

# Melt all, create climatologies
gfdlCHLlist <- reshape2::melt(chlGFDL, value.name = "gfdlHistCHL")
ipslCHLlist <- reshape2::melt(chlIPSL, value.name = "ipslHistCHL")
hadlCHLlist <- reshape2::melt(chlHADL, value.name = "hadlHistCHL")
gfdlCHLlist$date <- ipslCHLlist$date <- hadlCHLlist$date <- as.Date(gfdlCHLlist$date)
gfdlCHLlist$year <- ipslCHLlist$year <- hadlCHLlist$year <- year(gfdlCHLlist$date)
gfdlCHLlist$date <- ipslCHLlist$date <- hadlCHLlist$date <- NULL
gcmsCHLlist <- data.frame(cbind(gfdlCHLlist, "ipslHistCHL" = ipslCHLlist$ipslHistCHL, 
                                "hadlHistCHL" = hadlCHLlist$hadlHistCHL))
# Then aggregate
chlClimGCMs <- aggregate(cbind(gfdlHistCHL, ipslHistCHL, hadlHistCHL) ~ lon + lat, 
                         gcmsCHLlist, FUN = mean, na.rm = TRUE, na.action = na.omit)

# Convert chl from kg/m3 to mg/m3
chlClimGCMs$gfdlHistCHL <- chlClimGCMs$gfdlHistCHL * 1e+6 
chlClimGCMs$ipslHistCHL <- chlClimGCMs$ipslHistCHL * 1e+6
chlClimGCMs$hadlHistCHL <- chlClimGCMs$hadlHistCHL * 1e+6

# Join with CHL observations
# Note that chlClim is by month, but chlClimGCMs is only available by year
chlObsGCMs <- left_join(chlClim, chlClimGCMs, by = c("lon", "lat")) 

################################################################################################
# Step 7: Divide ESM projection by ESM historical mean to get annual fractional change
# First load GFDL projections
gfdlCHLpr <- nc_open(
  "E:/projections/1x1_NorthPacificFlooded/NEP_surfchl_Oyr_GFDL-ESM2M_rcp85_r1i1p1_2006-2100_1x1_ext.nc")
timegcmCHLpr <- ncvar_get(gfdlCHLpr, "time")
dategcmCHLpr <- substr(as.POSIXlt((timegcmCHLpr) * 86400, origin = '2006-01-01', tz = "GMT"), 1, 10) 
chlGFDLpr <- ncvar_get(gfdlCHLpr, "chl", start = c(1, 1, 1, 1), 
                       count = c(length(longcm), length(latgcm), 1, length(timegcmCHLpr)))
nc_close(gfdlCHLpr)
dimnames(chlGFDLpr) <- list(lon = longcm, lat = latgcm, date = dategcmCHLpr)
# Now IPSL
ipslCHLpr <- nc_open(
  "E:/projections/1x1_NorthPacificFlooded/NEP_surfchl_Oyr_IPSL-CM5A-MR_rcp85_r1i1p1_2006-2100_1x1_ext.nc")
chlIPSLpr <- ncvar_get(ipslCHLpr, "chl", start = c(1, 1, 1, 1), 
                       count = c(length(longcm), length(latgcm), 1, length(timegcmCHLpr)))
nc_close(ipslCHLpr)
dimnames(chlIPSLpr) <- list(lon = longcm, lat = latgcm, date = dategcmCHLpr)
# Now HADL
hadlCHLpr <- nc_open(
  "E:/projections/1x1_NorthPacificFlooded/NEP_surfchl_Oyr_HadGEM2-ES_rcp85_r2i1p1_2006-2100_1x1_ext.nc")
chlHADLpr <- ncvar_get(hadlCHLpr, "chl", start = c(1, 1, 1, 1), 
                       count = c(length(longcm), length(latgcm), 1, length(timegcmCHLpr)))
nc_close(hadlCHLpr)
dimnames(chlHADLpr) <- list(lon = longcm, lat = latgcm, date = dategcmCHLpr)

# Melt all, create climatologies
gfdlCHLlistpr <- reshape2::melt(chlGFDLpr, value.name = "gfdlFutCHL")
ipslCHLlistpr <- reshape2::melt(chlIPSLpr, value.name = "ipslFutCHL")
hadlCHLlistpr <- reshape2::melt(chlHADLpr, value.name = "hadlFutCHL")
gfdlCHLlistpr$date <- ipslCHLlistpr$date <- hadlCHLlistpr$date <- as.Date(gfdlCHLlistpr$date)
gfdlCHLlistpr$year <- ipslCHLlistpr$year <- hadlCHLlistpr$year <- year(gfdlCHLlistpr$date)
gfdlCHLlistpr$date <- ipslCHLlistpr$date <- hadlCHLlistpr$date <- NULL
gcmsCHLlistpr <- data.frame(cbind(gfdlCHLlistpr, "ipslFutCHL" = ipslCHLlistpr$ipslFutCHL, 
                                  "hadlFutCHL" = hadlCHLlistpr$hadlFutCHL))
# Expand gcmsCHLlistpr 12 times to represent 12 months
# (Because BGC projections are only available at annual timestep)
gcmsCHLlistpr2 <- gcmsCHLlistpr[rep(seq_len(nrow(gcmsCHLlistpr)), 12), ]
gcmsCHLlistpr2$month <- rep(1:12, each = nrow(gcmsCHLlistpr))

# Convert chl from kg/m3 to mg/m3
gcmsCHLlistpr2$gfdlFutCHL <- gcmsCHLlistpr2$gfdlFutCHL * 1e+6
gcmsCHLlistpr2$ipslFutCHL <- gcmsCHLlistpr2$ipslFutCHL * 1e+6
gcmsCHLlistpr2$hadlFutCHL <- gcmsCHLlistpr2$hadlFutCHL * 1e+6

# Now join to chl observations, and historical ESM fields
gcmsCHLlistpr3 <- left_join(chlObsGCMs, gcmsCHLlistpr2, by = c("lon", "lat", "month"))
gcmsCHLlistpr3 <- subset(gcmsCHLlistpr3, !is.na(gfdlHistCHL)) # Just removing some points outside CCS

###############################################################################################
# Step 8: calculate deltas (fractional for chl) and add to observations
# First correct negative values (only in Hadley)
gcmsCHLlistpr3$hadlHistCHL <- ifelse(gcmsCHLlistpr3$hadlHistCHL < 0.001, 0.001, gcmsCHLlistpr3$hadlHistCHL)
gcmsCHLlistpr3$hadlFutCHL <- ifelse(gcmsCHLlistpr3$hadlFutCHL < 0.001, 0.001, gcmsCHLlistpr3$hadlFutCHL)
# Now calculate deltas. ("BC" denotes bias-corrected)
gcmsCHLlistpr3$gfdlBCchl <- (gcmsCHLlistpr3$gfdlFutCHL / gcmsCHLlistpr3$gfdlHistCHL) * gcmsCHLlistpr3$obsChl
gcmsCHLlistpr3$ipslBCchl <- (gcmsCHLlistpr3$ipslFutCHL / gcmsCHLlistpr3$ipslHistCHL) * gcmsCHLlistpr3$obsChl
gcmsCHLlistpr3$hadlBCchl <- (gcmsCHLlistpr3$hadlFutCHL / gcmsCHLlistpr3$hadlHistCHL) * gcmsCHLlistpr3$obsChl
# # Check the values. (4th root is a pretty common transformation for chl, just allows us to see better)
# hist(gcmsCHLlistpr3$gfdlBCchl ^ 0.25)
# hist(gcmsCHLlistpr3$ipslBCchl ^ 0.25)
# hist(gcmsCHLlistpr3$hadlBCchl ^ 0.25)
# # Note the actual chl values in Hadley are reasonable, but the increases in the SW corner of the ROMS
# domain can be proportionally very large. This shouldn't be such an issue for the coastal Atlantis

# Combine SST and CHL bias-corrected projections (chl extent is smaller)
gcmsListBoth <- left_join(gcmsCHLlistpr3, gcmsSSTlistpr, by = c("lon", "lat", "month", "year"))
# Use SST to remove a few points on land (usually inland water bodies w chl obs)
gcmsListBoth <- subset(gcmsListBoth, !is.na(obsHistSST))

# Export if you want
saveRDS(gcmsListBoth, paste0(extDrive, 
                             "projections/1x1_NorthPacificFlooded/biasCorr/biasCorrSST_CHL_3ESMs_start2006.rds"))