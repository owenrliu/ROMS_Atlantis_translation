## Convert CalCurrentV3_utm.bgm to a line shapefile for use in translating ROMS fluxes
# Using R package 'rbgm' from Michael Sumner
#https://cran.r-project.org/web/packages/rbgm/
# install.packages('rbgm')
library(rbgm)
library(dplyr)
library(sf)
