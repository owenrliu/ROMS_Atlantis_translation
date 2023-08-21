#!/bin/bash

# this script loops over monthly files from the ROMS translation and:
# 1. Concatenates them into annual files along the time dimension
# 2. Turns them to CDFs
# 3. Add the missing lines for the timestep size and the parameters
# 4. Converts back to .nc
# 5. Moves all the annual files into a folder to be updated to the server where Atlantis is ran

# Note that this requires cdo to be installed on the machine where the operation is performed

# mergetime for all folders
# hydro
for year in {1991..2020}; do
    folder="${year}"
    echo "Merging files in folder: $folder"
    
    # Navigate to the subdirectory within the folder
    cd "${folder}/forcings/hydro_hd/"
    
    # Merge NetCDF files using cdo mergetime
    cdo mergetime *.nc goa_hydro_${year}.nc
    
    # Navigate back to the parent directory
    cd ../../../  # Navigate back to the parent folder
done

for year in {1991..2020}; do
    folder="${year}"
    echo "Merging files in folder: $folder"
    
    # Navigate to the subdirectory within the folder
    cd "${folder}/forcings/temp/"
    
    # Merge NetCDF files using cdo mergetime
    cdo mergetime *.nc goa_roms_temp_${year}.nc
    
    # Navigate back to the parent directory
    cd ../../../  # Navigate back to the parent folder
done

# salt
for year in {1991..2020}; do
    folder="${year}"
    echo "Merging files in folder: $folder"
    
    # Navigate to the subdirectory within the folder
    cd "${folder}/forcings/salt/"
    
    # Merge NetCDF files using cdo mergetime
    cdo mergetime *.nc goa_roms_salt_${year}.nc
    
    # Navigate back to the parent directory
    cd ../../../  # Navigate back to the parent folder
done

# ncdump for all annual files to add missing variables
# hydro
for year in {1991..2020}; do
    folder="${year}"
    echo "Converting to CDF in folder: $folder"
    
    # Navigate to the subdirectory within the folder
    cd "${folder}/forcings/hydro_hd/"
    
    # Merge NetCDF files using cdo mergetime
    ncdump goa_hydro_${year}.nc > goa_hydro_${year}.cdf
    
    # Navigate back to the parent directory
    cd ../../../  # Navigate back to the parent folder
done

# temp
for year in {1991..2020}; do
    folder="${year}"
    echo "Converting to CDF in folder: $folder"
    
    # Navigate to the subdirectory within the folder
    cd "${folder}/forcings/temp/"
    
    # Merge NetCDF files using cdo mergetime
    ncdump goa_roms_temp_${year}.nc > goa_roms_temp_${year}.cdf
    
    # Navigate back to the parent directory
    cd ../../../  # Navigate back to the parent folder
done

# salt
for year in {1991..2020}; do
    folder="${year}"
    echo "Converting to CDF in folder: $folder"
    
    # Navigate to the subdirectory within the folder
    cd "${folder}/forcings/salt/"
    
    # Merge NetCDF files using cdo mergetime
    ncdump goa_roms_salt_${year}.nc > goa_roms_salt_${year}.cdf
    
    # Navigate back to the parent directory
    cd ../../../  # Navigate back to the parent folder
done

# add the lines we need. The line number will be different for hydro and temp/salt
# hydro
for file in $(find . -type f -name "goa_hydro_*.cdf"); do
   sed -i "13i\		t:dt = 43200. ;" "$file"
   sed -i "32i\		:parameters = "" ;" "$file"
done

# temp
for file in $(find . -type f -name "goa_roms_temp_*.cdf"); do
   sed -i "12i\		t:dt = 43200. ;" "$file"
   sed -i "24i\		:parameters = "" ;" "$file"
done

# salt
for file in $(find . -type f -name "goa_roms_salt_*.cdf"); do
   sed -i "12i\		t:dt = 43200. ;" "$file"
   sed -i "24i\		:parameters = "" ;" "$file"
done

# now turn back to .nc, all files at the same time
for file in $(find . -type f -name "*.cdf"); do
    output_file="${file%.cdf}.nc"  # Replace .cdf with .nc in the output filename
    ncgen -o "$output_file" "$file"
done

# now copy all to dedicated folders
for year in {1991..2020}; do
    folder="${year}"
    echo "moving files: $folder"
    
    # find files
    hydro_file=$(find . -type f -name "goa_hydro_${year}.nc")
    temp_file=$(find . -type f -name "goa_roms_temp_${year}.nc")
    salt_file=$(find . -type f -name "goa_roms_salt_${year}.nc")

    # move
    mv "$hydro_file" hydro_forcings_revised/hydro
    mv "$temp_file" hydro_forcings_revised/temp
    mv "$salt_file" hydro_forcings_revised/salt
done



