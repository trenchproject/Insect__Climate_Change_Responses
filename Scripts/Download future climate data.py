###############################################################################
########## This script downloads CMIP6 climate data for a population ##########
###############################################################################

# Install required packages
import cdsapi # please read "ReadMe install python.docx" to learn how to download the CDS API module
import sys # needed to download netCDF4
import subprocess # needed to download netCDF4
if not 'netCDF4' in sys.modules:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'netCDF4']) # install netCDF4 if not installed
import netCDF4 # For more info: https://unidata.github.io/netcdf4-python/
from netCDF4 import num2date
import numpy as np
import os
import pandas as pd
from zipfile import ZipFile


# USER: please set path to desired folder for downloaded climate files
dir_path = '/Users/c.johnson/Documents/GitHub/Johnson_Insect_Responses/Climate data'

# USER: please enter desired latitude and longitude and name of location from "Climate station data.xlsx"
lat = 6.45 # from "Lat" column of "Climate station data.xlsx" (NOT "Latitude" column)
lon = 2.35 # from "Lon" column of "Climate station data.xlsx"
loc = "Benin" # from "Location" column of "Climate station data.xlsx"


# Set working directory to the same as the python code
cwd = os.getcwd()
if cwd != dir_path:
    os.chdir(dir_path)

# Get names of downloaded files
nc_max = "tasmax_day_CESM2_ssp370_r4i1p1f1_gn_20150101-21010101_v20200528.nc"
nc_min = nc_max.replace("tasmax","tasmin")


# FOR DAILY MAXIMUM TEMPERATURE
# Retrieve data and store as netCDF4 file
c = cdsapi.Client()
file_name = "Future Tmax " + loc

c.retrieve(
    'projections-cmip6',
    {
        'format': 'zip',
        'temporal_resolution': 'daily',
        'experiment': 'ssp3_7_0',
        'level': 'single_levels',
        'variable': 'daily_maximum_near_surface_air_temperature',
        'model': 'cesm2',
        'area': [ lat+0.5, lon-1.25/2, lat-0.5, lon+1.25/2, ], # create margins around desired lat/lon
        # NOTE: WARNING Recovering from HTTP error [500 Internal Server Error], try refreshing browser or logging into https://cds.climate.copernicus.eu/#!/home
    },
    file_name + '.zip')

# Open zip file, extract files, and then delete zip file
with ZipFile(file_name + '.zip', 'r') as zip:
    zip.extractall()
os.remove(file_name + '.zip')

# Open netCDF4 file
f = netCDF4.Dataset(nc_max) # naming convention: <variable_id>_<table_id>_<source_id>_<experiment_id>_<variant_label>_<grid_label>_<time_range>_<version updata>.nc

# Extract variable
var = f.variables['tasmax']
 
# Get dimensions: time, latitude, longitude
time_dim, lat_dim, lon_dim = var.get_dims()
time_var = f.variables[time_dim.name]
times = num2date(time_var[:], time_var.units)
latitudes = f.variables[lat_dim.name][:]
longitudes = f.variables[lon_dim.name][:]

# Extract data and save to CSV
times_grid, latitudes_grid, longitudes_grid = [
    x.flatten() for x in np.meshgrid(times, latitudes, longitudes, indexing='ij')]
df = pd.DataFrame({
    'time': [t.isoformat() for t in times_grid],
    'latitude': latitudes_grid,
    'longitude': longitudes_grid,
    'Tmax': var[:].flatten()})
df.to_csv(file_name + ".csv")
f.close()


# FOR DAILY MINIMUM TEMPERATURE
# Retrieve data and store as netCDF4 file
c = cdsapi.Client()
file_name2 = "Future Tmin " + loc

c.retrieve(
    'projections-cmip6',
    {
        'format': 'zip',
        'temporal_resolution': 'daily',
        'experiment': 'ssp3_7_0',
        'level': 'single_levels',
        'variable': 'daily_minimum_near_surface_air_temperature',
        'model': 'cesm2',
        'area': [ lat+0.5, lon-1.25/2, lat-0.5, lon+1.25/2, ], # create margins around desired lat/lon
        # NOTE: WARNING Recovering from HTTP error [500 Internal Server Error], try refreshing browser or logging into https://cds.climate.copernicus.eu/#!/home
    },
    file_name2 + '.zip')

# Open zip file, extract files, and then delete zip file
with ZipFile(file_name2 + '.zip', 'r') as zip:
    zip.extractall()
os.remove(file_name2 + '.zip')

# Open netCDF4 file
f = netCDF4.Dataset(nc_min) # naming convention: <variable_id>_<table_id>_<source_id>_<experiment_id>_<variant_label>_<grid_label>_<time_range>_<version updata>.nc

# Extract variable
var = f.variables['tasmin']
 
# Get dimensions: time, latitude, longitude
time_dim, lat_dim, lon_dim = var.get_dims()
time_var = f.variables[time_dim.name]
times = num2date(time_var[:], time_var.units)
latitudes = f.variables[lat_dim.name][:]
longitudes = f.variables[lon_dim.name][:]

# Extract data and save to CSV
times_grid, latitudes_grid, longitudes_grid = [
    x.flatten() for x in np.meshgrid(times, latitudes, longitudes, indexing='ij')]
df = pd.DataFrame({
    'time': [t.isoformat() for t in times_grid],
    'latitude': latitudes_grid,
    'longitude': longitudes_grid,
    'Tmin': var[:].flatten()})
df.to_csv(file_name2 + ".csv")
f.close()

# Delete Climate Store data files
os.remove(nc_max)
os.remove(nc_min)
for item in os.listdir(dir_path): # remove JSON and PNG files associated with NC file
    if item.endswith(".json"):
        os.remove(os.path.join(dir_path, item))
    if item.endswith(".png"):
        os.remove(os.path.join(dir_path, item))