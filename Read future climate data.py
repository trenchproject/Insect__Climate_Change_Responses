#!/usr/bin/python3
 
import cdsapi
import netCDF4 # For more info: https://unidata.github.io/netcdf4-python/
from netCDF4 import num2date
import numpy as np
import os
import pandas as pd
from zipfile import ZipFile


# USER: enter desired latitude and longitude and name of location
lat = 6.45
lon = 2.35
loc = "Benin"
# USER: must also update line 51 after netCDF file downloaded


# Set working directory to the same as the python code
cwd = os.getcwd()
if cwd != '/Users/johnson/Documents/Christopher/Washington/Research/Climate data analyses':
    os.chdir('/Users/johnson/Documents/Christopher/Washington/Research/Climate data analyses')


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


# USER: enter netcdf file name:
nc_max = "tasmax_day_CESM2_ssp370_r4i1p1f1_gn_20150101-21010101_v20200528.nc"
nc_min = nc_max.replace("tasmax","tasmin")


# Open netCDF4 file (USER: must update file name)
# Naming convention: <variable_id>_<table_id>_<source_id>_<experiment_id>_<variant_label>_<grid_label>_<time_range>_<version updata>.nc
f = netCDF4.Dataset(nc_max)

# Extract variable
var = f.variables['tasmax']
 
# Get dimensions assuming 3D: time, latitude, longitude
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

# Delete Climate Store data files
os.rename(nc_max, file_name + ".nc")
#os.remove("adaptor.esgf_wps.retrieve-1632876778.2203593-19762-16-63dfaacc-2e6c-4f07-87b5-975b6152e747_provenance.json")
#os.remove("adaptor.esgf_wps.retrieve-1632876778.2203593-19762-16-63dfaacc-2e6c-4f07-87b5-975b6152e747_provenance.png")



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

# Open netCDF4 file (USER: must update file name)
# Naming convention: <variable_id>_<table_id>_<source_id>_<experiment_id>_<variant_label>_<grid_label>_<time_range>_<version updata>.nc
f = netCDF4.Dataset(nc_min)

# Extract variable
var = f.variables['tasmin']
 
# Get dimensions assuming 3D: time, latitude, longitude
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
os.rename(nc_min, file_name2 + ".nc")
#os.remove("adaptor.esgf_wps.retrieve-1633011791.616639-15237-12-8103cbb3-333e-44f0-84ff-0e5a6f994590_provenance.json")
#os.remove("adaptor.esgf_wps.retrieve-1633011791.616639-15237-12-8103cbb3-333e-44f0-84ff-0e5a6f994590_provenance.png")



'''
# SUPPLEMENTARY CODES: DO NOT RUN THIS CODE OR IT WILL CREATE 1000s OF .CSV FILES AND CRASH THE COMPUTER!!!!!!
output_dir = './'

# =============================== METHOD 1 ============================
# Extract each time as a 2D pandas DataFrame and write it to CSV
# =====================================================================
os.makedirs(output_dir, exist_ok=True)
for i, t in enumerate(times):
    filename = os.path.join(output_dir, f'{t.isoformat()}.csv')
    print(f'Writing time {t} to {filename}')
    df = pd.DataFrame(var[i, :, :], index=latitudes, columns=longitudes)
    df.to_csv(filename)
print('Done')

# =============================== METHOD 2 ============================
# Write data as a table with 4 columns: time, latitude, longitude, value
# =====================================================================
filename = os.path.join(output_dir, 'table.csv')
print(f'Writing data in tabular form to {filename} (this may take some time)...')
times_grid, latitudes_grid, longitudes_grid = [
    x.flatten() for x in np.meshgrid(times, latitudes, longitudes, indexing='ij')]
df = pd.DataFrame({
    'time': [t.isoformat() for t in times_grid],
    'latitude': latitudes_grid,
    'longitude': longitudes_grid,
    'var': var[:].flatten()})
df.to_csv(filename, index=False)
print('Done')
'''