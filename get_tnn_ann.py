import numpy as np
import pandas as pd
import xarray as xr
from xarray import DataArray
from numpy import meshgrid, deg2rad, gradient, cos
import warnings
import timeit
import sys
import getopt
import argparse
import datetime

# Start stopwatch
start = timeit.default_timer()

def readOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="The parsing commands lists.")
    parser.add_argument("-f", "--file", help="Type your file")
    parser.add_argument("-o", "--outdir", help="Type your output directory")
    opts = parser.parse_args(args)
    return opts

# Call the function to read the argument values
options = readOptions(sys.argv[1:])

# import file(s)
ds = xr.open_dataset(options.file)

# Calculate TNn
ds_tnn = ds.resample(time='1Y').min()

# Change date format to just year
ds_tnn['time'] = ds_tnn.time.dt.strftime("%Y")

# Rename variable to TNn
ds_tnn = ds_tnn.rename({'tasmin':'tnnETCCDI'})

# Convert kelvin to celcius
ds_tnn['tnnETCCDI'] = ds_tnn['tnnETCCDI'] - 273.15

# Calculation for grid cell area
def area_grid(lat, lon):
    
    """
    Calculate the area of each grid cell
    (lat, lon)--> (lat, lon), grid-cell area (m**2)
    """

    x, y = meshgrid(lon, lat)
    R = radius(y)
    rlon = deg2rad(gradient(x, axis=1))
    rlat = deg2rad(gradient(y, axis=0))
    rx = R*rlon*cos(deg2rad(y))
    ry = R*rlat
  
    area = ry*rx

    data_array = DataArray(
        area, dims=["lat", "lon"], 
        coords={"lat": lat, "lon": lon})
    return data_array


# Note: this script **does not** make the assumption that Earth is a perfect sphere. It is more accurately defined as a spheroid below
def radius(lat):
    
    '''
    calculate radius of Earth using lat (degrees)
    lat--> radius(m)
    '''
    
    # Since Earth isn't a perfect sphere, must define the spheroid (sphere-like object)
    a = 6378137.0 #equatorial radius (m)
    b = 6356752.3142 #radius at poles (m)
    e2 = 1 - (b**2/a**2)
    
    # convert to geocentric
    lat = deg2rad(lat)
    lat_gc = np.arctan( (1-e2)*np.tan(lat) )

    # calculate radius
    r = ((a * (1 - e2)**0.5)/(1 - (e2 * np.cos(lat_gc)**2))**0.5)

    return r


# Compute weighting where 'ds' is the name of the dataset. replace ds with your dataset name. 
# Area dataArray
array_area = area_grid(ds_tnn['lat'], ds_tnn['lon'])
# Total area
total_area = array_area.sum(['lat','lon'])
# Area-weighted temp
temp_weighted = (ds_tnn['tnnETCCDI']*array_area)/total_area
# weighting
weight_tot = (array_area/total_area).sum(('lat','lon'))
# global average. If you are only interested in the global mean temp, the calculation is simple. Just sum all weighted values to get the average.
global_mean = temp_weighted.sum(dim=('lat','lon'))

# if you want to add the global means to the original dataset

# Add global mean to dataset
ds_out = ds_tnn.assign({'gm':global_mean})

# Timestamp for metadata
creation_date = str('Created by NASA GISS ' + 'on ' +  str(datetime.datetime.now()))

# Add in metadata
ds_out.tnnETCCDI.attrs = {'units': 'degrees_C',
 'long_name': 'Annual Minimum of Daily Minimum Temperature',
 'comment': 'includes both liquid and solid phases',
 'cell_methods': 'area: time: mean',
 'cell_measures': 'area: areacella',
 'history': creation_date}

ds_out.gm.attrs = {'units': 'degrees_C',
'long_name': 'Global Mean of Annual Minimum of Daily Minimum Temperature',
'history': creation_date}


ds_out.tnnETCCDI.encoding = {'unlimited_dims': {'time'},
 'dtype': 'float32'}

ds_out.gm.encoding = {'unlimited_dims': {'time'},
 'dtype': 'float32'}

# convert time to datetime
try:
  ds_out['time'] = pd.to_datetime(ds_out.time, format='%Y')
except:
  None

# Create new dataset with correct name
model = options.file.split('tasmin_day_')[1].split('_')[0]
scenario = options.file.split('tasmin_day_')[1].split('_')[1].split('_')[0]
variant = options.file.split('tasmin_day_')[1].split('_')[2]
years = options.file.split('tasmin_day_')[1].split('_')[4].split('.nc')[0]
ds_out.to_netcdf(options.outdir + '/tnnETCCDI_yr_' + model + '_' + scenario + '_' + variant + '_' + years + '.nc')

# Pause stopwatch
stop = timeit.default_timer()
# Get runtime
print('Runtime: ' + str(stop - start) + ' seconds')
