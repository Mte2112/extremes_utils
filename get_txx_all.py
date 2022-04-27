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
    parser.add_argument("-f", "--file", help="Type your file (input should be monthly)")
    parser.add_argument("-o", "--outdir", help="Type your output directory")
    opts = parser.parse_args(args)
    return opts

# Call the function to read the argument values
options = readOptions(sys.argv[1:])

# import file(s)
ds = xr.open_dataset(options.file)
ds_ann = ds
ds_djf = ds.where(ds['time.season']=='DJF', drop=True)
ds_jja = ds.where(ds['time.season']=='JJA', drop=True)

# Calculate txx
txx_ann = ds_ann.resample(time='1YS').max()
# XArray does not have an easy was to resample season, choose 'DJF', and drop others so my workaround is to take the monthly maximums, calculate the rolling 3month max (centered), and drop all non-january months 
txx_djf = ds_djf.resample(time='1MS').max().rolling(time=3, center=True).max()
txx_djf = txx_djf.where(txx_djf['time.month']==1, drop=True)
# The 1-year resampling should work here bc JJA does not span years, but do the same procedure as DJF for consistency
txx_jja = ds_jja.resample(time='1MS').max().rolling(time=3, center=True).max()
txx_jja = txx_jja.where(txx_jja['time.month']==7, drop=True)

# Drop time_bnds for seasonal
txx_djf = txx_djf.drop_vars('time_bnds')
txx_jja = txx_jja.drop_vars('time_bnds')

# Change date format to just year
txx_ann['time'] = txx_ann.time.dt.strftime("%Y")
txx_djf['time'] = txx_djf.time.dt.strftime("%Y")
txx_jja['time'] = txx_jja.time.dt.strftime("%Y")

# Rename variable to txx
txx_ann = txx_ann.rename({'tasmax':'txxETCCDI'})
txx_djf= txx_djf.rename({'tasmax':'txxETCCDI'})
txx_jja = txx_jja.rename({'tasmax':'txxETCCDI'})

# Convert kelvin to celcius
for dataset in [txx_ann, txx_djf, txx_jja]:
    if dataset['txxETCCDI'].max() > 200:
        txx_ann['txxETCCDI'] = txx_ann['txxETCCDI'] - 273.15
        txx_djf['txxETCCDI'] = txx_djf['txxETCCDI'] - 273.15
        txx_jja['txxETCCDI'] = txx_jja['txxETCCDI'] - 273.15
    else:
        None

# Calculation for grid cell area. 
# This is a probably unecessary alternative to taking the cosine of latitude. Calculates area at each cell and then divdes by total area
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
array_area = area_grid(txx_ann['lat'], txx_ann['lon'])
# Total area
total_area = array_area.sum(['lat','lon'])

# ANN
# Area-weighted temp
temp_weighted_ann = (txx_ann['txxETCCDI']*array_area)/total_area
# weighting
weight_tot_ann = (array_area/total_area).sum(('lat','lon'))
# global average. If you are only interested in the global mean temp, the calculation is simple. Just sum all weighted values to get the average.
global_mean_ann = temp_weighted_ann.sum(dim=('lat','lon'))

# DJF--can use ds_ann weights bc same grid
# Area-weighted temp
temp_weighted_djf = (txx_djf['txxETCCDI']*array_area)/total_area
# global average. If you are only interested in the global mean temp, the calculation is simple. Just sum all weighted values to get the average.
global_mean_djf = temp_weighted_djf.sum(dim=('lat','lon'))

# JJA
# if you want to add the global means to the original dataset
# Area-weighted temp
temp_weighted_jja = (txx_jja['txxETCCDI']*array_area)/total_area
# global average. If you are only interested in the global mean temp, the calculation is simple. Just sum all weighted values to get the average.
global_mean_jja = temp_weighted_jja.sum(dim=('lat','lon'))


# Add global mean to dataset
ds_out_ann = txx_ann.assign({'txx_gm_ann':global_mean_ann})
ds_out_djf = txx_djf.assign({'txx_gm_djf': global_mean_djf})
ds_out_jja = txx_jja.assign({'txx_gm_jja': global_mean_jja})

# The first year in DJF will be zero (since no December), so cut off first year
ds_out_djf = ds_out_djf.sel(time=ds_out_djf.time[1:])

# Timestamp for metadata
creation_date = str('Created by NASA GISS ' + 'on ' +  str(datetime.datetime.now()))

# Add in metadata
for dataset in ["ds_out_ann", "ds_out_djf", "ds_out_jja"]:
    locals()[dataset].lat.attrs = {'bounds': 'lat_bnds',
    'units': 'degrees_north',
    'axis': 'Y',
    'long_name': 'latitude',
    'standard_name': 'latitude'}
    locals()[dataset].lon.attrs = {'bounds': 'lon_bnds',
    'units': 'degrees_east',
    'axis': 'X',
    'long_name': 'longitude',
    'standard_name': 'longitude'}

ds_out_ann.txxETCCDI.attrs = {'units': 'degrees_C',
 'long_name': 'Annual Maximum of Daily Maximum Temperature',
 'comment': 'includes both liquid and solid phases',
 'history': creation_date}
ds_out_djf.txxETCCDI.attrs = {'units': 'degrees_C',
 'long_name': 'Seasonal (DJF) Maximum of Daily Maximum Temperature',
 'comment': 'includes both liquid and solid phases',
 'history': creation_date}
ds_out_jja.txxETCCDI.attrs = {'units': 'degrees_C',
 'long_name': 'Seasonal (JJA) Maximum of Daily Maximum Temperature',
 'comment': 'includes both liquid and solid phases',
 'history': creation_date}

ds_out_ann.txx_gm_ann.attrs = {'units': 'degrees_C',
'long_name': 'Global Mean of Annual Maximum of Daily Maximum Temperature',
'history': creation_date}
ds_out_djf.txx_gm_djf.attrs = {'units': 'degrees_C',
'long_name': 'Global Mean of Seasonal (DJF) Maximum of Daily Maximum Temperature',
'history': creation_date}
ds_out_jja.txx_gm_jja.attrs = {'units': 'degrees_C',
'long_name': 'Global Mean of Seasonal (JJA) Maximum of Daily Maximum Temperature',
'history': creation_date}


ds_out_ann.txxETCCDI.encoding = {'unlimited_dims': {'time'},
 'dtype': 'float32'}
ds_out_ann.txx_gm_ann.encoding = {'unlimited_dims': {'time'},
 'dtype': 'float32'}
ds_out_djf.txxETCCDI.encoding = {'unlimited_dims': {'time'},
 'dtype': 'float32'}
ds_out_djf.txx_gm_djf.encoding = {'unlimited_dims': {'time'},
 'dtype': 'float32'}
ds_out_jja.txxETCCDI.encoding = {'unlimited_dims': {'time'},
 'dtype': 'float32'}
ds_out_jja.txx_gm_jja.encoding = {'unlimited_dims': {'time'},
 'dtype': 'float32'}

# convert time to datetime
try:
    for dataset in ["ds_out_ann", "ds_out_djf", "ds_out_jja"]:
        locals()[dataset]['time'] = pd.to_datetime(locals()[dataset].time, format='%Y')
except:
    None

# Create new dataset with correct name
model = options.file.split('tasmax_day_')[1].split('_')[0]
scenario = options.file.split('tasmax_day_')[1].split('_')[1].split('_')[0]
variant = options.file.split('tasmax_day_')[1].split('_')[2]
years = options.file.split('tasmax_day_')[1].split('_')[4].split('.nc')[0]

# Make datasets
ds_out_ann.to_netcdf(options.outdir + '/txxETCCDI' + '_ANN_' + model + '_' + scenario + '_' + variant + '_' + years + '.nc')
ds_out_djf.to_netcdf(options.outdir + '/txxETCCDI' + '_DJF_' + model + '_' + scenario + '_' + variant + '_' + years + '.nc')
ds_out_jja.to_netcdf(options.outdir + '/txxETCCDI' + '_JJA_' + model + '_' + scenario + '_' + variant + '_' + years + '.nc')

# Pause stopwatch
stop = timeit.default_timer()
# Get runtime
print('Runtime: ' + str(stop - start) + ' seconds')
