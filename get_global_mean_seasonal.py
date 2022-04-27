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
import random

# Start stopwatch
start = timeit.default_timer()

def readOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="The parsing commands lists.")
    parser.add_argument("-f", "--file", help="Type your file (must be monthly)")
    parser.add_argument("-v", "--variable", help="Type your variable name")
    parser.add_argument("-o", "--outdir", help="Type your output directory here")
    opts = parser.parse_args(args)
    return opts

# Call the function to read the argument values
options = readOptions(sys.argv[1:])

# Get variable for output directory string
outdir = options.outdir

# import file(s)
ds = xr.open_dataset(options.file)
ds_djf = ds.where(ds['time.season']=='DJF', drop=True)
ds_jja = ds.where(ds['time.season']=='JJA', drop=True)

if options.variable == 'txxETCCDI':
    var = 'txxETCCDI'
    # XArray does not have an easy was to resample season, choose 'DJF', and drop others so my workaround is to take the monthly min/max, calculate the rolling 3month min/max (centered), and drop all non-january months 
    djf = ds_djf.rolling(time=3, center=True).max()
    djf = djf.where(djf['time.month']==1, drop=True)

    # The 1-year resampling should work here bc JJA does not span years, but do the same procedure as DJF for consistency
    jja = ds_jja.rolling(time=3, center=True).max()
    jja = jja.where(jja['time.month']==7, drop=True)
elif options.variable == 'tnnETCCDI':
    var = 'tnnETCCDI'
    # XArray does not have an easy was to resample season, choose 'DJF', and drop others so my workaround is to take the monthly min/max, calculate the rolling 3month min/max (centered), and drop all non-january months 
    djf = ds_djf.rolling(time=3, center=True).min()
    djf = djf.where(djf['time.month']==1, drop=True)

    # The 1-year resampling should work here bc JJA does not span years, but do the same procedure as DJF for consistency
    jja = ds_jja.rolling(time=3, center=True).min()
    jja = jja.where(jja['time.month']==7, drop=True)
else:
    print("variable must be 'txxETCCDI' or 'tnnETCCDI'")
    print("if you would like to calculate the seasonal mean for a different variable contact maxwell.t.elling@nasa.gov")

# Drop time_bnds for seasonal
djf = djf.drop_vars('time_bnds')
jja = jja.drop_vars('time_bnds')

try:
  detectlat = djf.lat.mean()
  detectlat = []
except:
  djf = djf.rename({'latitude': 'lat', 'longitude': 'lon'})
  jja = jja.rename({'latitude': 'lat', 'longitude': 'lon'})

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
array_area = area_grid(djf['lat'], djf['lon'])
# Total area
total_area = array_area.sum(['lat','lon'])

# DJF--can use ds_ann weights bc same grid
# Area-weighted temp
temp_weighted_djf = (djf[options.variable]*array_area)/total_area
# global average. If you are only interested in the global mean temp, the calculation is simple. Just sum all weighted values to get the average.
global_mean_djf = temp_weighted_djf.sum(dim=('lat','lon'))

# JJA
# if you want to add the global means to the original dataset
# Area-weighted temp
temp_weighted_jja = (jja[options.variable]*array_area)/total_area
# global average. If you are only interested in the global mean temp, the calculation is simple. Just sum all weighted values to get the average.
global_mean_jja = temp_weighted_jja.sum(dim=('lat','lon'))

# Add global mean to dataset
ds_out_djf = xr.Dataset()
ds_out_djf = ds_out_djf.assign({str(var + '_gm_djf'): global_mean_djf})
ds_out_jja = xr.Dataset()
ds_out_jja = ds_out_jja.assign({str(var + '_gm_jja'): global_mean_jja})

# Add and update attributes to ds
units = ds[options.variable].attrs['units']

# Get specifics for attributes
if var == 'tnnETCCDI':
    extrema = 'Minimum'
elif var == 'txxETCCDI':
    extrema = 'Maximum'

try:
  ds_out_djf[str(var + '_gm_djf')].attrs = {'units': units,
                            'long_name': str('Global Mean of Seasonal (DJF) ' + extrema + ' of Daily ' + extrema + ' Temperature')}
  ds_out_jja[str(var + '_gm_jja')].attrs = {'units': units,
                            'long_name': str('Global Mean of Seasonal (JJA) ' + extrema + ' of Daily ' + extrema + ' Temperature')}
except:
  print('No "units" or "long_name" attributes detected')

# If time is not time dtype,  convert to datetime
try:
  ds_out_djf['time'] = pd.to_datetime(ds_out_djf.time)
  ds_out_jja['time'] = pd.to_datetime(ds_out_jja.time)
except:
  None

# The first year in DJF will be zero (since no December), so cut off first year
ds_out_djf = ds_out_djf.sel(time=ds_out_djf.time[1:])

# Create new file name
if '/' in options.file:
  if options.file[-1] == "/":
    strfile = (options.file)[0:-1]
    try:
        inputfile = (strfile).split('/')[-1]
        if '_mon_' in inputfile:
            outfile_name = str(inputfile.split('_mon_')[0] + '_time_' + inputfile.split('_mon_')[1])
        elif '_yr_' in inputfile:
            outfile_name = str(inputfile.split('_yr_')[0] + '_time_' + inputfile.split('_yr_')[1])
    except:
        outfile_name = str(options.variable + '_time_' + str(random.randint(0,10000000)) + '.nc')
  else:
    try:
        inputfile = (options.file).split('/')[-1]
        if '_mon_' in inputfile:
            outfile_name = str(inputfile.split('_mon_')[0] + '_time_' + inputfile.split('_mon_')[1])
        elif '_yr_' in inputfile:
            outfile_name = str(inputfile.split('_yr_')[0] + '_time_' + inputfile.split('_yr_')[1])
    except:
        outfile_name = str(options.variable + '_time_' +  str(random.randint(0,10000000)) + '.nc')
else:
  try:
      inputfile = (options.file)
      if '_mon_' in inputfile:
            outfile_name = str(inputfile.split('_mon_')[0] + '_time_' + inputfile.split('_mon_')[1])
      elif '_yr_' in inputfile:
            outfile_name = str(inputfile.split('_yr_')[0] + '_time_' + inputfile.split('_yr_')[1])
  except:
      outfile_name = str(options.variable + '_time_' +  str(random.randint(0,10000000)) + '.nc')

# Create file name for ANN, DJF, JJA
outfile_djf = str(outfile_name.split('_time_')[0] + '_djf_' + outfile_name.split('_time_')[1])
outfile_jja = str(outfile_name.split('_time_')[0] + '_jja_' + outfile_name.split('_time_')[1])

# Create new netcdf
if outdir[-1] == '/':

    out_file_djf = str(str(outdir) + outfile_djf)
    ds_out_djf.to_netcdf(outfile_djf)

    out_file_jja = str(str(outdir) + outfile_jja)
    ds_out_jja.to_netcdf(outfile_jja)

    print(out_file_djf + ' created')
    print(out_file_jja + ' created')
    print('The names of the files may be strange. Please change if needed')
else:
    out_file = str(str(outdir) + '/' + outfile_name)

    out_file_djf = str(str(outdir) + '/' + outfile_djf)
    ds_out_djf.to_netcdf(out_file_djf)

    out_file_jja = str(str(outdir) + '/' + outfile_jja)
    ds_out_jja.to_netcdf(out_file_jja)

    print(out_file_djf + ' created')
    print(out_file_jja + ' created')
    print('The names of the files may be strange. Please change if needed')

# Break up alerts for easier reading
print("...")

# Pause stopwatch
stop = timeit.default_timer()
# Calculate run time
print('Runtime: ' + str((stop - start)) + ' seconds')
