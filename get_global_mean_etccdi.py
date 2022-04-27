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
    parser.add_argument("-f", "--file", help="Type your file")
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
array_area = area_grid(ds['lat'], ds['lon'])
# Total area
total_area = array_area.sum(['lat','lon'])
# Area-weighted temp
temp_weighted = (ds[options.variable]*array_area)/total_area
# weighting
weight_tot = (array_area/total_area).sum(('lat','lon'))
# global average. If you are only interested in the global mean temp, the calculation is simple. Just sum all weighted values to get the average.
global_mean = temp_weighted.sum(dim=('lat','lon'))


# Add global mean to dataset
dsgm = xr.Dataset()
dsgm = dsgm.assign({'gm':global_mean})

# Add and update attributes to ds
long_name = ds[options.variable].attrs['long_name']
units = ds[options.variable].attrs['units']

try:
  dsgm['gm'].attrs = {'units': units,
                            'long_name': str(long_name + ' (Global Mean)')}
except:
  print('No "units" or "long_name" attributes detected')

# If time is not time dtype,  convert to datetime
try:
  dsgm['time'] = pd.to_datetime(dsgm.time)
except:
  None

# Create new file name
if '/' in options.file:
  if options.file[-1] == "/":
    strfile = (options.file)[0:-1]
    try:
        inputfile = (options.file).split('/')[-1]
        removethis = str(inputfile.split('_')[0] + '_')
        outfile_name = str('GM_' + str(options.variable).split('ETCCDI')[0] + inputfile.split(removethis)[1])
    except:
        outfile_name = str('GM_' + options.variable + str(random.randint(0,10000000)) + '.nc')
  else:
    try:
        inputfile = (options.file).split('/')[-1]
        removethis = str(inputfile.split('_')[0] + '_')
        outfile_name = str('GM_' + str(options.variable).split('ETCCDI')[0] + inputfile.split(removethis)[1])
    except:
        outfile_name = str('GM_' + options.variable + str(random.randint(0,10000000)) + '.nc')
else:
  try:
      inputfile = (options.file)
      removethis = str(inputfile.split('_')[0] + '_')
      outfile_name = str('GM_' + str(options.variable).split('ETCCDI')[0] + inputfile.split(removethis)[1])
  except:
      outfile_name = str('GM_' + options.variable + str(random.randint(0,10000000)) + '.nc')


# Create new netcdf
if outdir[-1] == '/':
    out_file = str(str(outdir) + outfile_name)
    dsgm.to_netcdf(out_file)
    print('Finished calculating global means')
    print(out_file + ' created')
    print('The name of the file may be strange. Please change if needed')
else:
    out_file = str(str(outdir) + '/' + outfile_name)
    dsgm.to_netcdf(out_file)
    print('Finished calculating return values')
    print(out_file + ' created')
    print('The name of the file may be strange. Please change if needed')

# Break up alerts for easier reading
print("...")

# Pause stopwatch
stop = timeit.default_timer()
# Calculate run time
print('Runtime: ' + str((stop - start)) + ' seconds')
                                                                
