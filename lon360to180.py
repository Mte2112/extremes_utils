import pandas as pd
import xarray as xr
import timeit
import getopt
import argparse
import sys

# Start stopwatch
start = timeit.default_timer()

def readOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="The parsing commands lists.")
    parser.add_argument("-f", "--file", help="Type your file (input should be monthly)")
    parser.add_argument("-o", "--outfile", help="Type your output file, include the output directory")
    opts = parser.parse_args(args)
    return opts

# Call the function to read the argument values
options = readOptions(sys.argv[1:])

# Import dataset
ds = xr.open_dataset(options.file)

# Adjust lon values to make sure they are (-180, 180)
# Allows for easier plotting and compatibility with certain python packages
try:
    lon = 'lon'
    ds['_longitude_adjusted'] = xr.where(
        ds[lon] > 180,
        ds[lon] - 360,
        ds[lon])

    # adjusted coords to main lon coords
    # sort by adjusted coords
    ds = (
        ds
        .swap_dims({lon: '_longitude_adjusted'})
        .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
        .drop(lon))

    ds = ds.rename({'_longitude_adjusted': lon})
except:
    try:
        lon = 'longitude'
        ds['_longitude_adjusted'] = xr.where(
            ds[lon] > 180,
            ds[lon] - 360,
            ds[lon])

        # adjusted coords to main lon coords
        # sort by adjusted coords
        ds = (
            ds
            .swap_dims({lon: '_longitude_adjusted'})
            .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
            .drop(lon))

        ds = ds.rename({'_longitude_adjusted': lon})
    except:
        print("Error: 'lon' or 'longitude' dimension not found in dataset")

# Write the netCDF
ds.to_netcdf(options.outfile)

# Pause stopwatch
stop = timeit.default_timer()

# print ds
print(str(options.outfile + ' created'))
# Get runtime
print('Runtime: ' + str(stop - start) + ' seconds')
