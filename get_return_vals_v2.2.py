# Import packages
from pyextremes import EVA, __version__
import xarray as xr
import pandas as pd
import numpy as np
import timeit
import scipy
import sys
import getopt
import argparse
from datetime import date
import random
import warnings

# Disable warnings
warnings.filterwarnings("ignore")

# Start stopwatch
start = timeit.default_timer()


# Declare function to define command-line arguments
def readOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="The parsing commands lists.")
    parser.add_argument("-f", "--file", help="Type your file")
    parser.add_argument("-p", "--period", help="Type your return period (default is 20)")
    parser.add_argument("-v", "--variable", help="Type your variable name")
    parser.add_argument("-o", "--outdir", help="Type your output directory (default is work dir)")
    parser.add_argument("-a", "--alpha", help="Type your alpha (default is 0.95)")
    parser.add_argument("-n", "--samplesize", type=int, help="Type your sample size (default is 1000)")
    parser.add_argument("-m", "--method", help="Type your method (default is 'BM')")
    parser.add_argument("-b", "--blocksize", help="Type your block size (default is 365.2425D (1Y))")
    parser.add_argument("-t", "--extremestype", help="Type your extremes type (default is 'high')")
    parser.add_argument("-y", "--years", nargs='+', type=int, help="Type the start years where you want to calculate return values type (i.e. '-y 1980 1990 2000') (default is every year from the first year of the dataset to last year - w)")
    parser.add_argument("-w", "--timewindow", type=int, help="Type your time window for GEV estimation **must be an odd number** (default is 31  (30 in Li et al. 2021.) Only change for good reason)")

    opts = parser.parse_args(args)
    return opts


# Call the function to read the argument values
options = readOptions(sys.argv[1:])
# Initialize array of changeable argument names to iterate over
options_iter = [options.period, options.outdir, options.alpha, options.samplesize, 
                options.method, options.blocksize, options.extremestype,
                options.years, options.timewindow]
# Initialize array of parameters as strings to set as variables
str_vars = ['period', 'outdir', 'alpha', 'samplesize', 'method', 
            'blocksize', 'extremestype', 'years', 'timewindow']
# Initialize array of default values
default_vals = [20, './', 0.95, 1000, 'BM', '365.2425D', 'high', 'all', 31]


# Check if minimum command-line arguments passed
if options.file == None:
    print("No file argument passed")
    print("Please use command line argument '-f file'")
    sys.exit()
elif options.variable == None:
    print("No variable argument passed")
    print("Please use command line argument '-v variable'")
    sys.exit()
else:
    None

# Import data, convert to datetime
ds = xr.open_dataset(options.file)
# If cftime index, change to datetime
try:
    datetimeindex = ds.indexes['time'].to_datetimeindex()
    ds['time'] = datetimeindex
except:
    print('Input not CFtime NOleap. Not converted to dateime. May already be in correct format.' )

# Set variable values for model parameters
print("-------- Model Parameters --------")
for option, default, var in zip(options_iter, default_vals, str_vars):
    if option == None:
        locals()[var] = default
        print(var + ' set to default: ' + str(default))
    else:
        locals()[var] = option 
        print(var + ' set to ' + str(option))

##### Defining the function depending on the given (or not given) arguments for 'years' #####
##### Defining the function depending on the given (or not given) arguments for 'years' #####
# Check that input years are valid. If year + window > last year of the dataset, exit
# set time window +- to a variable 'window_wings'
window_wings = int((timewindow-1)/2)

for y in years:
    endyr = y + window_wings
    if endyr > int(ds['time'][-1].dt.strftime("%Y")):
        sys.exit(str(y) + ' is not a valid input for years. ' + str(y) + ' + specified timewindow exceeds the last year of the dataset')

# If no argument passed, will calculate $timewindow year period return values for each year from the first year of the dataset to (last year - timewindow) year
if years == 'all':
    # Define first and last year of the dataset for slicing and N-year intervals
    firstyear = ds['time'][0].dt.strftime("%Y")
    lastyear = ds['time'][-1].dt.strftime("%Y")

    last_gev_startyear = int(lastyear) - window_wings + 1
    num_startyears = int(lastyear)-int(firstyear) - window_wings + 1 # -timewindow +1 because and +1 for the second argument in arange
    
    # Define compute_gev for all years from start to (end-timewindow)
    def compute_gev(x):
        mastersummary=pd.DataFrame(columns=['return period', 'year', 'return value', 'lower ci', 'upper ci']).set_index('return period')
        for n in np.arange(0, num_startyears, 1):
            start = str(int(firstyear) + n)
            endyear = str(int(firstyear) + n + window_wings) 

            # try except is an alternative to setting errors to 'ignore'
            ser = x.sel(time=slice(start, endyear, 1)).to_series()
            model = EVA(data=ser)

            model.get_extremes(
            method='BM',
            extremes_type='high',
            block_size='365.2425D',
            errors="coerce",) # errors{'ignore', 'raise', 'coerce'}, default 'raise' If 'raise', then invalid parsing will raise an exception. If 'coerce', then invalid parsing will be set as NaN. If 'ignore', then invalid parsing will return the input.

            model.fit_model()
            summary = model.get_summary(
            return_period=20,
            alpha=0.95,
            n_samples=samplesize,
            )

            # Add year tag to summary
            summary['year'] = start
            mastersummary = mastersummary.append(summary)
        
        # For now, instructing users to only use one return period (capability to use an array (i.e. 20 yr, 50yr, 100yr))
        # If in the future someone wants to use multiple periods, then need to update script
        # For now, can drop 'return period' column in df because having a multi index of year and return period when returned as DataArray does NOT work correctly
        mastersummary = mastersummary.reset_index()
        mastersummary = mastersummary.drop('return period', axis=1)
        mastersummary = mastersummary.set_index(['year'])
        return xr.DataArray(mastersummary)
else:
    def compute_gev(x):
        mastersummary=pd.DataFrame(columns=['return period', 'year', 'return value', 'lower ci', 'upper ci']).set_index('return period')
        for y in years:
            start = str(y)
            endyear = str(int(y) + window_wings) 

            # try except is an alternative to setting errors to 'ignore'
            ser = x.sel(time=slice(start, endyear, 1)).to_series()
            model = EVA(data=ser)

            model.get_extremes(
            method='BM',
            extremes_type='high',
            block_size='365.2425D',
            errors="coerce",) # errors{'ignore', 'raise', 'coerce'}, default 'raise' If 'raise', then invalid parsing will raise an exception. If 'coerce', then invalid parsing will be set as NaN. If 'ignore', then invalid parsing will return the input.

            model.fit_model()
            summary = model.get_summary(
            return_period=20,
            alpha=0.95,
            n_samples=samplesize,
            )

            # Add year tag to summary
            summary['year'] = start
            mastersummary = mastersummary.append(summary)

        # For now, instructing users to only use one return period (capability to use an array (i.e. 20 yr, 50yr, 100yr))
        # If in the future someone wants to use multiple periods, then need to update script
        # For now, can drop 'return period' column in df because having a multi index of year and return period when returned as DataArray does NOT work correctly
        mastersummary = mastersummary.reset_index()
        mastersummary = mastersummary.drop('return period', axis=1)
        mastersummary = mastersummary.set_index(['year'])
        return xr.DataArray(mastersummary)

# Break up alerts for easier reading
print("...")
        
# Define gev/return value function
print("Calculating return values. This may take a while depending on file size and arguments")
# Break up alerts for easier reading
print("...")

# Apply GEV calculation to all grid points
stacked = ds[options.variable].stack(allpoints=['lat','lon'])
# apply the function over allpoints to calculate the trend at each point
returnval = stacked.groupby('allpoints').apply(compute_gev)
# unstack back to lat lon coordinates
returnval_unstacked = returnval.unstack('allpoints')

# Pull data from the dim_01 coordinate (return val, lower ci, upper ci) and create variables for each
out_ds = xr.Dataset()
out_ds = out_ds.assign({'returnval': returnval_unstacked.sel(dim_1='return value', drop=True),
'uci': returnval_unstacked.sel(dim_1='upper ci', drop=True),
'lci': returnval_unstacked.sel(dim_1='lower ci', drop=True)})


# Convert 'year' to datetime and rename to time for combatibility
out_ds['year'] = pd.to_datetime(out_ds.year, format='%Y-%m-%d')

out_ds = out_ds.rename({'year':'time'})

# Add metadata
out_ds.encoding = {'unlimited_dims': {'time'}}
out_ds.lon.attrs = {'bounds': 'lon_bnds',
                    'units': 'degrees_east',
                    'axis': 'X',
                    'long_name': 'longitude',
                    'standard_name': 'longitude'}
out_ds.lat.attrs = {'bounds': 'lat_bnds', 
                    'units': 'degrees_north',
                    'axis': 'Y',
                    'long_name': 'latitude',
                    'standard_name': 'latitude'}
out_ds.attrs = {'institution': 'Goddard Institute for Space Studies, New York, NY 10025, USA',
                'institution_id': 'NASA-GISS',
                'use': 'for internal use only',
                'creation_date': str(date.today()),
                'EVA_software': 'pyextremes, https://github.com/georgebv/pyextremes',
                'EVA_software_version': 'pyextremes v2.2.4',
                'contact': 'Max Elling: maxwell.t.elling@nasa.gov'}
returnvalattr = str(str(period) + ' year return value')
uciattr1 = str('upper confidence interval, alpha = ' + str(alpha))
lciattr1 = str('lower confidence interval, alpha = ' + str(alpha))

try:
  unitatt = ds[options.file].attrs['units']

  out_ds.returnval.attrs = {'units': unitatt,
                            'description': returnvalattr,
                            'long_name': 'return value'}
  out_ds.uci.attrs = {'units':unitatt,
                            'description': uciattr1,
                            'long_name': 'upper_confidence_interval'}

  out_ds.lci.attrs = {'units': unitatt,
                            'description': lciattr1,
                            'long_name': 'lower_confidence_interval'}
except:
  out_ds.returnval.attrs = {'description': returnvalattr,
                            'long_name': 'return value'}
  out_ds.uci.attrs = {'description': uciattr1,
                            'long_name': 'upper_confidence_interval'}

  out_ds.lci.attrs = {'description': lciattr1,
                            'long_name': 'lower_confidence_interval'}
  

# Cummbersome code to create new file name. The program essentially parses out the input file to create and output filename. If not in climdex.pcic format, returns filename in format 'Nyreturn_variableRANDOMNUMBER.nc'  
# The alternative to this is adding another command line argument to specify what you want the file to be called. Maybe implement that next if needed

if '/' in options.file:
  if options.file[-1] == "/":
    strfile = (options.file)[0:-1]
    try:
        inputfile = (options.file).split('/')[-1]
        removethis = str(inputfile.split('_')[0] + '_')

        # Remove the time resolution from the string
        if 'yr_' in inputfile.split(removethis)[1]:
            laststr = inputfile.split(removethis)[1].split('yr_')[0] + inputfile.split(removethis)[1].split('yr_')[1]
        elif 'mon_' in inputfile.split(removethis)[1]:
            laststr = inputfile.split(removethis)[1].split('mon_')[0] + inputfile.split(removethis)[1].split('mon_')[1] 
        else:
            laststr = inputfile.split(removethis)[1]

        outfile_name = str(str(options.variable).split('ETCCDI')[0] + '_' + str(period) + 'yreturn_' + laststr)
    except:
        outfile_name = str(str(period) + 'yreturn_' + options.variable + str(random.randint(0,10000000)) + '.nc')
  else:
    try:
        inputfile = (options.file).split('/')[-1]
        removethis = str(inputfile.split('_')[0] + '_')

        # Remove the time resolution from the string
        if 'yr_' in inputfile.split(removethis)[1]:
            laststr = inputfile.split(removethis)[1].split('yr_')[0] + inputfile.split(removethis)[1].split('yr_')[1]
        elif 'mon_' in inputfile.split(removethis)[1]:
            laststr = inputfile.split(removethis)[1].split('mon_')[0] + inputfile.split(removethis)[1].split('mon_')[1] 
        else:
            laststr = inputfile.split(removethis)[1]

        outfile_name = str(str(options.variable).split('ETCCDI')[0] + '_' + str(period) + 'yreturn_' + laststr)
    except:
        outfile_name = str(str(period) + 'yreturn_' + options.variable + str(random.randint(0,10000000)) + '.nc')
else:
  try:
      inputfile = (options.file)
      removethis = str(inputfile.split('_')[0] + '_')

      # Remove the time resolution from the string
      if 'yr_' in inputfile.split(removethis)[1]:
          laststr = inputfile.split(removethis)[1].split('yr_')[0] + inputfile.split(removethis)[1].split('yr_')[1]
      elif 'mon_' in inputfile.split(removethis)[1]:
          laststr = inputfile.split(removethis)[1].split('mon_')[0] + inputfile.split(removethis)[1].split('mon_')[1] 
      else:
          laststr = inputfile.split(removethis)[1]

      outfile_name = str(str(options.variable).split('ETCCDI')[0] + '_' + str(period) + 'yreturn_' + laststr)
  except:
      outfile_name = str(str(period) + 'yreturn_' + options.variable + str(random.randint(0,10000000)) + '.nc')

# Create new netcdf
if outdir[-1] == '/':
    out_file = str(str(outdir) + outfile_name)
    out_ds.to_netcdf(out_file)
    print('Finished calculating return values')
    print(out_file + ' created')
    print('The name of the file may be off. Please change if needed')
else:
    out_file = str(str(outdir) + '/' + outfile_name)
    out_ds.to_netcdf(out_file)
    print('Finished calculating return values')
    print(out_file + ' created')
    print('The name of the file may be off. Please change if needed')
# Break up alerts for easier reading
print("...")

# Pause stopwatch
stop = timeit.default_timer()
# Calculate run time
print('Runtime: ' + str(stop - start) + ' seconds')

