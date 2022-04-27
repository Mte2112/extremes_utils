import xarray as xr
from matplotlib import pyplot as plt
from datetime import datetime as dt
import numpy as np
import pandas as pd
import timeit
import sys
import datetime
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from datetime import datetime, timedelta

# IMPORT THE SOURCE FILES
sst = xr.open_dataset('1980-2014.aijh6E200f10aF40oQ40_subddjr_sorted_repl-na.nc')
#tsavg = xr.open_dataset('E200f10aF40oQ40_subddjr_sorted.nc')
sst_thr = xr.open_dataset('sst_thresh.nc')
#tsavg_thr = xr.open_dataset('./thresholds_files/tsavg_thresh.nc')

# Getsst_thresh diags
sst_thr = sst_thr.assign({'sst':sst_thr.tx90thresh}).drop('tx10thresh')

# Convert time coordinate from cftime to datetime (datetime is easier to apply functions to)
datetimeindex = sst.indexes['time'].to_datetimeindex()
sst['time'] = datetimeindex
#tsavg['time'] = datetimeindex

# open sea land mask
mask = xr.open_dataset('Z2HX2fromZ1QX1N.BS1.nc')

sst_thr['time'] = np.arange(1,366,1)
#tsavg_thr['time'] = np.arange(1,366,1)

sst_thr = sst_thr.rename({'time':'dayofyear'})
#tsavg_thr = tsavg_thr.rename({'time':'dayofyear'})

# Mask land--> Needs more logic to be able to do this without problems over land
sst = sst.where(mask.focean, drop=True)
#tsavg = tsavg.where(mask.focean, drop=True)
sst_thr = sst_thr.where(mask.focean, drop=True)
#tsavg_thr = tsavg_thr.where(mask.focean, drop=True)
###############################################################################################################################################
# Create master and save it as netcdf
master = xr.Dataset()
master.to_netcdf('sst_master_hw.nc')

latlb = -39
for n in np.arange(-37, 89, 2):
    latlb = latlb + 2 
    latub = latlb + 2
    lonlb = (-178.75 - 2.5)
    for p in sst.lon:
        lonlb = lonlb + 2.5
        lonub = lonlb + 2.5
        
        # IMPORT THE SOURCE FILES
        sst_onept = sst.sel(lat=slice(latlb, latub), lon=slice(lonlb, lonub))
        #tsavg_onept = tsavg.sel(lat=slice(latlb, latub), lon=slice(lonlb, lonub))
        sst_thresh = sst_thr.sel(lat=slice(latlb, latub), lon=slice(lonlb, lonub))
        #tsavg_thresh = tsavg_thr.sel(lat=slice(latlb, latub), lon=slice(lonlb, lonub))
        
        # Get start and end year
        startyear = 1981 # First year in data
        endyear = 2014 # Last year in data
        #window #currently only using 5-day window, can update to change if needed (probably should)

        num_coord_pairs = len(sst_onept.lat)*len(sst_onept.lon)
        num_years = len(np.arange(startyear,(endyear + 1),1))
        num_365day_iter = num_coord_pairs*num_years
        dayofyear_array=np.asarray(np.arange(1,366,1).tolist()*num_365day_iter)

        # Import ds to compare to thresh
        dummy_ds = sst_onept

        # convert to pandas df, add dayofyear array, set index and convert back to xarray
        dummy_df = dummy_ds.to_dataframe().reset_index()
        dummy_df = dummy_df.sort_values(by=['lat','lon'], axis=0)
        dummy_df['dayofyear'] = dayofyear_array
        dummy_df = dummy_df.drop('axyp', axis=1)
        dummy_df = dummy_df.set_index(['time','dayofyear','lat','lon'])

        # note: may not need to convert back to xarray bc the code below operates in df
        # converting back and forth formats it correctly
        #dummy_ds3 = xr.Dataset().from_dataframe(dummy_df, sparse = True) # sparse uses less memory and is good for multiindex
        dummy_ds3 = dummy_df.to_xarray()
        # remove dummy_df for memory
        #dummy_df = []
        #dummy_df2 = dummy_ds2.to_dataframe()

        del dayofyear_array
        del dummy_ds 
        del dummy_df  
        del num_coord_pairs 
        del num_years 
        del num_365day_iter 

        # Compute difference between the temps and threshold temps. Positive values will be above threshold
        dummy_ds3_vs_thresh = dummy_ds3.groupby('dayofyear') - sst_thresh
        try:
            dummy_ds3_vs_thresh = dummy_ds3_vs_thresh.to_dataframe().reset_index().drop('dayofyear', axis=1).dropna().set_index(['time', 'lat','lon']).to_xarray()
            #dummy_ds3_vs_thresh_timeseries = dummy_ds3_vs_thresh
            #dummy_ds3_vs_thresh = []
            dummy_ds3_vs_thresh = dummy_ds3_vs_thresh.sortby('time')

            def get_annual_wsdi(x):
                # Calculate heat waves
                consec_hotdays = []
                heatwave_lengths = []
                heatwave_enddate = []
                bool_array = []
                for diff in x:
                    if diff > 0:
                        consec_hotdays.append(1)
                    elif diff <= 0: # the 'equal to' is important bc the above does not include values == the 90p threshold. If not here, values == threshold will be skipped 
                        sum_hotdays = sum(consec_hotdays)
                        if sum_hotdays >= 6:
                            heatwave_lengths.append(sum_hotdays)
                            heatwave_enddate.append((diff.time - np.timedelta64(1, "D")).values) # End date of heat wave. to get time series for the heatwave lengths (to average by year, etc)
                            bool_array.extend([True]*sum_hotdays)
                            bool_array.extend([False]) # This is the count for the current non-hot day. otherwise, only non-heatwave hot days would be counted as False
                            consec_hotdays = [] 
                        else:
                            bool_array.extend([False]*sum_hotdays)
                            bool_array.extend([False]) # This is the count for the current non-hot day. otherwise, only non-heatwave hot days would be counted as False
                            consec_hotdays = []
                    else:
                        #print('Hmmm... this value may be a nan or other non-float value. Please assess below')
                        #print(diff)
                        #sys.exit()
                        None
                # Now check to see if there are any heatwave days at the end of the time series that were not added due to the loop ending on a heatwave day
                # Append these to the array
                if len(x) != len(bool_array):
                    num_tailend_hotdays = len(x) - len(bool_array)
                    if num_tailend_hotdays >= 6:
                        bool_array.extend([True]*num_tailend_hotdays)
                        heatwave_lengths.append(num_tailend_hotdays)
                        heatwave_enddate.append((diff.time - np.timedelta64(1, "D")).values)
                    else:
                        bool_array.extend([False]*num_tailend_hotdays)
                else:
                    None

                # This is not as accurate as my handmade version above. For now, just use this as a placeholder for the timeseries dataset format and feed my boolean array into it
                condition = (x > 0) & \
                            (x.shift(time=-1) > 0) & \
                            (x.shift(time=-2) > 0) & \
                            (x.shift(time=-3) > 0) & \
                            (x.shift(time=-4) > 0) & \
                            (x.shift(time=-5) > 0) 

                consecutive_count = x[condition].count()

                heatwave_days = x.where(condition)
                heatwave_days_nonan = x.where(condition, drop=True)

                # Mask dataset by boolean array created above
                bool_array = np.asarray(bool_array)

                condition_self_computed = condition
                condition_self_computed.values = bool_array

                heatwave_days2 = x.where(condition_self_computed)

                # add diagnostic to ds maybe add more later
                return xr.DataArray(heatwave_days2.resample(time='1Y').count())

            # Calculate thresholds by stack the coordinates and applying the get_thresh() function
            # Note- as mentioned above, the function only returns a value when the day is in the centered of the window (thus returning only one percentile value per day of year)`
            stacked = dummy_ds3_vs_thresh.sst.stack(allpoints=['lat', 'lon'])
            # apply the function over allpoints to calculate the trend at each point
            wsdi = stacked.groupby('allpoints').apply(get_annual_wsdi)
            # unstack back to lat lon coordinates
            wsdi_unstacked = wsdi.unstack('allpoints')

            def get_average_annual_length(x):
                # Calculate heat waves
                consec_hotdays = []
                heatwave_lengths = []
                heatwave_enddate = []
                bool_array = []
                for diff in x:
                    if diff > 0:
                        consec_hotdays.append(1)
                    elif diff <= 0: # the 'equal to' is important bc the above does not include values == the 90p threshold. If not here, values == threshold will be skipped 
                        sum_hotdays = sum(consec_hotdays)
                        if sum_hotdays >= 6:
                            heatwave_lengths.append(sum_hotdays)
                            heatwave_enddate.append((diff.time - np.timedelta64(1, "D")).values) # End date of heat wave. to get time series for the heatwave lengths (to average by year, etc)
                            bool_array.extend([True]*sum_hotdays)
                            bool_array.extend([False]) # This is the count for the current non-hot day. otherwise, only non-heatwave hot days would be counted as False
                            consec_hotdays = [] 
                        else:
                            bool_array.extend([False]*sum_hotdays)
                            bool_array.extend([False]) # This is the count for the current non-hot day. otherwise, only non-heatwave hot days would be counted as False
                            consec_hotdays = []
                    else:
                        #print('Hmmm... this value may be a nan or other non-float value. Please assess below')
                        #print(diff)
                        #sys.exit()
                        None

                # Now check to see if there are any heatwave days at the end of the time series that were not added due to the loop ending on a heatwave day
                # Append these to the array
                if len(x) != len(bool_array):
                    num_tailend_hotdays = len(x) - len(bool_array)
                    if num_tailend_hotdays >= 6:
                        bool_array.extend([True]*num_tailend_hotdays)
                        heatwave_lengths.append(num_tailend_hotdays)
                        heatwave_enddate.append((diff.time - np.timedelta64(1, "D")).values)
                    else:
                        bool_array.extend([False]*num_tailend_hotdays)
                else:
                    None

                # This is not as accurate as my handmade version above. For now, just use this as a placeholder for the timeseries dataset format and feed my boolean array into it
                condition = (x > 0) & \
                            (x.shift(time=-1) > 0) & \
                            (x.shift(time=-2) > 0) & \
                            (x.shift(time=-3) > 0) & \
                            (x.shift(time=-4) > 0) & \
                            (x.shift(time=-5) > 0) 

                consecutive_count = x[condition].count()

                heatwave_days = x.where(condition)
                heatwave_days_nonan = x.where(condition, drop=True)

                # Mask dataset by boolean array created above
                bool_array = np.asarray(bool_array)

                condition_self_computed = condition
                condition_self_computed.values = bool_array

                heatwave_days2 = x.where(condition_self_computed)

                # create dataframe for heatwave length and end date
                hl_df = pd.DataFrame({'length': heatwave_lengths, 'enddate': heatwave_enddate}).set_index('enddate')
                # convert to xarray and resample by year, take mean
                hl_ds = hl_df.to_xarray()
                annual_avg_hl = hl_ds.resample(enddate='1Y').mean()
                annual_count_heatwaves = hl_ds.resample(enddate='1Y').count()

                return xr.DataArray(annual_avg_hl.length)

            # Calculate thresholds by stack the coordinates and applying the get_thresh() function
            # Note- as mentioned above, the function only returns a value when the day is in the centered of the window (thus returning only one percentile value per day of year)`
            stacked2 = dummy_ds3_vs_thresh.sst.stack(allpoints=['lat', 'lon'])
            # apply the function over allpoints to calculate the trend at each point
            len_hw = stacked2.groupby('allpoints').apply(get_average_annual_length)
            # unstack back to lat lon coordinates
            len_hw_unstacked = len_hw.unstack('allpoints')

            def get_annual_count_hw(x):
                # Calculate heat waves
                consec_hotdays = []
                heatwave_lengths = []
                heatwave_enddate = []
                bool_array = []
                for diff in x:
                    if diff > 0:
                        consec_hotdays.append(1)
                    elif diff <= 0: # the 'equal to' is important bc the above does not include values == the 90p threshold. If not here, values == threshold will be skipped 
                        sum_hotdays = sum(consec_hotdays)
                        if sum_hotdays >= 6:
                            heatwave_lengths.append(sum_hotdays)
                            heatwave_enddate.append((diff.time - np.timedelta64(1, "D")).values) # End date of heat wave. to get time series for the heatwave lengths (to average by year, etc)
                            bool_array.extend([True]*sum_hotdays)
                            bool_array.extend([False]) # This is the count for the current non-hot day. otherwise, only non-heatwave hot days would be counted as False
                            consec_hotdays = [] 
                        else:
                            bool_array.extend([False]*sum_hotdays)
                            bool_array.extend([False]) # This is the count for the current non-hot day. otherwise, only non-heatwave hot days would be counted as False
                            consec_hotdays = []
                    else:
                        #print('Hmmm... this value may be a nan or other non-float value. Please assess below')
                        #print(diff)
                        #sys.exit()
                        None

                # Now check to see if there are any heatwave days at the end of the time series that were not added due to the loop ending on a heatwave day
                # Append these to the array
                if len(x) != len(bool_array):
                    num_tailend_hotdays = len(x) - len(bool_array)
                    if num_tailend_hotdays >= 6:
                        bool_array.extend([True]*num_tailend_hotdays)
                        heatwave_lengths.append(num_tailend_hotdays)
                        heatwave_enddate.append((diff.time - np.timedelta64(1, "D")).values)
                    else:
                        bool_array.extend([False]*num_tailend_hotdays)
                else:
                    None

                # This is not as accurate as my handmade version above. For now, just use this as a placeholder for the timeseries dataset format and feed my boolean array into it
                condition = (x > 0) & \
                            (x.shift(time=-1) > 0) & \
                            (x.shift(time=-2) > 0) & \
                            (x.shift(time=-3) > 0) & \
                            (x.shift(time=-4) > 0) & \
                            (x.shift(time=-5) > 0) 

                consecutive_count = x[condition].count()

                heatwave_days = x.where(condition)
                heatwave_days_nonan = x.where(condition, drop=True)

                # Mask dataset by boolean array created above
                bool_array = np.asarray(bool_array)

                condition_self_computed = condition
                condition_self_computed.values = bool_array

                heatwave_days2 = x.where(condition_self_computed)

                # create dataframe for heatwave length and end date
                hl_df = pd.DataFrame({'length': heatwave_lengths, 'enddate': heatwave_enddate}).set_index('enddate')
                # convert to xarray and resample by year, take mean
                hl_ds = hl_df.to_xarray()
                annual_avg_hl = hl_ds.resample(enddate='1Y').mean()
                annual_count_heatwaves = hl_ds.resample(enddate='1Y').count()

                return xr.DataArray(annual_count_heatwaves.length)

            # Calculate thresholds by stack the coordinates and applying the get_thresh() function
            # Note- as mentioned above, the function only returns a value when the day is in the centered of the window (thus returning only one percentile value per day of year)`
            stacked3 = dummy_ds3_vs_thresh.sst.stack(allpoints=['lat', 'lon'])
            # apply the function over allpoints to calculate the trend at each point
            count_hw = stacked3.groupby('allpoints').apply(get_annual_count_hw)
            # unstack back to lat lon coordinates
            count_hw_unstacked = count_hw.unstack('allpoints')

            # Make master dataset
            heatwave_diagnostics = xr.Dataset()
            heatwave_diagnostics = heatwave_diagnostics.assign({'len':len_hw_unstacked, 'hwct':count_hw_unstacked, 'wsdi':wsdi_unstacked })
            
            # Concat to ds
            xr.concat(objs = [xr.open_dataset('sst_master_hw.nc'), heatwave_diagnostics], dim='concat_dim', data_vars='all', coords='minimal' ).to_netcdf('sst_master_hw.nc')

            # Delete variables 
            del count_hw_unstacked
            del count_hw
            del stacked3
            del stacked2
            del stacked
            del wsdi
            del heatwave_diagnostics


            # Show progress
            print('[' + str(latlb) + ', ' + str(lonlb) + ']' + ' done') 
            print('[' + str(latlb) + ', ' + str(lonub) + ']' + ' done') 
            print('[' + str(latub) + ', ' + str(lonlb) + ']' + ' done') 
            print('[' + str(latub) + ', ' + str(lonub) + ']' + ' done') 
        except:
            print('[' + str(latlb) + ', ' + str(lonlb) + ']' + ' error') 
            print('[' + str(latlb) + ', ' + str(lonub) + ']' + ' error') 
            print('[' + str(latub) + ', ' + str(lonlb) + ']' + ' error') 
            print('[' + str(latub) + ', ' + str(lonub) + ']' + ' error')   
 
