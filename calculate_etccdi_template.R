# load essential packages
# packages are now available in R/4.1.0 module on discover
library(climdex.pcic.ncdf)
library(udunits2)
library(ncdf4)

# USAGE
## To run this script, enter "nohup Rscript ./calculate_etccdi_template.R &"

# list of one to three input files. e.g. c("file_a.nc","file_b.nc","file_c.nc")
# input files (daily) = tasmin, tasmax, precip
# note: you do not have to input all variables (i.e. you can still compute temp indices while excluding precip file in infiles) just note the precip indices will be excluded
infiles=c("/gpfsm/dnb53/projects/p117/pub/CMIP6/ScenarioMIP/NASA-GISS/GISS-E2-1-G/ssp119/r1i1p1f2/day/pr/gn/v20200115/pr_day_GISS-E2-1-G_ssp119_r1i1p1f2_gn_20150101-21001231.nc",
          "/gpfsm/dnb53/projects/p117/pub/CMIP6/ScenarioMIP/NASA-GISS/GISS-E2-1-G/ssp119/r1i1p1f2/day/tasmax/gn/v20200115/tasmax_day_GISS-E2-1-G_ssp119_r1i1p1f2_gn_20150101-21001231.nc",
          "/gpfsm/dnb53/projects/p117/pub/CMIP6/ScenarioMIP/NASA-GISS/GISS-E2-1-G/ssp119/r1i1p1f2/day/tasmin/gn/v20200115/tasmin_day_GISS-E2-1-G_ssp119_r1i1p1f2_gn_20150101-21001231.nc")
#infiles=c("/path/to/precip/file.nc", "/path/to/tasmax/file.nc", "/path/to/tasmin/file.nc")

# list of variable names according to above file(s)
#vars=c(tmax="tasmax", tmin="tasmin")
vars=c(prec="pr", tmax="tasmax", tmin="tasmin")

# output directory. Will be created if it does not exist. Files will appear here when created
# Note: Each index will be exist in its own netcdf file. If you would like all the variables stored in one netcdf file, you will have to do it yourself or reach out to me for a sample script of how to do this (in python) 
outdir="/discover/nobackup/projects/giss_ana/users/melling"

# Output filename format. Must use CMIP5 filename convention. i.e. "var_timeresolution_model_scenario_run_starttime-endtime.nc"
file.template="var_annual_GISS-E2-1-G_scenario_r1i1p1f2_20150101-21001231.nc"
#file.template="var_annual_GISS-E2-1-G_historical_r1i1p1f2_18500101-19491231.nc"
#file.template="var_timeresolution_model_scenario_run_YYYYMMDD-YYYYMMDD.nc" # i.e. "var_annual_GISS-E2-1-G_historical_r1i1p1f1_19500101-20141231.nc"

# author data
author.data=list(institution="GISS", institution_id="id")

# reference period
base.range=c(2020,2060) # reference period for percentile-based indices 

# number of cores to use, or FALSE for single core.
cores=FALSE

# list of indices to calculate, or NULL to calculate all. for the full list of indices see https://www.climdex.org/learn/indices/
#indices=c("txx","tnn")	# alternatively, can specify indices--> i.e. c("txx","tnn", "wsdi")
indices=NULL	# alternatively, can specify indices--> i.e. c("txx","tnn", "wsdi")

# input threshold file to use, or NULL for none.
thresholds.files=NULL#"thresholds.test.1991-1997.nc"

#######################################################
# Esoterics below, do not modify without a good reason.

# definition used for Excess Heat Factor (EHF). "PA13" for Perkins and Alexander (2013), this is the default. "NF13" for Nairn and Fawcett (2013).
EHF_DEF = "PA13"

# axis to split data on. For chunking up of grid, leave this.
axis.name="Y"

# number of data values to process at once. If you receive "Error: rows.per.slice >= 1 is not TRUE", try increasing this to 20. You might have a large grid.
maxvals=10

# output compatible with FCLIMDEX. Leave this.
fclimdex.compatible=FALSE

# time resolution. i.e. if you want just annual values set timeres = "annual". The default ("all") will compute all annual files and monthly where applicable
#timeres = "annual" #can chose to do just annual or just monthly--> i.e. "monthly" 
timeres = "all" #can chose to do just annual or just monthly--> i.e. "monthly" 

# call the package. This will spit out individual netcdf files for each index
# it may take a bit to run, depending on the reference period and time resolution. If computing all indices for multiple decades at monthly and annual resolutions, it could take hours to complete
# you can check if there are files in your output directory to ensure it is working. It will continuously populate these files as it runs
create.indices.from.files(infiles,outdir,file.template,author.data,variable.name.map=vars,base.range=base.range,parallel=cores,axis.to.split.on=axis.name,
                          climdex.vars.subset=indices,thresholds.files=thresholds.files,fclimdex.compatible=fclimdex.compatible,
                          cluster.type="SOCK",max.vals.millions=maxvals,climdex.time.resolution = timeres,
                          thresholds.name.map=c(tx05thresh="tx05thresh",tx10thresh="tx10thresh", tx50thresh="tx50thresh", tx90thresh="tx90thresh",tx95thresh="tx95thresh", 
                                                tn05thresh="tn05thresh",tn10thresh="tn10thresh",tn50thresh="tn50thresh",tn90thresh="tn90thresh",tn95thresh="tn95thresh",
                                                tx90thresh_15days="tx90thresh_15days",tn90thresh_15days="tn90thresh_15days",tavg90thresh_15days="tavg90thresh_15days",
                                                tavg05thresh="tavg05thresh",tavg95thresh="tavg95thresh",
                                                txraw="txraw",tnraw="tnraw",precraw="precraw", 
                                                r95thresh="r95thresh", r99thresh="r99thresh"))


# NOTE: if you are getting "error in nc.create!" it is likely that you already have files with the same names in your output directory. 
## this is possible if you previously ran the code even if it was unsuccessful. Delete these files in the output directory and run again.
