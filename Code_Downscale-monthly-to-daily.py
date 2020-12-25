#!/usr/bin/env ipython
import xarray as xr
import numpy as np

start_year_temp = 2001 # For BESS dataset
year = input('year: ')

era5 = xr.open_dataset('era5_'+str(year)+'_evap_daily_ydaysum_global_0.25.nc')
era5.e.values = era5.e.values*(-1000)
era5 = era5.sel(latitude = slice(60,-59.99))
era5 = era5.e
era5 = era5.where(era5.values > 0) 
#print(era5.time[0], era5.time[-1])
# Monthly era5 data  
era5_monthly = era5.resample(time = '1MS').sum('time')
era5 = np.array(era5)[:,::-1,:]
era5_monthly = np.array(era5_monthly)[:,::-1,:]

## BESS data is in mm/day #180x360x720 #for 2001-2015
BESS = xr.open_dataset('BESS_Evap_NN.nc')
BESS = BESS.Evaporation.sel(lat = slice(-60,60), time = year)
BESS = BESS.where(BESS.values > 0)
#print(BESS.time[0], BESS.time[-1])
BESS = np.array(BESS)            #added

from time import sleep
from tqdm import tqdm
import numpy as np

if int(year)%4 == 0:
    print('leap')
    dayperyear = 366
    day_per_month = [31,29,31,30,31,30,31,31,30,31,30,31] 
else:
    print('not leap')
    dayperyear = 365
    day_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]

Upscaled_data = np.zeros((dayperyear,480,1440))

## The data processing
for lat in tqdm(range(440)):
    for lon in (range(1440)):
        if int(year)%4 == 0:
            BESS_monthly = BESS[:,lat,lon]*day_per_month
            era5_monthly_resampled = np.concatenate([(np.linspace(era5_monthly[i,lat,lon], era5_monthly[i,lat,lon], day_per_month[i])) for i in range(12)], axis = 0)
            BESS_monthly_resampled = np.concatenate([(np.linspace(BESS_monthly[i], BESS_monthly[i], day_per_month[i])) for i in range(12)], axis = 0)
            Upscaled_data[:,lat,lon] = BESS_monthly_resampled*(era5[:,lat,lon]/era5_monthly_resampled)
        else:
            BESS_monthly = BESS[:,lat,lon]*day_per_month
            era5_monthly_resampled = np.concatenate([(np.linspace(era5_monthly[i,lat,lon], era5_monthly[i,lat,lon], day_per_month[i])) for i in range(12)], axis = 0)
            BESS_monthly_resampled = np.concatenate([(np.linspace(BESS_monthly[i], BESS_monthly[i], day_per_month[i])) for i in range(12)], axis = 0)
            Upscaled_data[:,lat,lon] = BESS_monthly_resampled*(era5[:,lat,lon]/era5_monthly_resampled)
			
import datetime
from netCDF4 import Dataset,num2date,date2num
import numpy as np
# -------------------------------------------------
# Make date axis:
yystart = int(year)
nyears = 1

ntime = nyears*dayperyear
unout = 'days since '+str(year)+'-01-01 00:00:00'
# -------------------------------------------------
res = 0.25
lon = np.arange(0,360,res);nx=np.size(lon)
lat = np.arange(-59.875,60,res);ny=np.size(lat)

dataout = Upscaled_data
datesout = [datetime.datetime(yy,mm,dd,0) for yy in range(yystart,yystart+nyears) for mm in range(1,13) for dd in range(1,day_per_month[mm-1]+1)]
# =================================================
ncout = Dataset('BESS_Upscaled_ET_'+str(year)+'.nc','w','NETCDF4'); # using netCDF3 for output format
ncout.createDimension('lon',nx);
ncout.createDimension('lat',ny);
ncout.createDimension('time',ntime);
lonvar = ncout.createVariable('lon','float32',('lon'));lonvar[:] = lon;
latvar = ncout.createVariable('lat','float32',('lat'));latvar[:] = lat;
timevar = ncout.createVariable('time','float64',('time'));timevar.setncattr('units',unout);timevar[:]=date2num(datesout,unout);
myvar = ncout.createVariable('Evaporation','float32',('time','lat','lon'));myvar.setncattr('units','mm/day');myvar[:] = dataout;
ncout.close();
