import xarray as xr
import numpy as np
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")

## Calculate maximum precipiration month
P_CHIRPS_argmax = xr.open_mfdataset('/chirps-v2.0.*.days_p25.nc').precip.sel(time = slice('2001','2012'))
P_CHIRPS_argmax = P_CHIRPS_argmax.resample(time = '1MS').sum('time')
P_CHIRPS_argmax = (P_CHIRPS_argmax.groupby('time.month').mean('time'))
P_CHIRPS_argmax = P_CHIRPS_argmax.argmax(axis = 0).values

def data_analysis(year,lat_top, lat_bottom, lon_min, lon_max):
    global P_CHIRPS, E_Ensemble
    P_CHIRPS = xr.open_mfdataset('/chirps-v2.0.*.days_p25.nc').precip.sel(time = slice(str(year),str(year+1)))
    P_CHIRPS = P_CHIRPS.sel(latitude = slice(lat_bottom,lat_top), longitude = slice(lon_min,lon_max))
    ## Evaporation dataset
    E_Ensemble = xr.open_mfdataset('/Evaporation.Ensemble_(equal_wgt).FLUXCOM(RS).BESS.PML_0.25res_daily_*.nc').Evaporation.sel(time = slice(str(year),str(year+1)))
    E_Ensemble = E_Ensemble.sel(lat = slice(lat_bottom,lat_top), lon = slice(lon_min,lon_max))
    
############################################################################################################
def Ensemble_RZSC(year,lat_top = 50, lat_bottom = -50, lon_min = 0, lon_max = 360):
    data_analysis(year,lat_top, lat_bottom, lon_min, lon_max)
    
    sum_of_days = [0,31,59,90,120,151,181,212,243,273,304,334]
    global max_deficit_annual
    print('Calculation for year: '+ str(year))
    print('Dataset: Ensemble')

    deficit_all = np.array(E_Ensemble) - np.array(P_CHIRPS)
    max_deficit_annual = np.zeros((P_CHIRPS_argmax.shape[0],P_CHIRPS_argmax.shape[1]))
    for lat in tqdm(range(P_CHIRPS_argmax.shape[0])):
        for long in range(P_CHIRPS_argmax.shape[1]):
            cummul_deficit = np.zeros((730))
            for i in range(sum_of_days[P_CHIRPS_argmax[lat,long]],sum_of_days[P_CHIRPS_argmax[lat,long]]+365):
                cummul_deficit[i] = deficit_all[i,lat,long] + cummul_deficit[i-1]
                if cummul_deficit[i] < 0:
                    cummul_deficit[i] = 0
                else:
                    continue
            max_deficit_annual[lat,long] = np.nanmax(cummul_deficit)
            
    import datetime
    from netCDF4 import Dataset,num2date,date2num
    # -----------------------
    nyears = 1;
    unout = 'days since '+str(year)+'-01-01 00:00:00'
    # -----------------------
    ny, nx = (400, 1440)
    lon = np.linspace(0.125,359.875,nx);
    lat = np.linspace(-49.875,49.875,ny);

    dataout = max_deficit_annual; # create some random data
    datesout = [datetime.datetime(int(year)+iyear,1,1) for iyear in range(nyears)]; # create datevalues
    # =========================
    ncout = Dataset('/home/chandra/data/Max_RZSC_annual_Chirps_Ensemble(BESS+PML+FLUXCOM)/Simulation7 (Sensitivity dataset)/Ensemble/Max_Rootzone_Ensemble_(BESS+PML+FLUXCOM)_Chirps_0.25res'+str(year)+'.nc', 'w','NETCDF4'); 
    ncout.createDimension('lon',nx);
    ncout.createDimension('lat',ny);
    ncout.createDimension('time',nyears);
    lonvar = ncout.createVariable('lon','float32',('lon'));lonvar[:] = lon;
    latvar = ncout.createVariable('lat','float32',('lat'));latvar[:] = lat;
    timevar = ncout.createVariable('time','float64',('time'));timevar.setncattr('units',unout);timevar[:]=date2num(datesout,unout);
    # Do not add space between names
    myvar = ncout.createVariable('RootZone_SC','float32',('time','lat','lon'));myvar.setncattr('units','mm');myvar[:] = dataout;
    ncout.close();
    
for year in range(2001,2012):
    Ensemble_RZSC(year)