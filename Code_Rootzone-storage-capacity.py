import xarray as xr
import numpy as np
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")

def sk_extremes(method = 'mle'):
    RZSC_annual = xr.open_mfdataset('/Max_Rootzone_Ensemble_(BESS+PML+FLUXCOM)_Chirps_0.25res*.nc').RootZone_SC
    print('Calculating for Ensemble')
    import skextremes as ske
    global mle_gumbel, mle_CI_lower, mle_CI_upper
    
    if str(method) == 'mle':
        ci_method = 'delta'
    else:
        ci_method = 'bootstrap'
    print('Extreme method: '+str(method))
    mle_gumbel = np.zeros((RZSC_annual.values.shape[1],RZSC_annual.values.shape[2]))
    mle_gumbel[:] = np.nan
    mle_CI_lower = np.zeros((RZSC_annual.values.shape[1],RZSC_annual.values.shape[2]))
    mle_CI_lower[:] = np.NaN
    mle_CI_upper = np.zeros((RZSC_annual.values.shape[1],RZSC_annual.values.shape[2]))
    mle_CI_upper[:] = np.NaN
    WD_nanmean = np.nanmean(RZSC_annual.values, axis = 0)
    
    for lat in tqdm(range(RZSC_annual.values.shape[1])):
        for lon in range(RZSC_annual.values.shape[2]):
            WD_cal = (RZSC_annual[:,lat,lon].values)
            if np.isnan(WD_nanmean[lat,lon]) == True:
                continue
            if np.nanmean(WD_nanmean[lat,lon]) == 0:
                continue
            else:
                model1 = ske.models.classic.Gumbel(WD_cal, fit_method = str(method), return_periods=20, ci = 0.05, ci_method = str(ci_method));
                mle_gumbel[lat,lon] = (model1.return_values)
                mle_CI_upper[lat,lon] = (model1._ci_Td[200])
                mle_CI_lower[lat,lon] = (model1._ci_Tu[200])
    
    year = 2001
    import datetime
    from netCDF4 import Dataset,num2date,date2num
    # -----------------------
    nyears = 1;
    unout = 'days since '+str(year)+'-01-01 00:00:00'
        # -----------------------
    ny, nx = (400, 1440)
    lon = np.linspace(0.125,359.875,nx);
    lat = np.linspace(-49.875,49.875,ny);

    dataout = mle_gumbel; 
    datesout = [datetime.datetime(int(year)+iyear,1,1) for iyear in range(nyears)]; # create datevalues
        # =========================
    ncout = Dataset('Max_Rootzone_20yearsreturn_Ensemble_(BESS+PML+FLUXCOM)_Chirps_0.25res_2001-2012_'+str(method)+'(skextremes).nc', 'w','NETCDF4'); # using netCDF3 for output format 
    ncout.createDimension('lon',nx);
    ncout.createDimension('lat',ny);
    ncout.createDimension('time',nyears);
    lonvar = ncout.createVariable('lon','float32',('lon'));lonvar[:] = lon;
    latvar = ncout.createVariable('lat','float32',('lat'));latvar[:] = lat;
    timevar = ncout.createVariable('time','float64',('time'));timevar.setncattr('units',unout);timevar[:]=date2num(datesout,unout);
    # Do not add space between names
    myvar = ncout.createVariable('mle_gumbel','float32',('time','lat','lon'));myvar.setncattr('units','mm');myvar[:] = dataout;
    ncout.close();
    
    # -----------------------
    nyears = 1;
    unout = 'days since '+str(year)+'-01-01 00:00:00'
        # -----------------------
    ny, nx = (400, 1440)
    lon = np.linspace(0.125,359.875,nx);
    lat = np.linspace(-49.875,49.875,ny);

    dataout = mle_CI_lower; 
    datesout = [datetime.datetime(int(year)+iyear,1,1) for iyear in range(nyears)]; # create datevalues
        # =========================
    ncout = Dataset('Max_Rootzone_20yearsreturn_Ensemble_(BESS+PML+FLUXCOM)_Chirps_0.25res_2001-2012_'+str(method)+'(skextremes)_CI_0.05_lower.nc', 'w','NETCDF4'); # using netCDF3 for output format 
    ncout.createDimension('lon',nx);
    ncout.createDimension('lat',ny);
    ncout.createDimension('time',nyears);
    lonvar = ncout.createVariable('lon','float32',('lon'));lonvar[:] = lon;
    latvar = ncout.createVariable('lat','float32',('lat'));latvar[:] = lat;
    timevar = ncout.createVariable('time','float64',('time'));timevar.setncattr('units',unout);timevar[:]=date2num(datesout,unout);
    # Do not add space between names
    myvar = ncout.createVariable('mle_gumbel_CI_lower','float32',('time','lat','lon'));myvar.setncattr('units','mm');myvar[:] = dataout;
    ncout.close();

    # -----------------------
    nyears = 1;
    unout = 'days since '+str(year)+'-01-01 00:00:00'
        # -----------------------
    ny, nx = (400, 1440)
    lon = np.linspace(0.125,359.875,nx);
    lat = np.linspace(-49.875,49.875,ny);

    dataout = mle_CI_upper; 
    datesout = [datetime.datetime(int(year)+iyear,1,1) for iyear in range(nyears)]; # create datevalues
        # =========================
    ncout = Dataset('Max_Rootzone_20yearsreturn_Ensemble_(BESS+PML+FLUXCOM)_Chirps_0.25res_2001-2012_'+str(method)+'(skextremes)_CI_0.05_upper.nc', 'w','NETCDF4'); # using netCDF3 for output format 
    ncout.createDimension('lon',nx);
    ncout.createDimension('lat',ny);
    ncout.createDimension('time',nyears);
    lonvar = ncout.createVariable('lon','float32',('lon'));lonvar[:] = lon;
    latvar = ncout.createVariable('lat','float32',('lat'));latvar[:] = lat;
    timevar = ncout.createVariable('time','float64',('time'));timevar.setncattr('units',unout);timevar[:]=date2num(datesout,unout);
    # Do not add space between names
    myvar = ncout.createVariable('mle_gumbel_CI_upper','float32',('time','lat','lon'));myvar.setncattr('units','mm');myvar[:] = dataout;

sk_extremes()
