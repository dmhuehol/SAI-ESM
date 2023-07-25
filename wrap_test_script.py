# cdo yearmonmean adaptor.mars.internal-1688061955.4278178-2414-15-51321947-5c9a-4619-ba29-38546b8cc74f.nc era5_reanalysis_tas_199601-201512_annual.nc
# ncrename -d longitude,lon -d latitude,lat era5OrigDim_reanalysis_tas_199601-201512_annual.nc era5_reanalysis_tas_199601-201512_annual_LON.nc
from icecream import ic
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

dp = '/Users/dhueholt/Documents/ecology_data/monthly_t2m/Land_and_Ocean_LatLong1.nc'
d = xr.open_dataset(dp)
ic(d)
sys.exit('STOP')
time = d.time
da = d.temperature
annmn = da.rolling(time=12).mean()
newDset = xr.Dataset(
    {outKey: (("lat","lon"), landmask)},
    coords={
        "lat": d['latitude'],
        "lon": d['longitude']
    }
)
newDset[outKey].attrs['long_name'] = 'BEST temperature'
newDset[outKey].attrs['units'] = 'dimensionless'
# ic(time)
# ic(np.diff(time))
# ic(d.rolling(time=12).mean().time)

# cdo mergetime cru_ts4*.nc cruts4_station_tas_199101-202012_mergetime.nc
# cdo yearmonmean cruts4_station_tas_199101-202012_mergetime.nc cruts4_station_tas_199101-202012_annual.nc
# cdo yearmonmean Land_and_Ocean_LatLong1.nc best_augment_tas_1850-2023_annual.nc