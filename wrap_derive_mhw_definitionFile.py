'''wrap_derive_mhw_definitionFile
Saves an MHW definitions file for a given POINT location (typical), REGION (of
dubious use), or all points globally. An MHW definitions file is just the
timeseries of the ensemble mean SST at a location for the reference period
2010-2019.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import xarray as xr
import numpy as np
from datetime import date
import marineHeatWaves as mhws

import fun_process_data as fpd
import region_library as rlib

# Note that data ID keys are slightly different than most wrap_ scripts because
# these have raw-type filenames!
dataDict = {
    "dataPath": '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/GLENS/', #'/Users/dhueholt/Documents/GLENS_data/daily_SST/GLENS_defsPeriod/',
    "idGlensCntrl": '*control*', #'*control*' or None
    "idGlensFdbck": None, #'*feedback*' or None
    "idArise": None, #'*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": '*BWSSP245*', # d'*BWSSP245*' or None
    "idS245Hist": '*BWHIST*', #'*BWHIST*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc' #cesm_component_mask.nc
}
setDict = {
    "landmaskFlag": None, #None or 'land'
    "areaAvgBool": True, #True=mean, 'sum'=sum (e.g. ice extent), False=functionality that hasn't been tested in two years, and yes I know that's not a boolean
    "convert": None, #TUPLE of converter(s), or None if using default units
    "realization": 'mean',
    "regOfInt": (rlib.Globe(),),
}
outDict = {
    "outKeyMn": "mn_sst",
    "savePath": '/glade/work/dhueholt/definitionFiles/'#'/Users/dhueholt/Documents/GLENS_fig/20220421_mhwsGlobal/',
}

darrList, cmnDict = fpd.call_to_open(dataDict, setDict) #Long because of I/O overhead on daily data
for reg in setDict["regOfInt"]:
    for darr in darrList: #for each scenario
        darrFullTimes = darr.resample(time='1D').asfreq() #Add missing timesteps with NaN value
        rmn = darrFullTimes.data

        lats = darrFullTimes.lat.data
        lons = darrFullTimes.lon.data
        times = darrFullTimes.time.data
        ordArr = fpd.make_ord_array(times)

        rlzMnClimDset = xr.Dataset(
            {outDict["outKeyMn"]: (("time","lat","lon"), rmn)},
            coords={
                "time": (('time'), ordArr),
                "lat": lats,
                "lon": lons
            }
        )
        rlzMnClimDset[outDict["outKeyMn"]].attrs = darrFullTimes.attrs
        rlzMnClimDset[outDict["outKeyMn"]].attrs['long_name'] = 'Ensemble mean SST'
        rlzMnClimDset[outDict["outKeyMn"]].attrs['lat'] = reg["regLats"]
        rlzMnClimDset[outDict["outKeyMn"]].attrs['lon'] = reg["regLons"]

        if 'GLENS:Control' in darrFullTimes.scenario:
            scnStr = 'GLENS'
        elif 'ARISE:Control' in darrFullTimes.scenario:
            scnStr = 'ARISE'
        else:
            ic('Unknown scenario!')
            scnStr = 'unknown'

        strOutRlzMn = 'mhwDefsFile_' + scnStr + '_' + reg["regSaveStr"]
        outFileRlzMn = outDict["savePath"] + strOutRlzMn + '.nc'
        rlzMnClimDset.to_netcdf(outFileRlzMn) #Save data
        ic(strOutRlzMn, outFileRlzMn, rlzMnClimDset) #icecream all useful parts of the output
