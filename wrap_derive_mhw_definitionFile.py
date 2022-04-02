'''wrap_derive_mhw_definitionFile
Defines MHW baseline for reference period at a given POINT location (typical) or
REGION (of dubious use).

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
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/daily_SST/ARISE_defsPeriod/',
    "idGlensCntrl": None, #'*control*' or None
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
    "realization": 'ensplot',
    "regOfInt": (rlib.WesternAustraliaMHW_point(),),
}
outDict = {
    "outKeyMn": "mn_SST",
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20220401_mhwLastSteps/2_defHistAr/',
}

for reg in setDict["regOfInt"]:
    darrList, cmnDict = fpd.call_to_open(dataDict, setDict) #Long because of I/O overhead on daily data
    for darr in darrList: #for each scenario
        darrPoint, locStr, locTitleStr = fpd.manage_area(darr, reg, areaAvgBool=setDict["areaAvgBool"])
        darrPointFullTimes = darrPoint.resample(time='1D').asfreq() #Add missing timesteps with NaN value
        rlzMnInd = len(darrPointFullTimes.realization)-1 #last realization index is mean
        rmn = darrPointFullTimes.isel(realization=rlzMnInd).data.squeeze()

        times = darrPointFullTimes.time.data
        ordArr = fpd.make_ord_array(times)

        rlzMnClimDset = xr.Dataset(
            {outDict["outKeyMn"]: (("time"), rmn)},
            coords={
                "time": (('time'), ordArr),
            }
        )
        rlzMnClimDset[outDict["outKeyMn"]].attrs = darrPointFullTimes.attrs
        rlzMnClimDset[outDict["outKeyMn"]].attrs['long_name'] = 'Ensemble mean SST'
        rlzMnClimDset[outDict["outKeyMn"]].attrs['lat'] = reg["regLats"]
        rlzMnClimDset[outDict["outKeyMn"]].attrs['lon'] = reg["regLons"]

        if 'GLENS:Control' in darrPointFullTimes.scenario:
            scnStr = 'GLENS'
        elif 'ARISE:Control' in darrPointFullTimes.scenario:
            scnStr = 'ARISE'
        else:
            ic('Unknown scenario!')
            scnStr = 'unknown'

        strOutRlzMn = 'mhwDefsFile_' + scnStr + '_' + reg["regSaveStr"]
        outFileRlzMn = outDict["savePath"] + strOutRlzMn + '.nc'
        rlzMnClimDset.to_netcdf(outFileRlzMn) #Save data
        ic(strOutRlzMn, outFileRlzMn, rlzMnClimDset) #icecream all useful parts of the output
