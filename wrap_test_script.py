''' 
Remove duplicate timesteps from merged data and resave. Run cdo mergetime from
wrap_mproc_cdo_prep first.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import glob
import xarray as xr

dataDict = {
    "dataPath": '/Users/dhueholt/Documents/ecology_data/annual_TREFHT/DelayedStart/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": None,  # '*DEFAULT*' or None
    "idS245Cntrl": None,  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": None, #'*ssp245*' or None
    "idUkesmArise": None, #'*arise-sai-1p5*' or None
    "idPiControl": None, #'*piControl*' or None
    "idDelayedStart": '*DELAYED*', # '*DELAYED*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary0p01_landmask.nc' #Landmask file location (UKESM)
}
outDict = {
    "outPath": '/Users/dhueholt/Documents/ecology_data/annual_TREFHT/DelayedStart/nodup/',
    "saveId": 'nodup'
}

# Code adapted from call_to_open
for dky in dataDict.keys():
    if dataDict[dky] is None:
        continue
    if 'id' in dky:
        inPath = dataDict["dataPath"] + dataDict[dky]
        inGlobs = sorted(glob.glob(inPath))
        # ic(inGlobs)
        for inFil in inGlobs:
            inDs = xr.open_dataset(inFil)
            # ic(inDs)
            noDupDs = inDs.drop_duplicates(dim='time')
            # ic(inDs)
            
            inFn = inFil.split('/')[-1]
            strOut = inFn.replace('.nc', outDict["saveId"] + '.nc')
            noDupDs.to_netcdf(outDict["outPath"] + strOut)
            ic(strOut)