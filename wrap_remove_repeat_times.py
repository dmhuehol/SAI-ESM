''' wrap_remove_repeat_times
At times, pre-release CESM simulations can contain incorrect repeated timesteps in
the output files. These repeated timesteps will throw off code that expects time 
to proceed strictly linearly e.g., cdo yearmonmean. This script opens files,
removes the repeated timesteps from merged data, then resaves the data. 

Run cdo mergetime from wrap_mproc_cdo_prep before this script--otherwise, the code
cannot identify the repeated timesteps.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import glob
import xarray as xr

dataDict = {
    "dataPath": '/Users/dhueholt/Documents/ecology_data/annual_TREFHT/arise1p0/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": None,  # '*DEFAULT*' or None
    "idS245Cntrl": None,  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": None, #'*ssp245*' or None
    "idUkesmArise": None, #'*arise-sai-1p5*' or None
    "idPiControl": None, #'*piControl*' or None
    "idDelayedStart": None, # '*DELAYED*' or None
    "idArise1p0": '*ARISE1P0*' # '*ARISE1P0*' or None
}
outDict = {
    "outPath": '/Users/dhueholt/Documents/ecology_data/annual_TREFHT/arise1p0/nodup/',
    "saveId": 'nodup' # Bonus string to add to distinguish output files
}

for dky in dataDict.keys():
    if dataDict[dky] is None:
        continue
    if 'id' in dky:
        inPath = dataDict["dataPath"] + dataDict[dky]
        inGlobs = sorted(glob.glob(inPath))
        for inFil in inGlobs:
            inDs = xr.open_dataset(inFil)
            noDupDs = inDs.drop_duplicates(dim='time')
            
            inFn = inFil.split('/')[-1]
            strOut = inFn.replace('.nc', outDict["saveId"] + '.nc')
            noDupDs.to_netcdf(outDict["outPath"] + strOut)
            ic(strOut)