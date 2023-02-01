# Demonstrate robustness in an intuitive figure-driven way
# Written by Daniel Hueholt
from icecream import ic

import numpy as np
import xarray as xr

import fun_convert_unit as fcu
import fun_process_data as fpd
import fun_robustness as fr
import fun_special_plot as fsp
import region_library as rlib


dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_TREFHT/',
    "idGlensCntrl": None, #'control_*' or None
    "idGlensFdbck": None, #'feedback_*' or None
    "idArise": '*SSP245-TSMLT-GAUSS*', #'*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": '*BWSSP245*', #'*BWSSP245*' or None
    "idS245Hist": '*BWHIST*', #'*BWHIST*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc' #cesm_component_mask.nc
}
setDict = {
    "landmaskFlag": 'land', #None or 'land'
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "arisePoi": [2041], #pdf
    "s245CntrlPoi": [2041], #pdf
    "timePeriod": 20, #pdf
    "plotStyle": 'step', #pdf
    "areaAvgBool": True, #True=mean, 'sum'=sum (e.g. ice extent), False=functionality that hasn't been tested in two years, and yes I know that's not a boolean
    "dimOfVrblty": {'rlzBool':True,'timeBool':True,'spcBool':False}, #pdf
    "convert": (fcu.kel_to_cel,), #TUPLE of converter(s), or None if using default units
    "realization": 'ensplot',
    "styleFlag": 2, #ts 0 automatic, 1 lines only, 2 for paper, 3 for IPCC
    "mute": False, #ts True/False to use image muting on parts of time period
    "ylim": None, #ts None for automatic, [min,max] for manual
    "ylabel": '', #ts None for automatic
    "yticks": None, #ts None for automatic, np.arange(mn,mx,step) for manual
    "xticks": True #ts
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20221214_newSuppExp/',
    "dpiVal": 400
}
loopDict = {
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": ('global', ) #rlib.region() for single (rlib.region(), rlib.region()) to concatenate)
}

ic(dataDict, setDict, outDict) #Lowers chances of making the wrong plots by mistake

scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
for lev in loopDict["levels"]:
    setDict["levOfInt"] = lev
    for reg in loopDict["regions"]:
        setDict["regOfInt"] = reg
        fsp.plot_rob_spaghetti_demo(scnList, dataDict, setDict, outDict)
