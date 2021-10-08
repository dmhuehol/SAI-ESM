''' wrap_ensplots_script
Runs ensemble plotting functions in ens_plot_fun.

dataDict is for inputs
setDict sets settings related to plotting
outDict is for outputs
loopDict determines which images are made

Written by Daniel Hueholt | August 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import cmocean
import numpy as np

import ens_plot_fun as epf
import fun_convert_unit as fcu
import process_glens_fun as pgf
import region_library as rlib

# Call regions
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5() #ipccWg1Ar5["allRegions"]
seaIcyRegions = rlib.atlas_seaicy_regions()
planetary = ('global', rlib.Tropics(), rlib.NorthernHemisphere(), rlib.SouthernHemisphere(),)
otherRegions = (rlib.NorthAtlanticWarmingHole(),rlib.AustralianContinent())
insets = (rlib.Sahara(),rlib.NorthEurope(),rlib.Amazon(),rlib.WestAsia(),
          rlib.SouthernAfrica(),rlib.WestNorthAmerica(),rlib.EastAsia(),
          rlib.CentralAmericaMexico(),rlib.SoutheastAsia(),rlib.NorthAtlanticWarmingHole(),
          rlib.AustralianContinent(),'global',)

# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_500TEMP/regrid/',
    "idGlensCntrl": 'control_*', #'control_*'
    "idGlensFdbck": 'feedback_*', #'feedback_*'
    "idSciris": '*SSP245-TSMLT-GAUSS*', #'*SSP245*'
    "idS245Cntrl": '*BWSSP245*', #'*ssp245*'
    "idS245Hist": '*BWHIST*', #'*historical*'
    "idMask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc' #cesm_component_mask.nc
}
setDict = {
    "landmaskFlag": None, #None or 'land'
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "scirisPoi": [2041], #pdf
    "s245CntrlPoi": [2041], #pdf
    "timePeriod": 20, #pdf
    "plotStyle": 'step', #pdf
    "dimOfVrblty": {'rlzBool':True,'timeBool':True,'spcBool':False}, #pdf
    "convert": None, #TUPLE of converter(s), or None if using default units
    "realization": 'ensplot',
    "insetFlag": 2
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20211008_refactoring/',
    "dpiVal": 800
}
loopDict = {
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": ('global',)
}

# Verify inputs (troubleshooting)
ic(dataDict)
ic(setDict)
ic(outDict)
# ic(setDict["convert"])

# Make images
rlzList, cmnDict = pgf.call_to_open(dataDict, setDict)
# rlzList = pgf.period_month_avg(rlzList)
dataDict = {**dataDict, **cmnDict}
for lev in loopDict["levels"]:
    setDict["levOfInt"] = lev
    for reg in loopDict["regions"]:
        setDict["regOfInt"] = reg
        # epf.plot_ens_spaghetti_timeseries(rlzList, dataDict, setDict, outDict)
        epf.plot_ens_spread_timeseries(rlzList, dataDict, setDict, outDict)
        # setDict["plotStyle"] = 'step'
        # epf.plot_ens_pdf(rlzList, dataDict, setDict, outDict)
        # setDict["plotStyle"] = 'kde'
        # epf.plot_ens_pdf(rlzList, dataDict, setDict, outDict)
        # setDict["plotStyle"] = 'hist'
        # epf.plot_ens_pdf(rlzList, dataDict, setDict, outDict)
