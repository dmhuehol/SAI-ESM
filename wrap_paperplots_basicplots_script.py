''' wrap_paperplots_basicplots_script
Make difference globe figures for GLENS/ARISE summary paper in development.
Figures are produced in 1-panel format with no annotations. Keynote is used to
stitch panels together, add title, add colorbar, etc.

dataDict is for inputs
setDict sets settings related to plotting
outDict is for outputs
loopDict determines which images are made

Written by Daniel Hueholt | April 2022
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import cmocean
import cmasher
import seaborn
import numpy as np

import fun_basic_plot as fbp
import fun_convert_unit as fcu
import fun_process_data as fpd
import region_library as rlib

# Call regions
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5() #ipccWg1Ar5["allRegions"]
testAllTypes = rlib.atlas_all_types() #testAllTypes["allRegions"]
gnsht = ('global', rlib.Arctic(), rlib.HudsonBay(), rlib.NorthernHemisphere(), rlib.SouthernHemisphere(),)

# Specials
tropicalPal = seaborn.diverging_palette(133, 324, as_cmap=True)
indRedPal = seaborn.diverging_palette(16.8, 270.2, s=100, l=40, as_cmap=True)
precipPal = seaborn.diverging_palette(58, 162, s=100, l=45, as_cmap=True)
xtPrecipPal = seaborn.diverging_palette(58, 162, s=100, l=30, as_cmap=True)

# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_TREFHT/',
    "idGlensCntrl": 'control_*', #'control_*' or None
    "idGlensFdbck": 'feedback_*', #'feedback_*' or None
    "idArise": '*SSP245-TSMLT-GAUSS*', #'*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": '*BWSSP245*', #'*BWSSP245*' or None
    "idS245Hist": '*BWHIST*', #'*BWHIST*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc'
}
setDict = {
    "landmaskFlag": None,
    "startIntvl": [2015,2020,2030,2035], #dg [2015,2020,2030,2035]
    "endIntvl": [2025,2030,2040,2045], #dg [2025,2030,2040,2045]
    "convert": (fcu.kel_to_cel,), #TUPLE of converter(s), or None if using default units
    "cmap": cmocean.cm.balance, #None for default cmocean "balance" or choose colormap here
    "cbVals": [-2,2], #None for automatic or [min,max] to override #dg
    "addCyclicPoint": False,
    "areaAvgBool": False,
    "robustnessBool": True
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20221201_sig/',
    "dpiVal": 800 #High-res for paper
}
loopDict = {
    "realizations": ('mean',), #number for individual member, 'mean' for ens mean of all available members
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": ('global',),#('global',rlib.Arctic(),rlib.EastNorthAmerica()),
    "aaBools": (True,)
}

ic(dataDict, setDict, outDict) #Lowers chances of making the wrong plots by mistake

# # Make images
setDict["realization"] = 'mean'
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None

# Annual mean temperature
fbp.plot_paper_panels_globe(scnList, dataDict, setDict, outDict)

# Annual mean precipitation
dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_PRECT/'
setDict["convert"] = (fcu.m_to_cm, fcu.persec_peryr)
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None
setDict["cmap"] = precipPal
setDict["cbVals"] = [-25,25]
fbp.plot_paper_panels_globe(scnList, dataDict, setDict, outDict)

# Annual mean SDII
dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/extreme_sdii/'
setDict["areaAvgBool"] = True
setDict["convert"] = None
setDict["landmaskFlag"] = 'land'
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None

setDict["cmap"] = precipPal
setDict["cbVals"] = [-1,1]
fbp.plot_paper_panels_globe(scnList, dataDict, setDict, outDict)

# SUPPLEMENTAL: Robustness figures
# TODO: MUST have robustness dictionary as an input somewhere--otherwise cannot
# automatically generate figures for both GLENS and ARISE!
# dataDict["idArise"] = None
# dataDict["idS245Cntrl"] = None
# dataDict["idS245Hist"] = None
# setDict["realization"] = 'ensplot'
# setDict["robustnessBool"] = True
# setDict["cmap"] = None #cmasher ghostlight is used automatically
# setDict["cbVals"] = None #proper range set automatically
# # Robustness: Annual mean temperature (GLENS)
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_TREFHT/'
# setDict["convert"] = (fcu.kel_to_cel,)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# fbp.plot_paper_robust_globe(scnList, dataDict, setDict, outDict)
# # Robustness: Annual mean precipitation (GLENS)
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_PRECT/'
# setDict["convert"] = (fcu.m_to_cm, fcu.persec_peryr)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# fbp.plot_paper_robust_globe(scnList, dataDict, setDict, outDict)
# # Robustness: Annual simple intensity index (GLENS)
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/extreme_sdii/'
# setDict["convert"] = None
# setDict["landmaskFlag"] = 'land'
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# fbp.plot_paper_robust_globe(scnList, dataDict, setDict, outDict)
# # Robustness: Annual mean temperature (ARISE)
# dataDict["idArise"] = '*SSP245-TSMLT-GAUSS*'
# dataDict["idS245Cntrl"] = '*BWSSP245*'
# dataDict["idS245Hist"] = '*BWHIST*'
# dataDict["idGlensCntrl"] = None
# dataDict["idGlensFdbck"] = None
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_TREFHT/'
# setDict["convert"] = (fcu.kel_to_cel,)
# setDict["landmaskFlag"] = None
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# fbp.plot_paper_robust_globe(scnList, dataDict, setDict, outDict)
# # Robustness: Annual mean precipitation (ARISE)
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_PRECT/'
# setDict["convert"] = (fcu.m_to_cm, fcu.persec_peryr)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# fbp.plot_paper_robust_globe(scnList, dataDict, setDict, outDict)
# # Robustness: Annual simple intensity index (ARISE)
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/extreme_sdii/'
# setDict["convert"] = None
# setDict["landmaskFlag"] = 'land'
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# fbp.plot_paper_robust_globe(scnList, dataDict, setDict, outDict)
