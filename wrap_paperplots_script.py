''' wrap_paperplots_script
Runs ensemble plotting functions in fun_ens_plot. Runs all figures used in ongoing
GLENS/ARISE summary paper in development.

dataDict is for data inputs
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
import numpy as np

import fun_ens_plot as fep
import fun_convert_unit as fcu
import fun_process_data as fpd
import region_library as rlib

# Call regions
ipccWg1Ar5 = rlib.atlas_ipcc_wg1ar5() #ipccWg1Ar5["allRegions"]
seaIcyRegions = rlib.atlas_seaicy_regions()
planetary = ('global', rlib.Tropics(), rlib.NorthernHemisphere(), rlib.SouthernHemisphere())
insets = rlib.atlas_insets()

# Dictionaries
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/annual_TREFHT/',
    "idGlensCntrl": 'control_*', #'control_*' or None
    "idGlensFdbck": 'feedback_*', #'feedback_*' or None
    "idArise": '*SSP245-TSMLT-GAUSS*', #'*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": '*BWSSP245*', #'*BWSSP245*' or None
    "idS245Hist": '*BWHIST*', #'*BWHIST*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc' #cesm_component_mask.nc
}
setDict = {
    "landmaskFlag": None, #None or 'land'
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "cntrlPoi": [2011,2041], #pdf
    "fdbckPoi": [2041], #pdf
    "arisePoi": [2041], #pdf
    "s245CntrlPoi": [2041], #pdf
    "timePeriod": 20, #pdf
    "plotStyle": 'step', #pdf
    "areaAvgBool": True,
    "dimOfVrblty": {'rlzBool':True,'timeBool':True,'spcBool':False}, #pdf
    "convert": (fcu.kel_to_cel,), #TUPLE of converter(s), or None if using default units
    "realization": 'ensplot',
    "insetFlag": 2, #ts 0 default, 1 lines only, 2 AMS style
    "mute": False, #ts True/False to use image muting on parts of time period
    "ylim": [14, 18.25], #ts None for automatic, [min,max] for manual
    "ylabel": '', #ts
    "yticks": np.arange(10,25,1), #ts
    "xticks": True
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20220902_hi/testpaper/',
    "dpiVal": 400
}
loopDict = {
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": (rlib.EastAfrica(),)
}

ic(dataDict, setDict, outDict) #Lowers chances of making the wrong plots by mistake

# # Make images
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# # Annual mean temperature
# setDict["regOfInt"] = 'global'
# setDict["xticks"] = False
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# setDict["regOfInt"] = rlib.Amazon()
# setDict["ylim"] = [25.5,30.5]
# setDict["yticks"] = np.arange(20,40,1)
# setDict["xticks"] = False
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# setDict["regOfInt"] = rlib.NorthEurope()
# setDict["ylim"] = [2,10.5]
# setDict["yticks"] = np.arange(0,20,2)
# setDict["xticks"] = True
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # Annual mean precipitation
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_PRECT/'
# setDict["convert"] = (fcu.m_to_cm, fcu.persec_peryr)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["regOfInt"] = 'global'
# setDict["ylim"] = [105,120]
# setDict["yticks"] = np.arange(105,130,5)
# setDict["xticks"] = False
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# setDict["regOfInt"] = rlib.AlaskaNorthwestCanada()
# setDict["ylim"] = [43,75]
# setDict["yticks"] = np.arange(43,80,10)
# setDict["xticks"] = False
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # NH monsoon seasonal mean (JJAS) precipitation
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/eimnsn_PRECT/anmn/'
# setDict["convert"] = (fcu.m_to_cm, fcu.persec_peryr)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["regOfInt"] = rlib.GeenEtAl20AsianMonsoonRegion()
# setDict["ylim"] = [180,345]
# setDict["yticks"] = np.arange(180,500,50)
# setDict["xticks"] = True
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # Annual mean sea surface temperature
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_OCNTEMP500/'
# setDict["convert"] = None
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["regOfInt"] = 'global'
# setDict["ylim"] = [18,21.25]
# setDict["yticks"] = np.arange(10,30,1)
# setDict["xticks"] = False
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # September Arctic sea ice extent
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/sept_ICEEXTENT/'
# setDict["areaAvgBool"] = 'sum'
# setDict["convert"] = (fcu.km2_to_milkm2,)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["regOfInt"] = rlib.Arctic()
# setDict["ylim"] = [0,9.5]
# setDict["yticks"] = np.arange(0,10,2)
# setDict["xticks"] = False
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # February Antarctic sea ice extent
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/feb_ICEEXTENT/'
# setDict["areaAvgBool"] = 'sum'
# setDict["convert"] = (fcu.km2_to_milkm2,)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["regOfInt"] = rlib.Antarctica()
# setDict["ylim"] = [2.5,10]
# setDict["yticks"] = np.arange(2.5,40,2)
# setDict["xticks"] = True
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # Mid-latitude annual tropical nights
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/extreme_clxTR/'
# setDict["areaAvgBool"] = True
# setDict["convert"] = None
# setDict["landmaskFlag"] = 'land'
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["regOfInt"] = (rlib.NorthernHemisphereMidLat(), rlib.SouthernHemisphereMidLat())
# setDict["ylim"] = [28,76]
# setDict["yticks"] = np.arange(28,100,15)
# setDict["xticks"] = False
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # East African SDII
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/extreme_sdii/'
# setDict["areaAvgBool"] = True
# setDict["convert"] = None
# setDict["landmaskFlag"] = None
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["regOfInt"] = (rlib.EastAfricaAyugiEtAl(),)
# setDict["ylim"] = [5.5,9.5]
# setDict["yticks"] = np.arange(6,30,1)
# setDict["xticks"] = False
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # W Australian MHWs
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/extreme_MHW/roll5_centerFalseShift4/'
# setDict["areaAvgBool"] = True
# setDict["convert"] = (fcu.roll5_to_percentdays,)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["landmaskFlag"] = None
# setDict["regOfInt"] = (rlib.WesternAustraliaMHW_point(),)
# setDict["ylim"] = [0,100]
# setDict["yticks"] = np.arange(0,101,30)
# setDict["xticks"] = True
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# # September sea ice thickness
# dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/sept_hi/'
# setDict["areaAvgBool"] = True
# setDict["convert"] = (fcu.m_to_cm_ice,)
# scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
# dataDict = {**dataDict, **cmnDict}
# setDict["levOfInt"] = None
#
# setDict["landmaskFlag"] = None
# setDict["regOfInt"] = rlib.Arctic()
# setDict["ylim"] = [0,110]
# setDict["yticks"] = np.arange(0,110,20)
# setDict["xticks"] = True
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)
#
# setDict["regOfInt"] = rlib.SouthernOcean()
# setDict["ylim"] = [30,75]
# setDict["yticks"] = np.arange(0,80,10)
# setDict["xticks"] = True
# fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

# February sea ice thickness
dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/feb_hi/'
setDict["areaAvgBool"] = True
setDict["convert"] = (fcu.m_to_cm_ice,)
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None

setDict["landmaskFlag"] = None
setDict["regOfInt"] = rlib.Arctic()
setDict["ylim"] = [0,155]
setDict["yticks"] = np.arange(0,200,25)
setDict["xticks"] = True
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

setDict["regOfInt"] = rlib.SouthernOcean()
setDict["ylim"] = [6,36]
setDict["yticks"] = np.arange(0,50,6)
setDict["xticks"] = True
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

endMsg = ic('Completed! :D')
