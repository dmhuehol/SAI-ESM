''' wrap_paperplots_script
Runs ensemble plotting functions in fun_ens_plot. Runs all figures used in ongoing
GLENS/ARISE summary paper in development.

dataDict is for data inputs
setDict sets settings related to plotting
outDict is for outputs
loopDict determines which images are made

Written by Daniel Hueholt | February 2022
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
    "dimOfVrblty": {'rlzBool':True,'timeBool':True,'spcBool':False}, #pdf
    "convert": (fcu.kel_to_cel,), #TUPLE of converter(s), or None if using default units
    "realization": 'ensplot',
    "insetFlag": 2, #ts 0 default, 1 lines only, 2 AMS style
    "mute": False, #ts True/False to use image muting on parts of time period
    "ylim": [14, 19], #ts None for automatic, [min,max] for manual
    "yticks": np.arange(0,50,1),
    "xticks": True
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/GLENS_fig/20220210_designingFinalFigs/4_aspect/',
    "dpiVal": 400
}
loopDict = {
    "levels": (None,), #'stratosphere', 'troposphere', 'total', numeric level(s), or None for surface variable
    "regions": (rlib.EastAfrica(),)
}

ic(dataDict, setDict, outDict) #Lowers chances of making the wrong plots by mistake

# Make images
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None

# Annual mean temperature
setDict["regOfInt"] = 'global'
setDict["xticks"] = False
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

setDict["regOfInt"] = rlib.Amazon()
setDict["ylim"] = [25.5,30.5]
setDict["xticks"] = False
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

setDict["regOfInt"] = rlib.NorthEurope()
setDict["ylim"] = [3,10]
setDict["yticks"] = np.arange(0,20,2)
setDict["xticks"] = True
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

# Annual mean precipitation
dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_PRECT/'
setDict["convert"] = (fcu.m_to_cm, fcu.persec_peryr)
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None

setDict["regOfInt"] = 'global'
setDict["ylim"] = [105,120]
setDict["yticks"] = np.arange(105,130,5)
setDict["xticks"] = False
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

setDict["regOfInt"] = rlib.EastAfrica()
setDict["ylim"] = [70,190]
setDict["yticks"] = np.arange(70,500,35)
setDict["xticks"] = False
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

# NH monsoon seasonal mean (JJAS) precipitation
dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/eimnsn_PRECT/anmn/'
setDict["convert"] = (fcu.m_to_cm, fcu.persec_peryr)
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None

setDict["regOfInt"] = rlib.GeenEtAl20AsianMonsoonRegion()
setDict["ylim"] = [180,345]
setDict["yticks"] = np.arange(180,500,50)
setDict["xticks"] = True
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

# Annual mean sea surface temperature
dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/annual_OCNTEMP500/'
setDict["convert"] = None
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None

setDict["regOfInt"] = 'global'
setDict["ylim"] = [18,22]
setDict["yticks"] = np.arange(10,30,1)
setDict["xticks"] = False
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

# September Arctic sea ice thickness
dataDict["dataPath"] = '/Users/dhueholt/Documents/GLENS_data/sept_hi/'
setDict["convert"] = (fcu.m_to_cm,)
scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
dataDict = {**dataDict, **cmnDict}
setDict["levOfInt"] = None

setDict["regOfInt"] = rlib.Arctic()
setDict["ylim"] = [0,105]
setDict["yticks"] = np.arange(0,150,25)
setDict["xticks"] = False
fep.plot_ens_spread_timeseries(scnList, dataDict, setDict, outDict)

endMsg = ic('Completed! :D')
