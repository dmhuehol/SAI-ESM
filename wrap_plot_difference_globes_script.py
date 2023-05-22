''' wrap_basicplots_script
Runs map plotting functions in fun_basic_plot.

dataDict: defines input data
setDict: settings for analysis/visualization
    plotPanel: determines which panel to plot, valid entries given below
        with relevant figure numbers from Hueholt et al. 2023
        'snapR85': snapshot for RCP8.5, Fig. 1a
        'snapS245': snapshot for SSP2-4.5, Fig. 1b
        'snapGLENS': snapshot around deployment for GLENS, Fig. 3a
        'snapARISE15': snapshot around deployment for ARISE-SAI-1.5, Fig. 3b
        'snapARISEDS': snapshot around deployment for ARISE-SAI-1.5-DelayedStart
        'snapARISE10': snapshot around deployment for ARISE-SAI-1.0
        'intiGLENS': intervention impact for GLENS, Fig. 6a
        'intiARISE15': intervention impact for ARISE-SAI-1.5, Fig. 6b
        'intiARISEDS': intervention impact for ARISE-SAI-1.5-DelayedStart
        'intiARISE10': intervention impact for ARISE-SAI-1.0
        'snapUKS245': snapshot for SSP2-4.5 in UKESM-ARISE
        'snapUKARISE15': snapshot around deployment for UKESM-ARISE-SAI-1.5
        'intiUKARISE15': intervention impact for UKESM-ARISE-SAI-1.5
        any other value: make a blank map
outDict: output image settings
loopDict: determines which images are made

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

from matplotlib import cm
import cmasher, cmocean, seaborn  # Colormap packages

import fun_plot_difference_globes as fpdg
import fun_convert_unit as fcu
import fun_process_data as fpd
import region_library as rlib

# Special color palettes
tropicalPal = seaborn.diverging_palette(133, 324, as_cmap=True)
precipPal = seaborn.diverging_palette(58, 162, s=100, l=45, as_cmap=True)

# Dictionaries to define inputs
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/ecology_data/annual_2mTemp/',
    "idGlensCntrl": None,  # 'control_*' or None
    "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": '*DEFAULT*',  # '*DEFAULT*' or None
    "idS245Cntrl": '*BWSSP245*',  # '*BWSSP245*' or None
    "idS245Hist": None,  # '*BWHIST*' or None
    "idUkesmNoSai": None, #'*ssp245*' or None
    "idUkesmArise": None, #'*arise-sai-1p5*' or None
    "idDelayedStart": None, # '*DELAYED*' or None
    "idArise1p0": None, # '*ARISE1P0*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary_landmask.nc' #Landmask file location (UKESM)
}
setDict = {
    "landmaskFlag": None,  # None no mask, 'land' to mask ocean, 'ocean' to mask land
    "strtIntvl": { # Window years for starting interval [start,end)
        "GLENS": [2015,2020],
        "CESM2-ARISE": [2030,2035],
        "UKESM-ARISE": [2030,2035]
        },
    "endIntvl": { # Window years for ending interval
        "GLENS": [2025,2030],
        "CESM2-ARISE": [2065,2070],
        "UKESM-ARISE": [2065,2070]
        },
    "convert": (fcu.kel_to_cel,),  # TUPLE of converter(s), None for default units
    "cmap": None,  # None for default (cmocean balance) or choose colormap
    "cbVals": [-2,2],  # None for automatic or [min,max] to override,
    "addCyclicPoint": False,  # True for ocean data/False for others
    "areaAvgBool": False,  # ALWAYS FALSE: no area averaging for a map!
    "robustnessBool": False,  # True/False to run robustness
    "plotPanel": 'snapARISE15' # See docstring for valid inputs
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/ecology_fig/20230420_ukesm/',
    "dpiVal": 400
}
loopDict = {
    "rlzs": ('ensplot',),  # number(s) for member(s), 'mean' ens mean all members, 'ensplot' for both member information and mean (i.e. for robustness)
    "levels": (None,),  # None for single-level variable (as for all used in Hueholt et al. 2023)
    "regions": ('global',),  # 'global' only for maps
}
ic(dataDict, setDict, outDict, loopDict)  # Show input settings at command line

# Make images
for rlz in loopDict["rlzs"]:
    setDict["realization"] = rlz
    scnList, cmnDict = fpd.call_to_open(dataDict, setDict)
    # ic(scnList, cmnDict)
    dataDict = {**dataDict, **cmnDict}
    # ic(dataDict)
    # sys.exit('STOP')

    for lev in loopDict["levels"]:
        setDict["levOfInt"] = lev
        fpdg.plot_single_basic_difference_globe(
            scnList, dataDict, setDict, outDict)
        # fpdg.plot_single_robust_globe(
        #    scnList, dataDict, setDict, outDict)