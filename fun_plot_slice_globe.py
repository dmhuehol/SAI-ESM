''' fun_plot_slice_globe
Contains the globe plotting functions for single time slices. The same three
dictionary inputs (defining the input files, plot settings, and output,
respectively) are used by each function.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import numpy as np
import xarray as xr
xr.set_options(keep_attrs=True)
import cftime
import scipy.stats as stats
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn
import cmocean
import cmasher

import matplotlib.font_manager as fm
fontPath = '/Users/dhueholt/Library/Fonts/'  #Location of font files
for font in fm.findSystemFonts(fontPath):
    fm.fontManager.addfont(font)

import fun_convert_unit as fcu
import fun_process_data as fpd
import fun_plot_tools as fpt
import fun_robustness as fr
import fun_special_plot as fsp
import region_library as rlib

def plot_single_slice_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 1 panel slice globe'''
    # Choose the requested scenario
    scnDict = fpt.make_scenario_dict(rlzList, setDict)
    if setDict["plotPanel"] == 'RCP85':
        actData = scnDict["RCP8.5"]
    elif setDict["plotPanel"] == 'CESMS245':
        actData = scnDict['SSP2-4.5']
    elif setDict["plotPanel"] == 'GLENS':
        actData = scnDict['GLENS-SAI']
    elif setDict["plotPanel"] == 'CESMARISE15':
        actData = scnDict['ARISE-SAI-1.5']
    elif setDict["plotPanel"] == 'UKESMS245':
        actData = scnDict['UKESM-SSP2-4.5']
    elif setDict["plotPanel"] == 'UKESMARISE15':
        actData = scnDict['UKESM-ARISE-SAI-1.5']
    else:
        sys.exit('Check plotPanel input')

    # Obtain ens mean, pointwise max/min, or individual realization
    if setDict["plotEnsType"] == 'mean': #Ensemble mean
        panel = actData.mean(dim='realization')
    elif setDict["plotEnsType"] == 'max': #Pointwise max
        panel = actData.max(dim='realization')
    elif setDict["plotEnsType"] == 'min': #Pointwise min
        panel = actData.min(dim='realization')
    elif isinstance(setDict["plotEnsType"], int): #Single member
        panel = actData.isel(realization=setDict["plotEnsType"])
    else:
        sys.exit('Check plotEnsType input!')
        
    panelStr = setDict["plotPanel"]

    CL = 0.

    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    fig = plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance if setDict["cmap"] is None else setDict["cmap"]
    cbAuto = np.sort(
        [-np.nanquantile(panel.data,0.75), np.nanquantile(panel.data,0.75)])
    cbVals = cbAuto if setDict["cbVals"] is None else setDict["cbVals"]
    md = fpd.meta_book(setDict, dataDict, panel)
    lats = rlzList[0].lat
    lons = rlzList[0].lon
    plt.rcParams.update({'font.family': 'Open Sans'})
    plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
    fpt.drawOnGlobe(
        ax, panel, lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1],
        cbarBool=False, fastBool=True, extent='both',
        addCyclicPoint=setDict["addCyclicPoint"], alph=1)

    savePrfx = '' #Easy modification for unique filename
    saveStr = panelStr + '_' + md["varSve"] + '_' \
        + 'rlz' + str(setDict["plotEnsType"])
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

    return None
