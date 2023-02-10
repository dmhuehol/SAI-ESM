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

    activeData = rlzList[0]

    globeEnsMn = activeData.mean(dim='realization')
    globeEnsMax = activeData.max(dim='realization')
    globeEnsMin = activeData.min(dim='realization')
    rlzInd = 4
    globeEnsMember = activeData.isel(realization=rlzInd)
    panel = globeEnsMember
    panelStr = setDict["plotPanel"]
    # Need intelligent scenario checking
    # if setDict["plotPanel"] == 'snapR85':
    #     snapR85 = toiEnd['RCP8.5'] - toiStart['RCP8.5']
    #     panel = snapR85
    # elif setDict["plotPanel"] == 'snapS245':
    #     snapS245 = toiEnd['SSP2-4.5'] - toiStart['SSP2-4.5']
    #     panel = snapS245
    # elif setDict["plotPanel"] == 'snapGLENS':
    #     snapGLENS = toiEnd['GLENS-SAI'] - toiStart['RCP8.5']
    #     panel = snapGLENS
    # elif setDict["plotPanel"] == 'snapARISE15':
    #     snapARISE15 = toiEnd['ARISE-SAI-1.5'] - toiStart['SSP2-4.5']
    #     panel = snapARISE15
    # elif setDict["plotPanel"] == 'intiGLENS':
    #     intiGLENS = toiEnd['GLENS-SAI'] - toiEnd['RCP8.5']
    #     panel = intiGLENS
    # elif setDict["plotPanel"] == 'intiARISE15':
    #     intiARISE15 = toiEnd['ARISE-SAI-1.5'] - toiEnd['SSP2-4.5']
    #     panel = intiARISE15
    # elif setDict["plotPanel"] == 'snapUKS245':
    #     snapUKS245 = toiEnd['UKESM-SSP2-4.5'] - toiStart['UKESM-SSP2-4.5']
    #     panel = snapUKS245
    # elif setDict["plotPanel"] == 'snapUKARISE15':
    #     snapUKS245 = toiEnd['UKESM-ARISE-SAI-1.5'] \
    #         - toiStart['UKESM-SSP2-4.5']
    #     panel = snapUKS245
    # elif setDict["plotPanel"] == 'intiUKARISE15':
    #     intiUKARISE15 = toiEnd['UKESM-ARISE-SAI-1.5'] \
    #         - toiEnd['UKESM-SSP2-4.5']
    #     panel = intiUKARISE15
    # else:
    #     blank = toiEnd['GLENS-SAI'].copy()
    #     blank.data = toiEnd['GLENS-SAI'] - toiEnd['GLENS-SAI']
    #     panel = blank

    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    fig = plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance if setDict["cmap"] is None else setDict["cmap"]
    cbAuto = np.sort(
        [-np.nanquantile(panel.data,0.75), np.nanquantile(panel.data,0.75)])
    cbVals = cbAuto if setDict["cbVals"] is None else setDict["cbVals"]
    # md = fpd.meta_book(setDict, dataDict, rlzList[0])
    lats = rlzList[0].lat
    lons = rlzList[0].lon
    plt.rcParams.update({'font.family': 'Open Sans'})
    plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
    fpt.drawOnGlobe(
        ax, panel, lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1],
        cbarBool=True, fastBool=True, extent='max',
        addCyclicPoint=setDict["addCyclicPoint"], alph=1)

    savePrfx = 'rlz' + str(rlzInd) + '_' #Easy modification for unique filename
    # Need automatic filename generation
    # if 'CESM2-ARISE' in panel.scenario:
    #     saveStr = panelStr + '_' + md['varSve'] + '_' + md['levSve'] \
    #         + '_' + md['lstDcd']["CESM2-ARISE"] + '_' + md['ensStr'] \
    #         + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    # elif 'GLENS' in panel.scenario:
    #     saveStr = panelStr + '_' + md['varSve'] + '_' + md['levSve'] \
    #         + '_' + md['lstDcd']["GLENS"] + '_' + md['ensStr'] + '_' \
    #         + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    # elif 'UKESM-ARISE' in panel.scenario:
    #     saveStr = panelStr + '_' + md['varSve'] + '_' + md['levSve'] \
    #         + '_' + md['lstDcd']["UKESM-ARISE"] + '_' + md['ensStr'] \
    #         + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    # else:
    #     ic('Unable to generate saveStr for scenario')
    #     saveStr = panelStr

    saveStr = panelStr

    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

    return None
