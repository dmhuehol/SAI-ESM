# # Emphasize matching members in spaghetti timeseries
# def plot_ens_spaghetti_timeseries(darrList, dataDict, setDict, outDict):
#     ''' Make a simple timeseries of output variable. Ensemble members are
#     visualized in a familiar, basic spaghetti plot. '''
#     plotRlzMn = True
#     setYear = [2015, 2070]
#     timeSlice = slice(cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0), cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
#
#     # Plot timeseries
#     fig,ax = plt.subplots()
#     scnToPlot = list() #Make list of all scenarios to be plotted
#     firstYearEachRlz = list()
#     for scnDarr in darrList:
#         rlzInScn = scnDarr['realization'].data #Number of realizations in scenario
#         scnToPlot.append(scnDarr.scenario) #Add scenario to list
#         if setDict["areaAvgBool"] == True:
#             rlzMn = scnDarr[len(dataToPlot['realization'])-1] #Last member is ensemble mean
#         elif setDict["areaAvgBool"] == 'sum':
#             rlzMn = scnDarr.copy()
#         for rc in rlzInScn:
#             rlzToi = scnDarr.sel(realization=rc, time=timeSlice) #Single rlz at time of interest
#             rlzLoi = fpd.obtain_levels(rlzToi, setDict["levOfInt"]) #Level of interest
#             rlzToPlot, locStr, locTitleStr = fpd.manage_area(rlzLoi, setDict["regOfInt"], areaAvgBool=setDict["areaAvgBool"]) #Area of interest
#             md = fpd.meta_book(setDict, dataDict, rlzToPlot, labelsToPlot=None) #Extract metadata
#             activeColor, activeLabel = fpt.line_from_scenario(rlzToPlot.scenario, md)
#             yrsToPlot = rlzToPlot['time'].dt.year.data #bndDct['mtchYrs']
#             firstYearEachRlz.append(rlzToPlot.data[0])
#             ic(len(rlzInScn))
#             if rc == 4:
#                 n = rc
#                 if "ARISE:Control" in rlzToPlot.scenario:
#                     actLab = 'SSP2-4.5 rlz ' + str(rc+1)
#                     plt.plot(yrsToPlot, rlzToPlot, color='#7FFF00', label=actLab, linewidth=0.5)
#                 elif "ARISE:Feedback" in rlzToPlot.scenario:
#                     actLab = 'ARISE rlz ' + str(rc+1)
#                     plt.plot(yrsToPlot, rlzToPlot, color='#264D00', label=actLab, linewidth=1.2)
#             elif rc == 9:
#                 n2 = rc
#                 if "ARISE:Feedback" in rlzToPlot.scenario:
#                     actLab = 'ARISE rlz ' + str(rc+1)
#                     plt.plot(yrsToPlot, rlzToPlot, color='#FF0081', label=actLab, linewidth=1.2)
#                 else:
#                     ic()
#                     # plt.plot(yrsToPlot, rlzToPlot, color='#D3D3D3', linewidth=0.3) #Individual rlz
#             else:
#                 ic()
#                 # plt.plot(yrsToPlot, rlzToPlot, color='#D3D3D3', linewidth=0.3) #Individual rlz
#         if plotRlzMn:
#             rlzToiMn = rlzMn.sel(time=timeSlice)
#             rlzLoiMn = fpd.obtain_levels(rlzToiMn, setDict["levOfInt"])
#             rlzAoiMn, _, __ = fpd.manage_area(rlzLoiMn, setDict["regOfInt"], areaAvgBool=setDict["areaAvgBool"])
#             rlzToPlotMn = rlzAoiMn.mean(dim='realization')
#             # plt.plot(yrsToPlot, rlzToPlotMn, color='#D3D3D3', label=activeLabel, linewidth=1.5) #Ens mean
#
#     # Plot metadata and settings
#     b,t = plt.ylim() if setDict['ylim'] is None else setDict['ylim']
#     fpt.plot_metaobjects(scnToPlot, fig, b, t)
#     leg = plt.legend()
#     plt.ylabel(md['unit'])
#     plt.autoscale(enable=True, axis='x', tight=True)
#     plt.autoscale(enable=True, axis='y', tight=True)
#     plt.xlim(setYear[0], setYear[1])
#     plt.title(md['varStr'] + ' ' + md['levStr'] + ': ' + str(setYear[0]) + '-' + str(setYear[1]) + ' ' + locTitleStr  + ' ' + 'spaghetti')
#
#     # Save image
#     savePrfx = '6' + 'AF' + str(n+1) + 'AF' + str(n2+1) + 'AC' + str(n+1) + '_'
#     saveStr = md['varSve'] + '_' + md['levSve'] + '_' + str(setYear[0]) + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' + md['ensPid']['spg']
#     savename = outDict["savePath"] + savePrfx + saveStr + '.png'
#     plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
#     plt.close()
#     ic(savename)

from icecream import ic
import sys
import warnings

import matplotlib.font_manager as fm
fontPath = '/Users/dhueholt/Library/Fonts/'  #Location of font files
for font in fm.findSystemFonts(fontPath):
    fm.fontManager.addfont(font)

import cartopy as ct
import cartopy.crs as ccrs
import cmocean as cmocean
import numpy as np
import numpy.ma as ma
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import seaborn as sn

import fun_process_data as fpd

cbarDict = {
    "cmap": cmocean.cm.balance,
    "range": [-1.5,1.5],
    "direction": 'vertical',
    "label": 'percent'
}
savePath = '/Users/dhueholt/Documents/GLENS_fig/20220519_seaiceareaAndMute/'
saveName = 'balance.png'
import fun_plot_tools as fpt
fpt.save_colorbar(cbarDict, savePath, saveName, dpiVal=400)
