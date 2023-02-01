''' fun_zonal_plot
Plot zonal average
Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import numpy as np
import xarray as xr
xr.set_options(keep_attrs=True)
import cftime
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

import fun_process_data as fpd
import fun_plot_tools as fpt
import region_library as rlib

def plot_zonal_avg(rlzList, dataDict, setDict, outDict):
    ''' Plot zonal average of a variable '''
    toiStart, toiEnd = fpt.make_panels(rlzList, setDict)
    diffToiR85 = toiEnd['RCP8.5'] - toiStart['RCP8.5']
    diffToiS245 = toiEnd['SSP2-4.5'] - toiStart['SSP2-4.5']
    diffToiG12R85 = toiEnd['GLENS-SAI'] - toiStart['RCP8.5']
    diffToiG15S245 = toiEnd['ARISE-SAI-1.5'] - toiStart['SSP2-4.5']
    intiG12R85 = toiEnd['GLENS-SAI'] - toiEnd['RCP8.5']
    intiG15S245 = toiEnd['ARISE-SAI-1.5'] - toiEnd['SSP2-4.5']

    panels = (toiStart['RCP8.5'], toiStart['SSP2-4.5'], toiEnd['RCP8.5'], toiEnd['SSP2-4.5'])

    zonPanels = list()
    for pan in panels:
        panZon = pan.mean(dim='lon')
        zonPanels.append(panZon)

    plt.rcParams.update({'font.size': 9})
    plt.rcParams.update({'font.family': 'Lato'})
    ax = plt.subplot(2,2,1)
    ic('Panel 1')
    ic(zonPanels[0])
    plt.title(zonPanels[0].scenario)
    plt.plot(zonPanels[0].lat,zonPanels[0])
    plt.xlim(0,12)
    # plt.ylim(-10,10)
    # plt.xticks(np.arange(-12,14,2))

    ax = plt.subplot(2,2,2)
    ic('Panel 2')
    ic(zonPanels[1])
    plt.title(zonPanels[1].scenario)
    plt.plot(zonPanels[1].lat,zonPanels[1])
    plt.xlim(0,12)
    # plt.ylim(-10,10)
    # plt.xticks(np.arange(-12,14,2))

    ax = plt.subplot(2,2,3)
    ic('Panel 3')
    ic(zonPanels[2])
    plt.title(zonPanels[2].scenario)
    plt.plot(zonPanels[2].lat,zonPanels[2])
    plt.xlim(0,12)
    # plt.ylim(-10,10)
    # plt.xticks(np.arange(-12,14,2))

    ax = plt.subplot(2,2,4)
    ic('Panel 4')
    ic(zonPanels[3])
    plt.title(zonPanels[3].scenario)
    plt.plot(zonPanels[3].lat,zonPanels[3])
    plt.xlim(0,12)
    # plt.ylim(-10,10)
    # plt.xticks(np.arange(-12,14,2))

    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    savePrfx = 'ZON_'
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['frstDcd'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g4p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)
