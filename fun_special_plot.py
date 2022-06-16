''' fun_special_plot
Contains some special plotting functions for specific uses, e.g. plotting
colorbars alone or plotting quantiles vs. members for robustness.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys
import warnings

import matplotlib.font_manager as fm
fontPath = '/Users/dhueholt/Library/Fonts/'  #Location of font files
for font in fm.findSystemFonts(fontPath):
    fm.fontManager.addfont(font)

import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

def save_colorbar(cbarDict, savePath, saveName, dpiVal=400):
    ''' Plots and saves a colorbar
        cbarDict should have a cmap, range, direction (horizontal vs. vertical),
        and label. Ex:
        cbarDict = {
            "cmap": cmocean.cm.delta,
            "range": [-15,15],
            "direction": 'vertical',
            "label": 'percent'
        }
    '''
    plt.rcParams.update({'font.size': 0})
    plt.rcParams.update({'font.family': 'Lato'})

    cbarRange = np.array([cbarDict["range"]])

    if cbarDict["direction"] == 'horizontal':
        plt.figure(figsize=(9,2.5))
        img = plt.imshow(cbarRange, cmap=cbarDict["cmap"])
        plt.gca().set_visible(False)
        colorAx = plt.axes([0.1,0.2,0.8,0.3])
        cb = plt.colorbar(orientation='horizontal', cax=colorAx)
        for label in cb.ax.get_xticklabels():
            print(label)
            # label.set_fontproperties(FiraSansThin) #Set font
            label.set_fontsize(0) #Set font size
    elif cbarDict["direction"] == 'vertical':
        plt.figure(figsize=(2.5,9))
        img = plt.imshow(cbarRange, cmap=cbarDict["cmap"])
        plt.gca().set_visible(False)
        colorAx = plt.axes([0.1,0.2,0.2,0.6])
        cb = plt.colorbar(orientation='vertical', cax=colorAx)
        for label in cb.ax.get_yticklabels():
            print(label)
            # label.set_fontproperties(FiraSansThin)
            # label.set_fontsize(18)
    else:
        sys.exit('Direction must be either horizontal or vertical.')

    # cb.set_label(cbarDict["label"], size='large')
    # cb.set_label(cbarDict["label"], size='large', fontproperties=FiraSansThin) # Set font
    plt.savefig(savePath + saveName + '.png', dpi=dpiVal)

def quantiles_vs_members(robustness, savePath=None):
    ''' Plots quantiles vs number of members '''
    plt.rcParams.update({'font.family': 'Palanquin'})
    plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight

    qList = list()
    quantiles = np.linspace(0,1,500)
    for q in quantiles:
        qList.append(np.nanquantile(robustness,q))

    plt.figure(figsize=(9,6))
    plt.plot(quantiles, qList, color='#ff4d81', linewidth=2.5)

    plt.xlabel('Quantiles', fontsize=18, fontweight='light')
    plt.ylabel('Ensemble members', fontsize=18, fontweight='light')
    plt.xlim(-0.01,1.01)
    plt.xticks([0,0.2,0.4,0.6,0.8,1], fontsize=14)
    plt.ylim(0,22)
    plt.yticks([0,3,6,9,12,15,18,21], fontsize=14)
    plt.title('Quantiles vs members: annual mean 2m temp', fontsize=22, fontweight='light')

    if savePath == None:
        plt.show()
    else:
        plt.savefig(savePath + 'QuantileVsMembers_AnnualMeanSeaSurfaceTemp' + '.png')
