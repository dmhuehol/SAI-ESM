''' plotting_tools
Contains functions to plot GLENS data, e.g. drawing data on a globe, making kernel
density estimates, histograms, step plots, etc. Includes functions related to
plotting, as well, such as paint_by_numbers which chooses colors for perceptual
difference based on the number of objects being drawn.

Unless otherwise specified:
Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
drawOnGlobe written by Prof. Elizabeth Barnes at Colorado State University
    Lightly edited by Daniel Hueholt
add_cyclic_point copied from cartopy utils by Prof. Elizabeth Barnes at Colorado State University
    Edited by Daniel Hueholt
'''

from icecream import ic
import sys
import warnings

import cartopy as ct
import cartopy.crs as ccrs
import cmocean as cmocean
import numpy as np
import numpy.ma as ma
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import seaborn as sn

def drawOnGlobe(ax, data, lats, lons, cmap='coolwarm', vmin=None, vmax=None, inc=None, cbarBool=True, contourMap=[], contourVals = [], fastBool=False, extent='both'):
    ''' Draws geolocated data on a globe '''

    data_crs = ct.crs.PlateCarree()
    data_cyc, lons_cyc = add_cyclic_point(data, coord=lons) #fixes white line by adding point#data,lons#ct.util.add_cyclic_point(data, coord=lons) #fixes white line by adding point

    ax.set_global()
    ax.coastlines(linewidth = 1.2, color='black')
    if(fastBool):
        image = ax.pcolormesh(lons_cyc, lats, data_cyc, transform=data_crs, cmap=cmap)
    else:
        image = ax.pcolor(lons_cyc, lats, data_cyc, transform=data_crs, cmap=cmap)

    if(np.size(contourMap) !=0 ):
        contourMap_cyc, __ = add_cyclic_point(contourMap, coord=lons) #fixes white line by adding point
        ax.contour(lons_cyc,lats,contourMap_cyc,contourVals, transform=data_crs, colors='fuchsia')

    if(cbarBool):
        cb = plt.colorbar(image, shrink=.75, orientation="vertical", pad=.02, extend=extent)
        cb.ax.tick_params(labelsize=6)
        try:
            cb.set_label(data.attrs['units'],size='small')
        except:
            print('No units in attributes; colorbar will be unlabeled.')
    else:
        cb = None

    image.set_clim(vmin,vmax)

    return cb, image

def add_cyclic_point(data, coord=None, axis=-1):
    ''' EAB: had issues with cartopy finding utils so copied for myself '''

    reverseSlicerBool = False #DMH
    if coord is not None:
        if coord.ndim != 1:
            raise ValueError('The coordinate must be 1-dimensional.')
        if len(coord) != data.shape[axis]:
            raise ValueError('The length of the coordinate does not match '
                             'the size of the corresponding dimension of '
                             'the data array: len(coord) = {}, '
                             'data.shape[{}] = {}.'.format(
                                 len(coord), axis, data.shape[axis]))
        delta_coord = np.diff(coord) #DMH: calculate grid spacing, essentially
        if not np.allclose(delta_coord, delta_coord[0]): #DMH: if grid spacing is not nearly uniform
            # ic(delta_coord - delta_coord[0])
            # ic(delta_coord < 1)
            # ic(coord)
            # sys.exit('STOP')
            warnings.warn('The coordinate is not equally spaced. This could be '
                          'because multiple sub-regions making up a single '
                          'region are being plotted (as when a region crosses '
                          'a meridian), in which case this message can be '
                          'ignored. Or, the underlying grid may be bad, in '
                          'which case that is problematic. Check your data and '
                          'be sure which applies to you!') #DMH
            reverseSlicerBool = True #DMH

        new_coord = ma.concatenate((coord, coord[-1:] + delta_coord[0]))
    slicer = [slice(None)] * data.ndim
    try:
        if not reverseSlicerBool: #DMH
            slicer[axis] = slice(0, 1) #DMH: Default behavior
        else: #DMH
            slicer[axis] = slice(1, 0) #DMH
            ic('Slicer has been reversed') #DMH
    except IndexError:
        raise ValueError('The specified axis does not correspond to an '
                         'array dimension.')
    new_data = ma.concatenate((data, data[tuple(slicer)]), axis=axis)
    if coord is None:
        return_value = new_data
    else:
        return_value = new_data, new_coord
    return return_value

def plot_pdf_kdeplot(handles, colors, labels, unit, savename, dpiVal=400):
    ''' Make kernel density estimates for several input handles '''

    plt.figure()
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            ax = sn.kdeplot(data=h, label=labels[ind], color=colors[ind], linewidth=2)
            ax.set(xlabel=unit, ylabel='Density')
    else:
        ax = sn.kdeplot(data=handles, label=labels, color=colors, linewidth=2)
        ax.set(xlabel=unit, ylabel='Density')

    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])

    ic(savename)
    plt.savefig(savename, dpi=dpiVal, bbox_inches='tight')
    plt.close()

def plot_pdf_hist(handles, colors, labels, unit, savename, binwidth, dpiVal=400):
    ''' Make histograms for several input handles'''

    plt.figure()
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            ax = sn.histplot(data=h, label=labels[ind], color=colors[ind], edgecolor='#3B3B3B', stat='probability', linewidth=0.8, kde=False, binwidth=binwidth)
            ax.set(xlabel=unit, ylabel='Density')
    else:
        ax = sn.histplot(data=handles, label=labels, color=colors, edgecolor='#3B3B3B', stat='probability', linewidth=0.8, kde=False, binwidth=binwidth)
        ax.set(xlabel=unit, ylabel='Density')

    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])

    ic(savename)
    plt.savefig(savename, dpi=dpiVal, bbox_inches='tight')
    plt.close()

def plot_pdf_step(handles, colors, labels, unit, savename, binwidth, dpiVal=400):
    ''' Make step plots for several input handles '''

    plt.figure()
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            bins = np.arange(np.nanmin(h),np.nanmax(h),binwidth)
            ax = plt.hist(h, label=labels[ind], color=colors[ind], bins=bins, density=True, histtype='step')
            plt.xlabel(unit)
            plt.ylabel("Density")
    else:
        ax = sn.histplot(data=handles, label=labels, color=colors, edgecolor='#3B3B3B', stat='density', linewidth=0.8, kde=False, binwidth=binwidth)
        ax.set(xlabel=unit, ylabel='Density')

    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])

    ic(savename)
    plt.savefig(savename, dpi=dpiVal, bbox_inches='tight')
    plt.close()

def select_colors(baselineFlag, nCntrl, nFdbck):
    '''Returns colors for a set number of feedback and control objects '''

    colorsToPlot = list()
    baselineColor = 'slategrey'
    cntrlColors = ["#F2BABA", "#E88989", "#DF5757", "#D93636", "#D32828", "#A21F1F", "#701515", "#3F0C0C"]
    fdbckColors = ["#D2BBE8", "#B48FDA", "#9763CB", "#8346C1", "#7A3DB6", "#5C2E8A", "#3F1F5E", "#211132"]

    if baselineFlag:
        colorsToPlot.append(baselineColor)

    colorsToPlot = paint_by_numbers(colorsToPlot, cntrlColors, nCntrl)
    colorsToPlot = paint_by_numbers(colorsToPlot, fdbckColors, nFdbck)
    return colorsToPlot

def paint_by_numbers(colorsToPlot, colList, nfc):
    ''' Choose colors for perceptual distinction given number of objects to be plotted '''

    if nfc == 0:
        pass
    elif nfc == 1:
        colorsToPlot.append(colList[3])
    elif nfc == 2:
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[4])
    elif nfc == 3:
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[5])
    elif nfc == 4:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[7])
    elif nfc == 5:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[6])
        colorsToPlot.append(colList[7])
    elif nfc == 6:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[6])
        colorsToPlot.append(colList[7])
    elif nfc == 7:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[6])
        colorsToPlot.append(colList[7])
    elif nfc == 8:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[4])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[6])
        colorsToPlot.append(colList[7])

    return colorsToPlot

def generate_labels(labelsList, setDict, ensPrp, baselineFlag):
    ''' Generate labels for figure titles and output filenames '''

    for cdc,cdv in enumerate(setDict["cntrlPoi"]):
        if cdv < 2020:
            continue #Do not auto-generate if interval starts during the 2010-2019 "Baseline" period
        startYearStr = str(cdv)
        endYearStr = str(cdv + setDict["timePeriod"] - 1)
        if cdv+setDict["timePeriod"] == 2100:
            endYearStr = str(2099)
        if (setDict["realization"] == 'mean') and cdv+setDict["timePeriod"]>ensPrp["dscntntyYrs"][0]:
            labelStr = startYearStr + '-' + endYearStr + ' ' + 'RCP8.5' + ' ' + '[r'+str(ensPrp["drc"][1])+']'
        elif (setDict["realization"] == 'mean') and cdv+setDict["timePeriod"]<ensPrp["dscntntyYrs"][0]:
            labelStr = startYearStr + '-' + endYearStr + ' ' + 'RCP8.5' + ' ' + '[r'+str(ensPrp["drc"][0])+']'
        else:
            labelStr = startYearStr + '-' + endYearStr + ' ' + 'RCP8.5'
        labelsList.append(labelStr)

    for cdc,cdv in enumerate(setDict["fdbckPoi"]):
        if cdv < 2020:
            continue #Do not auto-generate if interval starts during the 2010-2019 "Baseline" period
        startYearStr = str(cdv)
        endYearStr = str(cdv + setDict["timePeriod"] - 1)
        if cdv+setDict["timePeriod"] == 2100:
            endYearStr = str(2099)
        if (setDict["realization"] == 'mean') and cdv+setDict["timePeriod"]>ensPrp["dscntntyYrs"][0]:
            labelStr = startYearStr + '-' + endYearStr + ' ' + 'SAI' + ' ' + '[r'+str(ensPrp["drf"][1])+']'
        elif (setDict["realization"] == 'mean') and cdv+setDict["timePeriod"]<ensPrp["dscntntyYrs"][0]:
            labelStr = startYearStr + '-' + endYearStr + ' ' + 'SAI' + ' ' + '[r'+str(ensPrp["drf"][0])+']'
        else:
            labelStr = startYearStr + '-' + endYearStr + ' ' + 'SAI'
        labelsList.append(labelStr)

    if baselineFlag:
        labelsList.insert(0,'2011-2030 Baseline')
        if (setDict["realization"] == 'mean'):
            labelsList[0] = labelsList[0] + ' ' + '[r21]'

    return labelsList

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
    cbarRange = np.array([cbarDict["range"]])

    if cbarDict["direction"] == 'horizontal':
        plt.figure(figsize=(9,3))
        img = plt.imshow(cbarRange, cmap=cbarDict["cmap"])
        plt.gca().set_visible(False)
        colorAx = plt.axes([0.1,0.2,0.8,0.6])
        cb = plt.colorbar(orientation='horizontal', cax=colorAx)
        for label in cb.ax.get_xticklabels():
            print(label)
            # label.set_fontproperties(FiraSansThin) #Set font
            # label.set_fontsize(18) #Set font size
    elif cbarDict["direction"] == 'vertical':
        plt.figure(figsize=(3,9))
        img = plt.imshow(cbarRange, cmap=cbarDict["cmap"])
        plt.gca().set_visible(False)
        colorAx = plt.axes([0.1,0.2,0.5,0.6])
        cb = plt.colorbar(orientation='vertical', cax=colorAx)
        for label in cb.ax.get_yticklabels():
            print(label)
            # label.set_fontproperties(FiraSansThin)
            # label.set_fontsize(18)
    else:
        sys.exit('Direction must be either horizontal or vertical.')

    cb.set_label(cbarDict["label"], size='large')
    # cb.set_label(cbarDict["label"], size='large', fontproperties=FiraSansThin) # Set font
    plt.savefig(savePath + saveName + '.png', dpi=dpiVal)

def find_widest_quantile(darr):
    ''' Figures out whether [0.01Q,-0.01Q] or [-0.99Q,0.99Q] is more appropriate '''
    absQuant99 = np.abs(darr.quantile(0.99))
    absQuant001 = np.abs(darr.quantile(0.01))
    bigVal = np.max([absQuant99,absQuant001])
    widestRange = np.array([-bigVal, bigVal])

    return widestRange
