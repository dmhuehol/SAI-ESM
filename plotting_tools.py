''' plotting_tools
Contains functions to plot GLENS data, e.g. drawing data on a globe, making kernel
density estimates, histograms, step plots, etc. Includes functions related to
plotting, as well, such as paint_by_numbers which chooses colors for perceptual
difference based on the number of objects being drawn.

Unless otherwise specified:
Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import cartopy as ct
import cartopy.crs as ccrs
import cmocean as cmocean
import numpy as np
import numpy.ma as ma
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import seaborn as sn

# drawOnGlobe written by Prof. Elizabeth Barnes at Colorado State University, lightly edited by Daniel Hueholt
# add_cyclic_point copied from cartopy utils by Prof. Elizabeth Barnes at Colorado State University

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
    ''' had issues with cartopy finding utils so copied for myself -EAB '''

    if coord is not None:
        if coord.ndim != 1:
            raise ValueError('The coordinate must be 1-dimensional.')
        if len(coord) != data.shape[axis]:
            raise ValueError('The length of the coordinate does not match '
                             'the size of the corresponding dimension of '
                             'the data array: len(coord) = {}, '
                             'data.shape[{}] = {}.'.format(
                                 len(coord), axis, data.shape[axis]))
        delta_coord = np.diff(coord)
        if not np.allclose(delta_coord, delta_coord[0]):
            raise ValueError('The coordinate must be equally spaced.')
        new_coord = ma.concatenate((coord, coord[-1:] + delta_coord[0]))
    slicer = [slice(None)] * data.ndim
    try:
        slicer[axis] = slice(0, 1)
    except IndexError:
        raise ValueError('The specified axis does not correspond to an '
                         'array dimension.')
    new_data = ma.concatenate((data, data[tuple(slicer)]), axis=axis)
    if coord is None:
        return_value = new_data
    else:
        return_value = new_data, new_coord
    return return_value

def plot_pdf_kdeplot(handles, colors, labels, unit, saveName, dpiVal=400):
    ''' Make kernel density estimates for several input handles '''

    plt.figure()
    print('Plotting!')
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            ax = sn.kdeplot(data=h, label=labels[ind], color=colors[ind], linewidth=2)
            ax.set(xlabel=unit, ylabel='Density')
    else:
        ax = sn.kdeplot(data=handles, label=labels, color=colors, linewidth=2)
        ax.set(xlabel=unit, ylabel='Density')

    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])

    ic(saveName)
    plt.savefig(saveName + '.png', dpi=dpiVal, bbox_inches='tight')

def plot_pdf_hist(handles, colors, labels, unit, saveName, binwidth, dpiVal=400):
    ''' Make histograms for several input handles'''

    plt.figure()
    print('Plotting!')
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            ax = sn.histplot(data=h, label=labels[ind], color=colors[ind], edgecolor='#3B3B3B', stat='probability', linewidth=0.8, kde=False, binwidth=binwidth)
            ax.set(xlabel=unit, ylabel='Density')
    else:
        ax = sn.histplot(data=handles, label=labels, color=colors, edgecolor='#3B3B3B', stat='probability', linewidth=0.8, kde=False, binwidth=binwidth)
        ax.set(xlabel=unit, ylabel='Density')

    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])

    ic(saveName)
    plt.savefig(saveName + '.png', dpi=dpiVal, bbox_inches='tight')

def plot_pdf_step(handles, colors, labels, unit, saveName, binwidth, dpiVal=400):
    ''' Make step plots for several input handles '''

    plt.figure()
    print('Plotting!')
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            bins = np.arange(np.min(h),np.max(h),binwidth)
            ax = plt.hist(h, label=labels[ind], color=colors[ind], bins=bins, density=True, histtype='step')
            plt.xlabel(unit)
            plt.ylabel("Density")
    else:
        ax = sn.histplot(data=handles, label=labels, color=colors, edgecolor='#3B3B3B', stat='density', linewidth=0.8, kde=False, binwidth=binwidth)
        ax.set(xlabel=unit, ylabel='Density')

    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])

    ic(saveName)
    plt.savefig(saveName + '.png', dpi=dpiVal, bbox_inches='tight')

def select_colors(baselineFlag, nFdbck, nCntrl):
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

    if nfc == 1:
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

def generate_labels(labelsList, intervalsToPlot, timePeriod, type):
    ''' Generate labels for figure titles and output filenames '''

    for cdc,cdv in enumerate(intervalsToPlot):
        if cdv == 2010:
            continue
        startYearStr = str(cdv)
        endYearStr = str(cdv + timePeriod - 1)
        if cdv+timePeriod == 2100:
            endYearStr = str(2099)
        labelStr = startYearStr + '-' + endYearStr + ' ' + type
        labelsList.append(labelStr)

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
