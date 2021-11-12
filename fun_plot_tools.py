''' fun_plot_tools
Contains functions for plotting, e.g. drawing data on a globe, making kernel
density estimates, histograms, or step plots. Also includes functions related to
plotting, such as functions to choose colors or generate labels.

Unless otherwise specified:
Written by Daniel Hueholt
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
    data_cyc, lons_cyc = add_cyclic_point(data, coord=lons) #fixes white line by adding point

    ax.set_global()
    ax.coastlines(linewidth = 1.2, color='black')
    if(fastBool):
        image = ax.pcolormesh(lons_cyc, lats, data_cyc, transform=data_crs, cmap=cmap)
        # lonCirc = np.arange(0,360)
        # latCirc = np.zeros(np.shape(lonCirc)) + 75
        # plt.plot(lonCirc, latCirc, color='r', linewidth=1, transform=data_crs)
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
            # ic(delta_coord - delta_coord[0], delta_coord < 1, coord) #troubleshooting
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
    slicedData = data[tuple(slicer)] #DMH: assigned to var for easy access
    # DMH: manually assign ocean data (otherwise will be NaNs and output fails)
    if np.isnan(slicedData).all().data:
        sliceShape = np.shape(slicedData)
        merData = data.sel(lon=358.75).data
        slicedData = np.zeros(sliceShape)
        for sd,sv in enumerate(slicedData):
            slicedData[sd,0] = merData[sd]
        # ic(slicedData)
    new_data = ma.concatenate((data, slicedData), axis=axis) #DMH
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
        ax = sn.kdeplot(data=handles[0], label=labels[0], color=colors[0], linewidth=2) #Force index 0 or seaborn will be confused by 1-element "array" input
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
        ax = sn.histplot(data=handles[0], label=labels[0], color=colors[0], edgecolor='#3B3B3B', stat='probability', linewidth=0.8, kde=False, binwidth=binwidth) #Force index 0 or seaborn will be confused by 1-element "array" input
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
        ax = sn.histplot(data=handles[0], label=labels[0], color=colors[0], edgecolor='#3B3B3B', stat='density', linewidth=0.8, kde=False, binwidth=binwidth) #Force index 0 or seaborn will be confused by 1-element "array" input
        ax.set(xlabel=unit, ylabel='Density')
    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])

    ic(savename)
    plt.savefig(savename, dpi=dpiVal, bbox_inches='tight')
    plt.close()

def paint_by_numbers(colorsToPlot, scnId, nColors):
    ''' Choose colors for perceptual distinction given number of plotted objects
        GLENS Cntrl/Fdbck luminances: 80.484,67.169,55.274,49.081,46.343,35.522,23.790,11.597
        ARISE/SSP2-4.5 luminances: 90,85,80,75,70,65,60,55
    '''
    if scnId == 'RCP8.5':
        colList = ["#F2BABA", "#E88989", "#DF5757", "#D93636", "#D32828", "#A21F1F", "#701515", "#3F0C0C"]
    elif scnId == 'G1.2(8.5)':
        colList = ["#D2BBE8", "#B48FDA", "#9763CB", "#8346C1", "#7A3DB6", "#5C2E8A", "#3F1F5E", "#211132"]
    elif scnId == 'SSP2-4.5':
        colList = ["#FFCF65", "#FFC158", "#FFB34A", "#F8A53D", "#E8982E", "#D98B1F", "#CA7E0C", "#BB7100"]
    elif scnId == 'G1.5(2-4.5)':
        colList = ["#5BFCDC", "#48EDCE", "#33DFC0", "#12D0B2", "#00C2A5", "#00B498", "#00A68B", "#00997E"]
    else:
        ic('Unknown scenario!')
        colList = None

    if nColors == 0:
        pass
    elif nColors == 1:
        colorsToPlot.append(colList[3])
    elif nColors == 2:
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[4])
    elif nColors == 3:
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[5])
    elif nColors == 4:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[7])
    elif nColors == 5:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[6])
        colorsToPlot.append(colList[7])
    elif nColors == 6:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[6])
        colorsToPlot.append(colList[7])
    elif nColors == 7:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[6])
        colorsToPlot.append(colList[7])
    elif nColors == 8:
        colorsToPlot.append(colList[0])
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[4])
        colorsToPlot.append(colList[5])
        colorsToPlot.append(colList[6])
        colorsToPlot.append(colList[7])

    return colorsToPlot

def write_labels(labelsList, dataScnPoiKey, scnPoiKey, scnId, ensPrpKey, setDict):
    ''' Writes the label text given loop information '''
    for cdc,cdv in enumerate(scnPoiKey):
        if dataScnPoiKey is None:
            continue #No label for no data
        if cdv < 2020:
            continue #Do not auto-generate if interval starts during the 2010-2019 reference period
        startYearStr = str(cdv)
        endYearStr = str(cdv + setDict["timePeriod"] - 1)
        if cdv+setDict["timePeriod"] == 2100:
            endYearStr = str(2099)
        if (setDict["realization"] == 'mean') and cdv+setDict["timePeriod"]>ensPrp["dscntntyYrs"][0]:
            labelStr = startYearStr + '-' + endYearStr + ' ' + scnId + ' ' + '[r'+str(ensPrpKey[1])+']'
        elif (setDict["realization"] == 'mean') and cdv+setDict["timePeriod"]<ensPrp["dscntntyYrs"][0]:
            labelStr = startYearStr + '-' + endYearStr + ' ' + scnId + ' ' + '[r'+str(ensPrpKey[0])+']'
        else:
            labelStr = startYearStr + '-' + endYearStr + ' ' + scnId
        labelsList.append(labelStr)

    return labelsList

def generate_labels_colors(labelsList, colorsToPlot, dataDict, setDict, ensPrp, rfrncFlag, nGlensCntrlPoi, nGlensFdbckPoi, nAriseFdbckPoi, nAriseCntrlPoi):
    ''' Generate labels and colors for figure titles and output filenames.
    Clumsy to require so many inputs--the code is tidy but the logic is not!
    '''
    if rfrncFlag: #For the reference period beginning pre-2020
        if (setDict["cntrlPoi"][0] < 2020) and (dataDict["idGlensCntrl"] is not None): #Use GLENS if present as n(rlz) larger
            startRef = setDict["cntrlPoi"][0]
            endRef = startRef + setDict["timePeriod"] - 1
            # endRef = startRef + 5 - 1
            nGlensCntrlPoi = nGlensCntrlPoi - 1
            labelsList.append(str(startRef) + '-' + str(endRef) + ' Ref. RCP8.5')
        elif (setDict["s245CntrlPoi"][0] < 2020) and (dataDict["idS245Cntrl"] is not None): #Try ARISE if GLENS not available
            startRef = setDict["s245CntrlPoi"][0]
            endRef = startRef + setDict["timePeriod"] - 1
            # endRef = startRef + 5 - 1
            nAriseCntrlPoi = nAriseCntrlPoi - 1
            labelsList.append(str(startRef) + '-' + str(endRef) + ' Ref. SSP2-4.5')
        else:
            raise NoDataError('No data from reference period! Check inputs and try again.')
        rfrncColor = '#788697'
        colorsToPlot.append(rfrncColor)

    dataScnPoiKey = dataDict["idGlensCntrl"]
    scnPoiKey = setDict["cntrlPoi"]
    scnId = 'RCP8.5'
    ensPrpKey = ensPrp["drc"]
    labelsList = write_labels(labelsList, dataScnPoiKey, scnPoiKey, scnId, ensPrpKey, setDict)
    colorsToPlot = paint_by_numbers(colorsToPlot, scnId, nGlensCntrlPoi)

    dataScnPoiKey = dataDict["idGlensFdbck"]
    scnPoiKey = setDict["fdbckPoi"]
    scnId = 'G1.2(8.5)'
    ensPrpKey = ensPrp["drf"]
    labelsList = write_labels(labelsList, dataScnPoiKey, scnPoiKey, scnId, ensPrpKey, setDict)
    colorsToPlot = paint_by_numbers(colorsToPlot, scnId, nGlensFdbckPoi)

    dataScnPoiKey = dataDict["idS245Cntrl"]
    scnPoiKey = setDict["s245CntrlPoi"]
    scnId = 'SSP2-4.5'
    ensPrpKey = ensPrp["drs245"]
    labelsList = write_labels(labelsList, dataScnPoiKey, scnPoiKey, scnId, ensPrpKey, setDict)
    colorsToPlot = paint_by_numbers(colorsToPlot, scnId, nAriseCntrlPoi)

    dataScnPoiKey = dataDict["idArise"]
    scnPoiKey = setDict["arisePoi"]
    scnId = 'G1.5(2-4.5)'
    ensPrpKey = ensPrp["drari"]
    labelsList = write_labels(labelsList, dataScnPoiKey, scnPoiKey, scnId, ensPrpKey, setDict)
    colorsToPlot = paint_by_numbers(colorsToPlot, scnId, nAriseFdbckPoi)

    return labelsList, colorsToPlot

def line_from_scenario(scn, md):
    ''' Get line color and label from scenario information '''
    if 'GLENS:Control' in scn:
        activeColor = '#D93636' #Red
        activeLabel = md['cntrlStr']
    elif 'GLENS:Feedback' in scn:
        activeColor = '#8346C1' #Purple
        activeLabel = md['fdbckStr']
    elif 'ARISE:Feedback' in scn:
        activeColor = '#12D0B2' #Turquoise
        activeLabel = md['ariseStr']
    elif 'ARISE:Control' in scn:
        activeColor = '#F8A53D' #Orange
        activeLabel = md['s245Cntrl']
    else:
        activeColor = '#000000'
        activeLabel = 'Unknown'
        ic('Unknown scenario! Plotting with black line and unknown label.')

    return activeColor, activeLabel

def plot_metaobjects(scnToPlot, fig):
    ''' Determines which metaobjects to plot based on scenario '''
    b,t = plt.ylim()
    if any('ARISE:Control' in scn for scn in scnToPlot): #Triangle for change in number of rlzs
        plt.plot(2015, b+(abs(b-t))*0.01, color='#F8A53D', marker='v')
    if any('GLENS:Control' in scn for scn in scnToPlot):
        plt.plot(2030, b+(abs(b-t))*0.01, color='#D93636', marker='v')
    if any('GLENS:Feedback' in scn for scn in scnToPlot): #Dashed line for model SAI initiation
        plt.plot([2020,2020], [b,t], color='#8346C1', linewidth=0.7, linestyle='dashed')
    if any('ARISE:Feedback' in scn for scn in scnToPlot):
        plt.plot([2035,2035], [b,t], color='#12D0B2', linewidth=0.7, linestyle='dashed')

    return

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
