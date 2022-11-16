''' fun_plot_tools
Contains plotting functions, e.g. drawing data on a globe. Also includes
functions for related tasks, such as generating labels.

Unless otherwise specified:
Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University

drawOnGlobe written by Prof. Elizabeth Barnes at Colorado State University
    Lightly edited by Daniel Hueholt
add_cyclic_point copied from cartopy utils by Prof. Elizabeth Barnes at Colorado State University
    Modified by Daniel Hueholt to add edge cases and documentation
'''

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

def make_panels(rlzList, setDict):
    ''' Extract periods of interest, average, & store by scenario for panels '''
    toiStart = dict()
    toiEnd = dict()
    for rc,rDarr in enumerate(rlzList):
        rlzLoi = fpd.obtain_levels(rDarr, setDict["levOfInt"])
        if 'realization' in rlzLoi.dims: # If 'ensplot' specified in wrap_basicplots_script
            mnInd = len(rlzLoi.realization)-1 # Final index is ensemble mean
            rlzLoi = rlzLoi.isel(realization=mnInd)
        shrtScn = rlzLoi.scenario.split('/')[len(rlzLoi.scenario.split('/'))-1]
        if 'Control' in rlzLoi.attrs['scenario']:
            if 'GLENS' in rlzLoi.attrs['scenario']:
                # ic('GLENS Control')
                toiStartLp = fpd.average_over_years(
                    rlzLoi, setDict["startIntvl"][0], setDict["startIntvl"][1])
                toiEndLp = fpd.average_over_years(
                    rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
                toiStart[shrtScn] = toiStartLp
                toiEnd[shrtScn] = toiEndLp
            elif 'ARISE' in rlzLoi.attrs['scenario']:
                # ic('ARISE Control')
                toiStartLp = fpd.average_over_years(
                    rlzLoi, setDict["startIntvl"][2], setDict["startIntvl"][3])
                toiEndLp = fpd.average_over_years(
                    rlzLoi, setDict["endIntvl"][2], setDict["endIntvl"][3])
                toiStart[shrtScn] = toiStartLp
                toiEnd[shrtScn] = toiEndLp
        elif 'Feedback' in rlzLoi.attrs['scenario']:
            if 'GLENS' in rlzLoi.attrs['scenario']:
                # ic('GLENS Feedback')
                toiEndLp = fpd.average_over_years(
                    rlzLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
                toiEnd[shrtScn] = toiEndLp
            elif 'ARISE' in rlzLoi.attrs['scenario']:
                # ic('ARISE Feedback')
                toiEndLp = fpd.average_over_years(
                    rlzLoi, setDict["endIntvl"][2], setDict["endIntvl"][3])
                toiEnd[shrtScn] = toiEndLp
        else:
            ic('This should not occur, but does it?')

    return toiStart, toiEnd

def drawOnGlobe(
        ax, data, lats, lons, cmap='coolwarm', vmin=None, vmax=None, inc=None,
        cbarBool=True, contourMap=[], contourVals = [], fastBool=False,
        extent='both', addCyclicPoint=False, alph=1):
    ''' Draws geolocated data on a globe. Written by Prof. Elizabeth Barnes at
        Colorado State University, lightly edited by Daniel Hueholt '''
    data_crs = ct.crs.PlateCarree()
    if addCyclicPoint: # Add cyclic point to prime meridian for ocean data
        data_cyc, lons_cyc = add_cyclic_point(data, coord=lons, axis=-1)
    else:
        data_cyc = data
        lons_cyc = lons
    ax.set_global()
    ax.coastlines(linewidth = 1.2, color='black')

    if(fastBool):
        image = ax.pcolormesh(
            lons_cyc, lats, data_cyc, transform=data_crs, cmap=cmap, alpha=alph)
        # lonCirc = np.arange(0,360)
        # latCirc = np.zeros(np.shape(lonCirc)) + 75
        # plt.plot(lonCirc, latCirc, color='r', linewidth=1, transform=data_crs)
    else:
        image = ax.pcolor(
            lons_cyc, lats, data_cyc, transform=data_crs, cmap=cmap)
    if(np.size(contourMap) !=0 ):
        if addCyclicPoint:
            contourMap_cyc, __ = add_cyclic_point(
                contourMap, coord=lons, axis=-1)
        else:
            contourMap_cyc = contourMap
        ax.contour(
            lons_cyc, lats, contourMap_cyc, contourVals, transform=data_crs,
            colors='fuchsia')
    if(cbarBool):
        cb = plt.colorbar(
            image, shrink=.75, orientation="vertical", pad=.02, extend=extent)
        cb.ax.tick_params(labelsize=6) #def: labelsize=6
        try:
            # cb.set_label(data.attrs['units'],size='small')
            cb.set_label('\u00B0C',size='medium')
            # cb.set_label('num ensemble members',size='small')
        except:
            print('No units in attributes; colorbar will be unlabeled.')
    else:
        cb = None
    image.set_clim(vmin,vmax)

    return cb, image

def add_cyclic_point(data, coord=None, axis=-1):
    ''' EAB: had issues with cartopy finding utils so copied for myself. Edited
    by Daniel Hueholt to deal with various edge cases. '''
    reverseSlicerBool = False #DMH: If you still have a white stripe at Prime Meridian, try flipping truth of this Bool

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
    # If plotting non-ocean data and the process fails with an obscure error,
    #   try commenting this block back out!
    if np.isnan(slicedData).all().data:
        sliceShape = np.shape(slicedData)
        merData = data.sel(lon=358.75).data #TODO: should work with numpy array
        slicedData = np.zeros(sliceShape)
        for sd,sv in enumerate(slicedData):
            slicedData[sd,0] = merData[sd]
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

def mute_by_numbers(thresh):
    ''' Mute a colorbar. Colors in comments are HSLuv unless otherwise noted'''
    # TODO: currently only works for GLENS (21 colors)
    grayList = ['#000000',
                '#111111',
                '#1b1b1b',
                '#262626',
                '#303030',
                '#3b3b3b',
                '#474747',
                '#525252',
                '#5e5e5e',
                '#6a6a6a',
                '#777777', #HSL=0,0,50
                '#848484',
                '#919191',
                '#9e9e9e',
                '#ababab',
                '#b9b9b9',
                '#c6c6c6',
                '#d4d4d4',
                '#e2e2e2',
                '#f1f1f1',
                '#ffffff']

    pinkList = ['#000000',
                '#2c000d',
                '#3f0016',
                '#52001e',
                '#660028',
                '#7b0031',
                '#90003b',
                '#a60045',
                '#bc004f',
                '#d3005a',
                '#ff1470',
                '#ea0064',
                '#ff4d81', #HSL=0,100,60
                '#ff6c91',
                '#ff86a1',
                '#ff9cb1',
                '#ffb1c0',
                '#ffc6d0',
                '#ffd9e0',
                '#ffecef',
                '#ffffff']

    ghostgrayList = ['#000000',
                     '#050505',
                     '#131313',
                     '#1d1d1d',
                     '#282828',
                     '#333333',
                     '#3c3c3c',
                     '#484848',
                     '#535353',
                     '#5e5e5e',
                     '#6a6a6a',
                     '#757575',
                     '#818181',
                     '#8e8e8e',
                     '#9a9a9a',
                     '#a8a8a8',
                     '#b6b6b6',
                     '#c5c5c5',
                     '#d2d2d2',
                     '#dfdfdf',
                     '#eeeeee']

    ghostlightList = ['#000000',
                      '#05050a',
                      '#12111f',
                      '#1f1933',
                      '#2f204a', #5
                      '#43265c',
                      '#532d61',
                      '#603b60', #8
                      '#6a4961',
                      '#725763', #10
                      '#7b6566', #11
                      '#827269', #12
                      '#89816c', #13
                      '#8f916f',
                      '#93a070',
                      '#97b171',
                      '#9bc26e', #17
                      '#a2d464',
                      '#bce058', #19
                      '#dbe955',
                      '#fef255'] #H,S,L=80.8,93.0,94.0

    grayList = ghostgrayList
    colorList = ghostlightList

    muteList = grayList[:thresh] + colorList[thresh:]

    return muteList

def mute_by_numbers_arise(thresh):
    ''' Mute a colorbar. Colors in comments are HSLuv unless otherwise noted'''

    ghostgrayList = ['#131313',
                     '#282828',
                     '#3c3c3c',
                     '#535353',
                     '#6a6a6a',
                     '#818181',
                     '#9a9a9a',
                     '#b6b6b6',
                     '#d2d2d2',
                     '#eeeeee']

    ghostlightList = ['#12111f',
                      '#2f204a', #5
                      '#532d61',#8
                      '#6a4961',#10
                      '#7b6566', #11#12
                      '#89816c', #13
                      '#93a070',
                      '#9bc26e', #17
                      '#bce058', #19
                      '#fef255'] #H,S,L=80.8,93.0,94.0

    grayList = ghostgrayList
    colorList = ghostlightList

    muteList = grayList[:thresh] + colorList[thresh:]

    return muteList

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

def plot_metaobjects(scnToPlot, fig, b, t):
    ''' Determines which metaobjects to plot based on scenario '''
    ic('Automatic metaobjects disabled!')
    # if any('ARISE:Control' in scn for scn in scnToPlot): #Triangle for change in number of rlzs
    #     plt.plot(2015, b+(abs(b-t))*0.01, color='#F8A53D', marker='v')
    #     plt.plot(2070, b+(abs(b-t))*0.01, mfc='#F8A53D', mec='#12D0B2', marker='v')
    # if any('GLENS:Control' in scn for scn in scnToPlot):
    #     plt.plot(2030, b+(abs(b-t))*0.01, color='#D93636', marker='v')
    # if any('GLENS:Feedback' in scn for scn in scnToPlot): #Dashed line for model SAI initiation
    #     plt.plot([2020,2020], [b,t], color='#8346C1', linewidth=0.7, linestyle='dashed')
    # if any('ARISE:Feedback' in scn for scn in scnToPlot):
    #     plt.plot([2035,2035], [b,t], color='#12D0B2', linewidth=0.7, linestyle='dashed')

    plt.plot([2020,2020], [b,t], color='#8346C1', linewidth=1.2, linestyle='dashed')
    plt.plot([2035,2035], [b,t], color='#12D0B2', linewidth=1.2, linestyle='dashed')

    return

def find_widest_quantile(darr):
    ''' Figures out whether [0.01Q,-0.01Q] or [-0.99Q,0.99Q] is more appropriate '''
    absQuant99 = np.abs(darr.quantile(0.99))
    absQuant001 = np.abs(darr.quantile(0.01))
    bigVal = np.max([absQuant99,absQuant001])
    widestRange = np.array([-bigVal, bigVal])

    return widestRange
