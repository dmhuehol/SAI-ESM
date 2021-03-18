import cartopy as ct
import cartopy.crs as ccrs
import cmocean as cmocean
import numpy as np
import numpy.ma as ma
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import seaborn as sn


# drawOnGlobe and add_cyclic_point written by Prof. Elizabeth Barnes at Colorado State University
# add_cyclic_point copied from cartopy utils

def drawOnGlobe(ax, data, lats, lons, cmap='coolwarm', vmin=None, vmax=None, inc=None, cbarBool=True, contourMap=[], contourVals = [], fastBool=False, extent='both'):

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
    else:
        cb = None

    image.set_clim(vmin,vmax)

    return cb, image

def add_cyclic_point(data, coord=None, axis=-1):

    # had issues with cartopy finding utils so copied for myself

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

def plot_pdf_kdeplot(handles, colors, labels, savePath, saveName, dpiVal=400):
# Plot kde pdfs for several input handles
    plt.figure()
    print('Plotting!')
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            # print(ind)
            # print(h)
            # print(labels[ind])
            # print(colors[ind])
            ax = sn.kdeplot(data=h, label=labels[ind], color=colors[ind], linewidth=2)
            # ax = sn.kdeplot(data=h.data.flatten(), label=labels[ind], color=colors[ind], linewidth=2)
            ax.set(xlabel='Celsius', ylabel='Density')
    else:
        ax = sn.kdeplot(data=handles, label=labels, color=colors, linewidth=2)
        ax.set(xlabel='Celsius', ylabel='Density')

    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])
    # plt.show()
    plt.savefig(savePath + saveName + '.png', dpi=dpiVal, bbox_inches='tight')

def plot_pdf_hist(handles, colors, labels, savePath, saveName, binwidth, dpiVal=400):
# Plot histogram pdfs for several input handles
    plt.figure()
    print('Plotting!')
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            # print(ind)
            # print(h)
            # print(labels[ind])
            # print(colors[ind])
            ax = sn.histplot(data=h, label=labels[ind], color=colors[ind], edgecolor='#3B3B3B', linewidth=0.8, kde=False, binwidth=binwidth)
            # ax = sn.histplot(data=h.data.flatten(), label=labels[ind], color=colors[ind], edgecolor='#3B3B3B', linewidth=0.8, kde=False, binwidth=binwidth)
            ax.set(xlabel='Celsius', ylabel='Density')
    else:
        ax = sn.histplot(data=handles, label=labels, color=colors, edgecolor='#3B3B3B', stat='density', linewidth=0.8, kde=False, binwidth=binwidth)
        ax.set(xlabel='Celsius', ylabel='Density')

    plt.legend(bbox_to_anchor=(0.83,-0.1), ncol=2, fontsize=8)
    plt.title(labels[np.size(labels)-1])
    # plt.show()
    plt.savefig(savePath + saveName + '.png', dpi=dpiVal, bbox_inches='tight')

def select_colors(baselineFlag, nFdbck, nCntrl):
# Returns colors for baseline and a set number of feedback and control objects
# to plot. Aimed at the pdf plots but may be more broadly useful.
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
# Fill colors to plot by number
    if nfc == 1:
        colorsToPlot.append(colList[3])
    elif nfc == 2:
        colorsToPlot.append(colList[2])
        colorsToPlot.append(colList[5])
    elif nfc == 3:
        colorsToPlot.append(colList[1])
        colorsToPlot.append(colList[3])
        colorsToPlot.append(colList[6])
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
# Generate labels for figure titles and output filenames
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
        # handlesToPlot[cdv]["data"] = darr[startInd:endInd]
        # handlesToPlot.append(darr[startInd:endInd])

# labelsToPlot = ['2010-2019 Baseline', '2050-2059 RCP8.5', '2090-2099 RCP8.5', '2050-2059 SAI', '2090-2099 SAI', 'Global SST PDFs in GLENS']
#
