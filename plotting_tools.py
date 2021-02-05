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

def plot_pdf_kdeplot(handles, colors, bins, labels, savePath, saveName):
    plt.figure()
    print('Plotting!')
    if np.size(colors) > 1:
        for ind, h in enumerate(handles):
            hHist = np.histogram(h, bins[ind])
            xVals = hHist[1][:-1]
            pdfVals = hHist[0].astype(float)/np.size(h)
            ax = sn.kdeplot(data=h, label=labels[ind], color=colors[ind], linewidth=2)
            ax.set(xlabel='Standard deviations', ylabel='Density')
    else:
        hHist = np.histogram(handles,bins)
        xVals = hHist[1][:-1]
        pdfVals = hHist[0].astype(float)/np.size(handles)
        ax = sn.kdeplot(data=handles, label=labels, color=colors, linewidth=2)
        ax.set(xlabel='Standard deviations',ylabel='Density')

    plt.legend()
    plt.title(labels[np.size(labels)-1])
    # plt.show()
    plt.savefig(savePath + saveName + '.png', dpi=300)
