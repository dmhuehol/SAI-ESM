# Try spatial gradient with CRU_TS4
from icecream import ic
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as sndimg
import xarray as xr

import cartopy as ct
import cmasher
import seaborn

import fun_calc_var as fcv

import matplotlib.colors

def drawOnGlobe(
        ax, data, lats, lons, cmap='coolwarm', vmin=None, vmax=None, inc=None,
        cbarBool=True, contourMap=[], contourVals = [], fastBool=True,
        extent='both', addCyclicPoint=False, alph=1):
    ''' Draws geolocated data on a globe. Originally written by Prof. Elizabeth Barnes at
    Colorado State University, edited by Daniel Hueholt '''
    data_crs = ct.crs.PlateCarree()
    if addCyclicPoint: # Add cyclic point to prime meridian for ocean data
        data_cyc, lons_cyc = add_cyclic_point(data, coord=lons, axis=-1)
    else:
        data_cyc = data
        lons_cyc = lons
    ax.set_global()
    ax.coastlines(linewidth = 1.2, color='black')
    import matplotlib.colors as colors
    if(fastBool):
        image = ax.pcolormesh(
            lons_cyc, lats, data_cyc, transform=data_crs, cmap=cmap, alpha=alph, norm=colors.LogNorm(vmin=0.0001, vmax=0.1))
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
            cb.set_label('degC/km', size='medium')
            # cb.set_label('\u00B0C',size='medium')
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

outDict = {
    "outPath": "/Users/dhueholt/Documents/ecology_fig/20230222_lineByLine/",
    "dpi": 400
}

def calc_spatial_grad(darr):
    ''' Calculate spatial gradient '''
    latNew = darr.lat.data
    lonNew = darr.lon.data
    earthRad = 6371000 / 1000 #Earth's radius in km

    lonGrid,latGrid = np.meshgrid(lonNew, latNew)
    latGridRad = np.deg2rad(latGrid)
    dLat = np.gradient(latGrid)
    ic(dLat[0])
    dLon = np.gradient(lonGrid)
    ic(dLon)
    distY = dLat[0] * earthRad * (np.pi / 180)
    distX = dLon[1] * (np.pi/180) * earthRad * np.cos(latGridRad)

    if 'realization' in darr.dims:
        nsGradRlzList = list()
        ewGradRlzList = list()
        for r in darr.realization.data:
            nsGradRlz = sndimg.sobel( # Calculate north/south spatial gradient
                darr.isel(realization=r), axis=0, mode='reflect')
            nsGradRlzList.append(nsGradRlz)
            ewGradRlz = sndimg.sobel( # Calculate east/west spatial gradient
                darr.isel(realization=r), axis=1, mode='reflect')
            ewGradRlzList.append(ewGradRlz)
        nsGrad = np.mean(nsGradRlzList, axis=0)
        ewGrad = np.mean(ewGradRlzList, axis=0)
    else:
        nsGrad = sndimg.sobel(darr, axis=0, mode='reflect')
        ewGrad = sndimg.sobel(darr, axis=1, mode='reflect')
    nsGradSc = nsGrad / (8 * distY)
    ewGradSc = ewGrad / (8 * distX)
    # Total spatial gradient
    totGrad = np.sqrt((nsGradSc ** 2) + (ewGradSc ** 2))#np.sqrt((nsGradSc ** 2) + (ewGradSc ** 2))
    spatGrad = totGrad #* 57.29578

    return spatGrad

inPath = '/Users/dhueholt/Documents/ecology_data/cru_ts4/'
fname = 'cru_ts4.06.1961.1970.tmp.dat.nc'

cruts = xr.open_dataset(inPath+fname)
crutmp = cruts['tmp']
# ic(crutmp)

crutmpToi = crutmp.mean(dim='time')
ic(crutmpToi)

# a1 = np.array([50, 45, 50])
# a2 = np.array([30, 30, 30])
# a3 = np.array([8, 10, 10])
# crutmpToi = np.vstack((a1, a2, a3))
# ic(crutmpToi)
# ns = sndimg.sobel(crutmpToi, axis=0)
# ew = sndimg.sobel(crutmpToi, axis=1)
# ic(ns, ew)
# sys.exit('STOP')

spatGrad = calc_spatial_grad(crutmpToi)
ic(spatGrad)
ic(fcv.check_stats(spatGrad))
# spatGrad10 = spatGrad/10
# ic(fcv.check_stats(spatGrad10))

# Plotting
CL = 0.
proj = ct.crs.EqualEarth(central_longitude = CL)
lat = crutmp.lat
lon = crutmp.lon

import fun_plot_tools as fpt
# tri_map = fpt.get_wacky_colormap('cmr.flamingo', 'Oranges', 'Reds')
sCm = 'turbo'
sFig = plt.figure()
plt.rcParams.update({'font.family': 'Open Sans'})
sAx = plt.subplot(1, 1, 1, projection=proj)
sCl = [0.0001, 0.1]
sCb, sIm = drawOnGlobe(
    sAx, spatGrad, lat, lon, cmap=sCm,
    vmin=sCl[0], vmax=sCl[1])
plt.title('Spatial grad of temperature 1961-1970 CRU_TS4')
sSavename = outDict["outPath"] + "1961-1970" + "_sGradno57_cru_globe.png"
plt.savefig(sSavename, dpi=outDict["dpi"], bbox_inches='tight')


ic('Hooray!')