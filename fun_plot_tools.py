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
import cmasher
import cmocean as cmocean
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
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
        ic(shrtScn)
        if 'No-SAI' in rlzLoi.attrs['scenario']:
            if 'GLENS' in rlzLoi.attrs['scenario']:
                # ic('GLENS No-SAI')
                toiStartLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["strtIntvl"]["GLENS"][0],
                    setDict["strtIntvl"]["GLENS]"][1])
                toiEndLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["endIntvl"]["GLENS"][0],
                    setDict["endIntvl"]["GLENS"][1])
                toiStart[shrtScn] = toiStartLp
                toiEnd[shrtScn] = toiEndLp
            elif 'CESM2-ARISE' in rlzLoi.attrs['scenario']:
                # ic('ARISE No-SAI')
                toiStartLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["strtIntvl"]["CESM2-ARISE"][0],
                    setDict["strtIntvl"]["CESM2-ARISE"][1])
                toiEndLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["endIntvl"]["CESM2-ARISE"][0],
                    setDict["endIntvl"]["CESM2-ARISE"][1])
                toiStart[shrtScn] = toiStartLp
                toiEnd[shrtScn] = toiEndLp
            elif 'UKESM-ARISE' in rlzLoi.attrs['scenario']:
                # ic('UKESM-ARISE No-SAI')
                toiStartLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["strtIntvl"]["UKESM-ARISE"][0],
                    setDict["strtIntvl"]["UKESM-ARISE"][1])
                toiEndLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["endIntvl"]["UKESM-ARISE"][0],
                    setDict["endIntvl"]["UKESM-ARISE"][1])
                toiStart[shrtScn] = toiStartLp
                toiEnd[shrtScn] = toiEndLp
        elif '/SAI/' in rlzLoi.attrs['scenario']:
            if 'GLENS' in rlzLoi.attrs['scenario']:
                # ic('GLENS SAI')
                toiEndLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["endIntvl"]["GLENS"][0],
                    setDict["endIntvl"]["GLENS"][1])
                toiEnd[shrtScn] = toiEndLp
            elif 'CESM2-ARISE' in rlzLoi.attrs['scenario']:
                # ic('ARISE SAI')
                toiEndLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["endIntvl"]["CESM2-ARISE"][0],
                    setDict["endIntvl"]["CESM2-ARISE"][1])
                toiEnd[shrtScn] = toiEndLp
            elif 'UKESM-ARISE' in rlzLoi.attrs['scenario']:
                # ic('UKESM ARISE SAI')
                toiEndLp = fpd.average_over_years(
                    rlzLoi,
                    setDict["endIntvl"]["UKESM-ARISE"][0],
                    setDict["endIntvl"]["UKESM-ARISE"][1])
                toiEnd[shrtScn] = toiEndLp
        else:
            ic('This should not occur, but does it?')

    return toiStart, toiEnd

def make_scenario_dict(rlzList, setDict):
    ''' Make dictionary accessible by scenario '''
    # Think about whether it makes sense to do this here? Maybe it would be
    # better to make a dictionary object back in bind_scenario
    # This would require substantial rewriting of EVERYTHING, so maybe a long-term
    # thought
    scnDict = dict()
    for darr in rlzList:
        rlzLoi = fpd.obtain_levels(darr, setDict["levOfInt"])
        shrtScn = rlzLoi.scenario.split('/')[len(rlzLoi.scenario.split('/'))-1]
        scnDict[shrtScn] = rlzLoi
        
    return scnDict     
    
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
            image, shrink=.75, orientation="vertical", 
            pad=.02, extend=extent)
            #0-2, 2-10, 10-30, 30-50, 50+
            # ticks=[-50,-30,-10,-2,2,10,30,50])
            # ticks=[-2,-1,0,1,2])
        cb.ax.tick_params(labelsize=6) #def: labelsize=6
        try:
            cb.set_label(data.attrs['units'],size='small')
            # cb.set_label('\u00B0C/yr',size='medium')
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
    
def palette_lists(paletteKey):
    ''' Just some nice colors '''
    zmzm = [ 
        [113/255, 0, 165/255, 1], #286,100,30
        [149/255, 0, 214/255, 1], #285,100,40
        [180/255, 40/255, 255/255, 1], #285,100,50
        [191/255, 104/255, 255/255, 1], #285,100,60
        [204/255, 146/255, 255/255, 1], #285,100,70
        [219/255, 184/255, 255/255, 1], #285,100,80
        [236/255, 220/255, 255/255, 1], #285,100,90
        [255/255, 255/255, 255/255, 1], #285,100,100
        [255/255, 218/255, 216/255, 1], #16,100,90
        [255/255, 179/255, 174/255, 1], #16,100,80
        [255/255, 137/255, 127/255, 1], #16,100,70
        [255/255, 85/255, 59/255, 1], #16,100,60
        [224/255, 54/255, 0, 1], #16,100,50
        [179/255, 41/255, 0, 1], #16,100,40
        [137/255, 29/255, 0, 1]  #16,100,30
    ]
    
    zmzmGray = [ 
        [113/255, 0, 165/255, 1], #286,100,30
        [149/255, 0, 214/255, 1], #285,100,40
        [180/255, 40/255, 255/255, 1], #285,100,50
        [191/255, 104/255, 255/255, 1], #285,100,60
        [204/255, 146/255, 255/255, 1], #285,100,70
        [219/255, 184/255, 255/255, 1], #285,100,80
        [236/255, 220/255, 255/255, 1], #285,100,90
        [227/255, 227/255, 227/255, 1], #285,0,90
        [255/255, 218/255, 216/255, 1], #16,100,90
        [255/255, 179/255, 174/255, 1], #16,100,80
        [255/255, 137/255, 127/255, 1], #16,100,70
        [255/255, 85/255, 59/255, 1], #16,100,60
        [224/255, 54/255, 0, 1], #16,100,50
        [179/255, 41/255, 0, 1], #16,100,40
        [137/255, 29/255, 0, 1]  #16,100,30
    ]
    
    colors = {
        "zmzm": zmzm,
        "zmzmGray": zmzmGray
    }
    palette = colors[paletteKey]
    
    return palette

def line_from_scenario(scn, md):
    ''' Get line color and label from scenario information '''
    if 'GLENS:Control' in scn:
        activeColor = '#D93636' #Red
        activeLabel = md['cntrlStr']
    elif 'GLENS:Feedback' in scn:
        activeColor = '#8346C1' #Purple
        activeLabel = md['fdbckStr']
    elif 'CESM2-ARISE:Feedback' in scn:
        activeColor = '#12D0B2' #Turquoise
        activeLabel = md['ariseStr']
    elif 'ARISE-DelayedStart:Feedback' in scn:
        activeColor = '#12D0B2'
        activeLabel = md["arisedsStr"]
    elif 'ARISE-1.0:Feedback' in scn:
        activeColor = '#12D0B2'
        activeLabel = md["arise1p0Str"]
    elif 'CESM2-ARISE:Control' in scn:
        activeColor = '#F8A53D' #Orange
        activeLabel = md['s245Cntrl']
    elif 'UKESM-ARISE:Feedback' in scn:
        activeColor = '#12D0B2'
        activeLabel = md['ukAriseStr']
    elif 'UKESM-ARISE:Control' in scn:
        activeColor = '#F8A53D'
        activeLabel = md['ukS245Str']
    else:
        activeColor = '#000000'
        activeLabel = 'Unknown'
        ic('Unknown scenario! Plotting with black line and unknown label.')

    return activeColor, activeLabel

def plot_metaobjects(scnToPlot, fig, b, t, lw=1.2):
    ''' Determines which metaobjects to plot based on scenario '''
    # ic('Automatic metaobjects disabled!')
    # Commented block plots little triangles denoting change in ensemble size
    # if any('ARISE:Control' in scn for scn in scnToPlot):
    #     plt.plot(2015, b+(abs(b-t))*0.01, color='#F8A53D', marker='v')
    #     plt.plot(2070, b+(abs(b-t))*0.01, mfc='#F8A53D', mec='#12D0B2', marker='v')
    # if any('GLENS:Control' in scn for scn in scnToPlot):
    #     plt.plot(2030, b+(abs(b-t))*0.01, color='#D93636', marker='v')

    # Commented block only plots vertical line denoting deployment for scenarios present
    # if any('GLENS:Feedback' in scn for scn in scnToPlot): #Dashed line for model SAI initiation
    #     plt.plot([2020,2020], [b,t], color='#8346C1', linewidth=0.7, linestyle='dashed')
    # if any('ARISE:Feedback' in scn for scn in scnToPlot):
    #     plt.plot([2035,2035], [b,t], color='#12D0B2', linewidth=0.7, linestyle='dashed')

    # Always plot vertical lines denoting deployment in 2020, 2035
    plt.plot([2020,2020], [b,t], color='#8346C1', linewidth=lw, linestyle='dashed')
    plt.plot([2035,2035], [b,t], color='#12D0B2', linewidth=lw, linestyle='dashed')
    if ('ARISE-DelayedStart:Feedback' in scn for scn in scnToPlot):
        plt.plot([2045,2045], [b,t], color='#12D0B2', linewidth=lw, linestyle='dashed')

    return

def find_widest_quantile(darr):
    ''' Figures out whether [0.01Q,-0.01Q] or [-0.99Q,0.99Q] is more appropriate '''
    absQuant99 = np.abs(darr.quantile(0.99))
    absQuant001 = np.abs(darr.quantile(0.01))
    bigVal = np.max([absQuant99,absQuant001])
    widestRange = np.array([-bigVal, bigVal])

    return widestRange
    
def get_trifurcate_colormap(c1, c2, c3):
    ''' Create trifurcated colormap for decadal climate distance '''
    # Range parameters below are hard-coded for 0-1, 1-10, 10+
    pt1 = cmasher.get_sub_cmap(c1, 0.2, 0.2)
    pt2 = cmasher.get_sub_cmap(c2, 0.5, 1)
    pt3 = c3
    
    map1 = pt1(np.linspace(0, 1, 100)) #0 to 1
    map2 = pt2(np.linspace(0, 1, 900)) #1 to 10
    map3 = pt3
    
    colors = np.vstack((map1, map2, map3))
    tri_cmap = mcolors.LinearSegmentedColormap.from_list('tri_map', colors)

    return tri_cmap
    
def get_trifurcate_div_colormap(c1, c2, c3):
    ''' Create trifurcated colormap for decadal climate distance '''
    # Range parameters below are hard-coded for 0-1, 1-10, 10+
    pt1 = c1
    pt2 = cmasher.get_sub_cmap(c2, 0, 1)
    pt3 = c3
    
    map1 = pt1 #0 to 1
    map2 = pt2(np.linspace(0, 1, 100)) #1 to 10
    map3 = pt3
    
    colors = np.vstack((map1, map2, map3))
    tri_cmap = mcolors.ListedColormap(colors)

    return tri_cmap
    
def get_cspd_colormap(palKey):
    ''' Create custom colormap for climate velocity '''
    #What we really want is 0-2, 2-10, 10-30, 30-50, 50+
    zmzm = palette_lists(palKey)
    ic(len(zmzm))
    map1 = np.reshape(
        np.tile(zmzm[0], 1), (1,4))
    map2 = np.reshape(
        np.tile(zmzm[1], 20), (20,4))
    map3 = np.reshape(
        np.tile(zmzm[3], 20), (20,4))
    map4 = np.reshape(
        np.tile(zmzm[5], 8), (8,4))
    map5 = np.reshape(
        np.tile(zmzm[7], 4), (4,4))
    map6 = np.reshape(
        np.tile(zmzm[9], 8), (8,4))
    map7 = np.reshape(
        np.tile(zmzm[11], 20), (20,4))
    map8 = np.reshape(
        np.tile(zmzm[13], 20), (20,4))
    map9 = np.reshape(
        np.tile(zmzm[14], 1), (1,4))
    
    colors = np.vstack(
        (map1, map2, map3, map4, map5, map6, map7, map8, map9))
    custom_cmap = mcolors.ListedColormap(colors)

    return custom_cmap