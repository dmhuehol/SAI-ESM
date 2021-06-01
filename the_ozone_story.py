"""theOzoneStory
This makes figures showing change in stratospheric ozone in maxed out aesthetics
format.

Written by Daniel Hueholt | June 2021
Graduate Research Assistant at Colorado State University
"""

from icecream import ic
import sys

import xarray as xr
xr.set_options(keep_attrs=True)
import matplotlib
import matplotlib.font_manager as fm
matplotlib.font_manager._rebuild()
fontPath = '/Users/dhueholt/Library/Fonts/'  # the location of the font file
firaSansPath = {"Thin": 'FiraSans-Thin.ttf',
                "Medium": 'FiraSans-Medium.ttf'}
robotoPath = {"Regular": 'Roboto-Regular.ttf'}
Roboto = fm.FontProperties(fname=fontPath+robotoPath['Regular'])  # get the font based on the font_path
FiraSansMed = fm.FontProperties(fname=fontPath+firaSansPath['Medium'])
FiraSansThin = fm.FontProperties(fname=fontPath+firaSansPath['Thin'])
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy
import cartopy.crs as ccrs
import numpy as np

import difference_over_time as dot
import process_glens_fun as pgf
import plotting_tools as plt_tls
import fun_convert_unit as fcu

# Inputs
dataPath = '/Users/dhueholt/Documents/GLENS_data/annual_o3/'
filenameCntrl = 'control_003_O3_201001-201912_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
filenameFdbck = 'feedback_003_O3_202001-202912_203001-203912_204001-204912_205001-205912_206001-206912_207001-207912_208001-208912_209001-209912_annual.nc'
cntrlPath = dataPath + filenameCntrl
fdbckPath = dataPath + filenameFdbck

startInt = [2010,2019]
finalInt = [2090,2099]

# Open data
glensDsetCntrl = xr.open_dataset(cntrlPath)
glensDsetFdbck = xr.open_dataset(fdbckPath)
dataKey = pgf.discover_data_var(glensDsetCntrl)
glensDarrCntrl = glensDsetCntrl[dataKey]
glensDarrFdbck = glensDsetFdbck[dataKey]

savePath='/Users/dhueholt/Documents/GLENS_fig/20210601_theOzoneStory/'
dpi_val = 800

def to3s_dg(levOfInt,tStr):
    ''' Makes difference globe '''

    savePrfx='NORM_globe_1p_FdbckCntrl_'

    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensDarrCntrl, levOfInt)
    glensFdbckLoi = pgf.obtain_levels(glensDarrFdbck, levOfInt)

    # Average over years
    toiStart = dot.average_over_years(glensCntrlLoi, startInt[0], startInt[1]) # 2010-2019 is baseline, injection begins 2020
    toiEndFdbck = dot.average_over_years(glensFdbckLoi, finalInt[0], finalInt[1])
    diffToiFdbck =  toiEndFdbck - toiStart

    # Unit conversion
    diffToiFdbckPlot = fcu.molmol_to_ppm(diffToiFdbck)
    toiStartUnit = fcu.molmol_to_ppm(toiStart)
    diffToiFdbckPlotNorm = diffToiFdbckPlot / np.max(toiStartUnit) * 100

    # Plotting
    lastDcd = str(finalInt[0]) + '-' + str(finalInt[1])
    levStr = pgf.make_level_string(glensCntrlLoi, levOfInt)

    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = 'PuOr'
    minVal = -15
    maxVal = 15

    plt_tls.drawOnGlobe(ax, diffToiFdbckPlotNorm, glensDarrFdbck.lat, glensDarrFdbck.lon, cmap, vmin=minVal, vmax=maxVal, cbarBool=True, fastBool=True, extent='max')
    plt.title(tStr,fontproperties=FiraSansMed,fontsize=18) #Override automatic title generation here
    saveStr = savePrfx + dataKey + '_' + levStr + '_' + lastDcd
    savename = savePath + saveStr + '.png'
    plt.savefig(savename, dpi=dpi_val, bbox_inches='tight')
    ic(savename)

to3s_dg(levOfInt='stratosphere',tStr="Percent change in stratospheric ozone 2010-2019 to 2090-2099")
to3s_dg(levOfInt=[250,50],tStr="Percent change in 250mb to 50mb ozone 2010-2019 to 2090-2099")

def to3s_ts():
    ''' Make stratospheric ozone timeseries '''
    levOfInt = 'stratosphere'
    regionToPlot = 'global'
    saveFile = 'timeseries_O3_'
    saveName = savePath + saveFile

    bndDct = pgf.find_matching_year_bounds(glensDarrCntrl, glensDarrFdbck)
    glensCntrlPoi = glensDarrCntrl[bndDct['cntrlStrtMtch']:bndDct['cntrlEndMtch']+1] #RANGES IN PYTHON ARE [)
    glensFdbckPoi = glensDarrFdbck[bndDct['fdbckStrtMtch']:bndDct['fdbckEndMtch']+1]

    # Obtain levels
    glensCntrlPoi = pgf.obtain_levels(glensCntrlPoi, levOfInt)
    glensFdbckPoi = pgf.obtain_levels(glensFdbckPoi, levOfInt)

    # Deal with area
    cntrlToPlot, locStr, locTitleStr = pgf.manage_area(glensCntrlPoi, regionToPlot, areaAvgBool=True)
    fdbckToPlot, locStr, locTitleStr = pgf.manage_area(glensFdbckPoi, regionToPlot, areaAvgBool=True)

    # Unit conversion
    cntrlToPlot = fcu.molmol_to_ppm(cntrlToPlot)
    fdbckToPlot = fcu.molmol_to_ppm(fdbckToPlot)

    # Plotting
    yStr = 'parts per million'
    levStr = pgf.make_level_string(glensCntrlPoi, levOfInt)

    # Make timeseries
    fig, ax = plt.subplots()
    plt.plot(bndDct['mtchYrs'],cntrlToPlot.data,color='#DF8C20',label='RCP8.5') #These are the cuckooColormap colors
    plt.plot(bndDct['mtchYrs'],fdbckToPlot.data,color='#20DFCC',label='SAI')
    ax.annotate('Global stratospheric ozone concentration',(0.02,0.935),(0.02,0.935),'axes fraction',fontproperties=FiraSansMed,fontsize=11,color='#000000')
    ax.annotate('RCP8.5',(0.236,0.482),(0.236,0.482),'axes fraction',fontproperties=FiraSansMed,fontsize=11,color='#DF8C20')
    ax.annotate('SAI',(0.419,0.28),(0.419,0.28),'axes fraction',fontproperties=FiraSansMed,fontsize=11,color='#20DFCC')
    plt.ylabel(yStr,fontproperties=FiraSansThin)
    for label in ax.get_xticklabels():
        label.set_fontproperties(FiraSansThin)
    for label in ax.get_yticklabels():
        label.set_fontproperties(FiraSansThin)
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.savefig(saveName + locStr + '_' + levStr + '.png',dpi=dpi_val,bbox_inches='tight')
    ic(saveName + locStr + '_' + levStr + '.png')

    print('Completed!')

to3s_ts()
