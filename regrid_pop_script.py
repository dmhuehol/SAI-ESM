''' regrid_pop_script
Regrids Parallel Ocean Model (POP) output from the B-grid to a standard lat/lon
using interpolation.

Originally written by Emily Gordon based on code provided by Zachary Labe
Modified by Daniel Hueholt
'''

from icecream import ic
import sys

import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import scipy.interpolate as interp

def extract_pop_latlons(popFile):
    ''' Extract lat/lons from the POP grid '''
    popgrid = xr.open_dataset(popFile)
    popLat = np.asarray(popgrid.TLAT.data)
    popLon = np.asarray(popgrid.TLONG.data)

    return popLat, popLon

def regrid(dataIn,latIn,lonIn,latOut,lonOut):
    ''' Takes POP model output on a B-grid (T-cells for scalars,
    U-cells for vectors) (latInxlonIn) and regrids to (latOutxlonOut) '''

    lonOut,latOut = np.meshgrid(lonOut,latOut) # make grid
    dataRg = np.ravel(dataIn) # move inputs to vectors
    lonRg = np.ravel(lonIn)
    latRg = np.ravel(latIn)
    dataInterpRg = interp.griddata((latRg,lonRg),dataRg,(latOut,lonOut), method='linear')

    return dataInterpRg

# Make list of data to index
dataPath = '/glade/scratch/dhueholt/sept_IFRAC/'
strList = sorted(glob.glob(dataPath + "*.nc"))
dataVar = 'IFRAC'

inPop = '/glade/work/dhueholt/grids/control_IFRAC_useForGrid.nc'
popLat,popLon = extract_pop_latlons(inPop)

for strc, strv in enumerate(strList):
    popDataArray = xr.open_dataset(strv,decode_times=False)
    varOfInt = popDataArray[dataVar]
    time = popDataArray.time
    latNew = np.arange(-90,91)
    lonNew = np.arange(0,360)

    dataRegrid = np.empty((varOfInt.shape[0],latNew.shape[0],lonNew.shape[0]))
    for bc in range(varOfInt.shape[0]):
        if bc % 10 == 0: # progress indicator
            prog = np.round(bc / varOfInt.shape[0],3)
            ic(prog)
        newDat = regrid(varOfInt[bc,:,:],popLat,popLon,latNew,lonNew)
        dataRegrid[bc,:,:] = newDat

    newDset = xr.Dataset(
        {dataVar: (("time","lat","lon"), dataRegrid)},
        coords={
            "time": time,
            "lat": latNew,
            "lon": lonNew
        }
    )
    newDset.attrs = popDataArray.attrs
    newDset[dataVar].attrs = popDataArray[dataVar].attrs

    originalFilename = strv[0:len(strv)-3]
    strOut =  originalFilename + '_RG' + '.nc' #originalfilename_RG.nc
    newDset.to_netcdf(strOut)
    ic(strOut)
