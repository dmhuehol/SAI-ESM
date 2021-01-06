# import netcdf4 as nc

import xarray as xr

def simple_diff_calc(filename):
    dataPath = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/'
    # filename = dataPath + 'control.001.T.r90x45.shift.annual.nc'
    # endFile = dataPath + 'feedback.020.T.r90x45.shift.annual.nc'

    filename = dataPath + filename
    glensDataset = xr.open_dataset(filename)
    # endDataset = xr.open_dataset(endFile)

    dataValues = glensDataset.T
    # endTemp = endDataset.T

    #Coordinates are time, lon, lat, lev
    #Time increment is one year
    #Available levels are 1000, 500, 300, 50 hPa
    startData = dataValues[0,:,:,:]
    endData = dataValues[len(dataValues)-1,:,:,:]

    differenceVal = endData - startData

    return glensDataset, differenceVal
