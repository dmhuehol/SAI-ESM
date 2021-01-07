import xarray as xr
import numpy as np

def import_glens_dataset(filename):
    dataPath = '/Users/dhueholt/Documents/ANN_GeoEng/data/GLENS/'

    filename = dataPath + filename
    glensDataset = xr.open_dataset(filename)

    return glensDataset

def average_over_years(glensDataset,startYear,endYear):
    # dataValues = glensDataset.T
    datasetYears = glensDataset['time'].dt.year.data
    startInd = int(np.where(datasetYears == startYear)[0])
    endInd = int(np.where(datasetYears == endYear)[0])
    intervalOfInterest = glensDataset.T[startInd:endInd]
    glensMeanToi = intervalOfInterest.mean(dim='time')

    return glensMeanToi


def simple_diff_calc(glensDataset):
    dataValues = glensDataset.T
    # endTemp = endDataset.T

    #Coordinates are time, lon, lat, lev
    #Time increment is one year
    #Available levels are 1000, 500, 300, 50 hPa
    startData = dataValues[0,:,:,:]
    endData = dataValues[len(dataValues)-1,:,:,:]

    differenceVal = endData - startData

    return glensDataset, differenceVal
