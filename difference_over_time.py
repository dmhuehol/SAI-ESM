import xarray as xr
import numpy as np

def average_over_years(glensDarr,startYear,endYear):
    # dataValues = glensDataset.T
    datasetYears = glensDarr['time'].dt.year.data
    # print(datasetYears)
    startInd = int(np.where(datasetYears == startYear)[0])
    endInd = int(np.where(datasetYears == endYear)[0])
    intervalOfInterest = glensDarr[startInd:endInd] #HERE IT IS restore with # intervalOfInterest = glensDataset.T[startInd:endInd]
    glensMeanToi = intervalOfInterest.mean(dim='time')

    return glensMeanToi

def average_over_years_temp(glensDataset,startYear,endYear):
    # dataValues = glensDataset.T
    datasetYears = glensDataset['time'].dt.year.data
    # print(datasetYears)
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
