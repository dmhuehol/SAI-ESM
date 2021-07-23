''' difference_over_time
Contains functions to calculate values related to changes over time. For example,
can average over a given number of years.

Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
'''
from icecream import ic

import xarray as xr
import numpy as np

def average_over_years(glensDarr,startYear,endYear):
    ''' Take average over a period of time '''
    datasetYears = glensDarr['time'].dt.year.data
    startInd = int(np.where(datasetYears == startYear)[0])
    endInd = int(np.where(datasetYears == endYear)[0])
    intervalOfInterest = glensDarr[startInd:endInd]
    glensMeanToi = intervalOfInterest.mean(dim='time')

    return glensMeanToi


def simple_diff_calc(glensDataset):
    ''' Calculate difference from start to end timestep '''
    dataValues = glensDataset.T
    startData = dataValues[0,:,:,:]
    endData = dataValues[len(dataValues)-1,:,:,:]
    differenceVal = endData - startData

    return glensDataset, differenceVal
