''' fun_calc_var
Contains functions for calculating variables which are more complicated than
unit conversions, but not complicated enough to require saving the data
separately.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import cftime
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as sndimg
from sklearn.linear_model import LinearRegression
import xarray as xr

def calc_temporal_grad(darr, years):
    ''' Calculate temporal gradient by estimating 'heights' between a start and
    end time. Vectorized linear regression is strongly based on
    hrishichandanpurkar.blogspot.com/2017/0 9/vectorized-functions-for-correlation.html
    implementation '''
    startTime = cftime.DatetimeNoLeap(years[0], 7, 15, 12, 0, 0, 0)
    endTime = cftime.DatetimeNoLeap(years[1], 7, 15, 12, 0, 0, 0)

    xData = darr.time.dt.year

    timeEntries = len(darr.time)
    xTimeMn = xData.mean(dim='time')
    yTimeMn = darr.mean(dim='time')
    xTimeStdev  = xData.std(dim='time')
    yTimeStdev  = darr.std(dim='time')
    # ic(timeEntries, xTimeMn, yTimeMn, xTimeStdev, yTimeStdev)
    cov = np.sum((xData - xTimeMn)*(darr - yTimeMn), axis=0)/(timeEntries)
    regSlp = cov/(xTimeStdev**2)
    regInt = yTimeMn - xTimeMn*regSlp
    # ic(cov, regSlp, regInt)

    temporalGrad = {
        "grad": regSlp,
        "intcpt": regInt
    }

    return temporalGrad

def calc_spatial_grad(darr):
    ''' Calculate spatial gradient '''
    latNew = darr.lat.data
    lonNew = darr.lon.data
    earthRad = 6371000 / 1000 #Earth's radius in km

    lonGrid,latGrid = np.meshgrid(lonNew, latNew)
    latGridRad = np.deg2rad(latGrid)
    dLat = np.gradient(latGrid)
    dLon = np.gradient(lonGrid)
    distY = dLat[0] * earthRad * np.pi / 180
    distX = (dLon[1]/180) * np.pi * earthRad * np.cos(latGridRad)

    # Calculate north/south spatial gradient
    nsGrad = sndimg.sobel(darr, axis=0, mode='reflect')
    nsGradSc = nsGrad / (8 * distY)
    # Calculate east/west spatial gradient
    ewGrad = sndimg.sobel(darr, axis=1, mode='reflect')
    ewGradSc = ewGrad / (8 * distX)
    # Total spatial gradient
    totGrad = np.sqrt((nsGradSc ** 2) + (ewGradSc ** 2))
    spatGrad = np.arctan(totGrad) * 57.29578

    return spatGrad

def calc_climate_speed(darr, setDict):
    ''' Calculate the climate speed of a variable '''
    if 'CESM2-ARISE' in darr.scenario:
        setYrs = setDict["calcIntvl"]["CESM2-ARISE"]
    elif 'UKESM-ARISE' in darr.scenario:
        setYrs = setDict["calcIntvl"]["UKESM-ARISE"]
    elif 'GLENS' in darr.scenario:
        setYrs = setDict["calcIntvl"]["GLENS"]
    timeSliceCesm = slice(
        cftime.DatetimeNoLeap(setYrs[0], 7, 15, 12, 0, 0, 0),
        cftime.DatetimeNoLeap(setYrs[1], 7, 15, 12, 0, 0, 0))
    timeSliceUkesm = slice(
        cftime.Datetime360Day(setYrs[0], 6, 30, 0, 0, 0, 0),
        cftime.Datetime360Day(setYrs[1], 6, 30, 0, 0, 0, 0))
    try:
        darrToi = darr.sel(time=timeSliceCesm)
    except: # This will get messy when another model runs ARISE :)
        darrToi = darr.sel(time=timeSliceUkesm)

    tGradDict = calc_temporal_grad(darrToi, years=setYrs)
    tGrad = tGradDict["grad"].compute()
    darrTimeMn = darrToi.mean(dim='time')
    sGrad = calc_spatial_grad(darrTimeMn)

    climSpd = tGrad / sGrad
    climSpd.attrs = darr.attrs
    # TODO: make attribute generation flexible by variable
    climSpd.attrs['long_name'] = 'Climate speed of 2m temperature'
    climSpd.attrs['units'] = 'deg C / km'

    return climSpd

def calc_decadal_climate_distance(darr, setDict):
    ''' Calculate the decadal climate distance of a variable '''
    climSpd = calc_climate_speed(darr, setDict)

    dcd = np.abs(climSpd) * 10
    dcd.attrs = darr.attrs
    # TODO: make attribute generation flexible by variable
    dcd.attrs['long_name'] = 'Decadal climate distance of 2m temperature'
    dcd.attrs['units'] = 'km'

    return dcd

### Helper functions
def check_stats(darr):
    # darr = darr.compute()

    statDict = {
        "med": np.median(darr),
        "mean": np.mean(darr),
        "max": np.max(darr),
        "min": np.min(darr)
    }

    return statDict
