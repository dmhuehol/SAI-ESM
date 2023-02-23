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
    ''' Calculate temporal gradient for an array. Vectorized linear
    regression is strongly based on implementation found at:
    hrishichandanpurkar.blogspot.com/2017/09/vectorized-functions-for-correlation.html
    '''
    startTime = cftime.DatetimeNoLeap(years[0], 7, 15, 12, 0, 0, 0)
    endTime = cftime.DatetimeNoLeap(years[1], 7, 15, 12, 0, 0, 0)
    
    xData = darr.time.dt.year
    timeEntries = len(darr.time)
    xTimeMn = xData.mean(dim='time')
    yTimeMn = darr.mean(dim='time')
    xTimeStdev  = xData.std(dim='time')
    yTimeStdev  = darr.std(dim='time')
    cov = np.sum(
        (xData - xTimeMn) * (darr - yTimeMn), axis=0) / (timeEntries)
    regSlp = cov/(xTimeStdev**2)
    regInt = yTimeMn - xTimeMn*regSlp
    # ic(cov, regSlp, regInt) # Troubleshooting

    temporalGrad = {
        "grad": regSlp,
        "intcpt": regInt
    }

    return temporalGrad

def calc_spatial_grad(darr):
    ''' Calculate spatial gradient '''
    ic()
    lats = darr.lat.data
    lons = darr.lon.data
    earthRad = 6371000 / 1000 #Earth's radius in km

    lonGrid, latGrid = np.meshgrid(lons, lats)
    latGridRad = np.deg2rad(latGrid)
    dLat = np.gradient(latGrid)[0]
    dLon = np.gradient(lonGrid)[1]
    boxLngth = earthRad * np.pi / 180 # Box length in km
    distY = dLat * boxLngth # Lat distance in km
    distX = dLon * boxLngth * np.cos(latGridRad) # Lon distance in km
    
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
        nsGrad = np.nanmean(nsGradRlzList, axis=0)
        ewGrad = np.nanmean(ewGradRlzList, axis=0)
    else:
        nsGrad = sndimg.sobel(darr, axis=0, mode='reflect')
        ewGrad = sndimg.sobel(darr, axis=1, mode='reflect')
    nsGradSc = nsGrad / (8 * distY)
    ewGradSc = ewGrad / (8 * distX)
    # Total spatial gradient
    totGrad = np.sqrt((nsGradSc ** 2) + (ewGradSc ** 2))
    spatGrad = totGrad

    return spatGrad

def calc_climate_speed(darr, setDict):
    ''' Calculate the climate speed of a variable '''
    if 'CESM2' in darr.scenario:
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
    climSpd.attrs['long_name'] = 'Climate speed of 2m temperature' \
        + ' ' + str(setYrs[0]) + '-' + str(setYrs[1])
    climSpd.attrs['units'] = 'km/yr'
    
    # Display for troubleshooting
    ic(check_stats(tGrad.data))
    # ic(check_stats(sGrad.data))
    # ic(check_stats(climSpd.data))
    # ic(climSpd.attrs)
    return climSpd

def calc_decadal_climate_distance(darr, setDict):
    ''' Calculate the decadal climate distance of a variable '''
    climSpd = calc_climate_speed(darr, setDict)

    dcd = np.abs(climSpd) * 10
    dcd.attrs = darr.attrs
    # TODO: make attribute generation flexible by variable
    dcd.attrs['long_name'] = climSpd.attrs['long_name'].replace(
        'Climate speed', 'Decadal climate distance')
    dcd.attrs['units'] = 'km'

    return dcd

### Helper functions
def check_stats(darr):
    # darr = darr.compute()

    statDict = {
        "shape": np.shape(darr),
        "med": np.nanmedian(darr),
        "mean": np.nanmean(darr),
        "max": np.nanmax(darr),
        "min": np.nanmin(darr)
    }

    return statDict
