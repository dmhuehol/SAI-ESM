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

import fun_process_data as fpd
import fun_derive_data as fdd

def time_slice_darr(darr, setDict):
    '''Make darr over times of interest from input years'''
    # The order of these if statements is important!
    if 'CESM2-ARISE-DelayedStart' in darr.scenario:
        setYrs = setDict["calcIntvl"]["CESM2-ARISE-DelayedStart"]
    elif 'CESM2-ARISE:Control' in darr.scenario:
        setYrs = setDict["calcIntvl"]["CESM2-SSP245"]
    elif 'CESM2-ARISE' in darr.scenario:
        setYrs = setDict["calcIntvl"]["CESM2-ARISE"]
    elif 'UKESM-ARISE' in darr.scenario:
        setYrs = setDict["calcIntvl"]["UKESM-ARISE"]
    elif 'GLENS' in darr.scenario:
        setYrs = setDict["calcIntvl"]["GLENS"]
    elif 'PreindustrialControl' in darr.scenario:
        setYrs = setDict["calcIntvl"]["piControl"]

    darrToiItvlList = list()
    for sy in setYrs:
        ic(sy)
        if 'PreindustrialControl' in darr.scenario:
            timeSliceCesm = slice(
                cftime.DatetimeNoLeap(sy[0], 6, 29, 6, 0, 0, 0),
                cftime.DatetimeNoLeap(sy[1], 6, 29, 6, 0, 0, 0))
        else:
            timeSliceCesm = slice(
                cftime.DatetimeNoLeap(sy[0], 7, 15, 12, 0, 0, 0),
                cftime.DatetimeNoLeap(sy[1], 7, 15, 12, 0, 0, 0))
        timeSliceUkesm = slice(
            cftime.Datetime360Day(sy[0], 6, 30, 0, 0, 0, 0),
            cftime.Datetime360Day(sy[1], 6, 30, 0, 0, 0, 0))    
        try:
            darrToi = darr.sel(time=timeSliceCesm)
        except: # This will get messy when another model runs ARISE :)
            darrToi = darr.sel(time=timeSliceUkesm)
        darrToiItvlList.append(darrToi)
    
    return darrToiItvlList

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

def calc_spatial_grad_comp(darr):
    ''' Calculate spatial gradient components '''
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
            nsGradRlzList.append(
                sndimg.sobel( # North/south spatial gradient
                    darr.isel(realization=r), axis=0, mode='reflect'))
            ewGradRlzList.append(
                sndimg.sobel( # East/west spatial gradient
                    darr.isel(realization=r), axis=1, mode='reflect'))
        nsGrad = np.nanmean(nsGradRlzList, axis=0)
        ewGrad = np.nanmean(ewGradRlzList, axis=0)
    else:
        nsGrad = sndimg.sobel(darr, axis=1, mode='reflect')
        ewGrad = sndimg.sobel(darr, axis=2, mode='reflect')
    nsGradSc = nsGrad / (8 * distY)
    ewGradSc = ewGrad / (8 * distX)

    return nsGradSc, ewGradSc
    
def calc_spat_temp_grad(darrToi, setDict):
    ''' Calculate spatial and temporal gradients of temperature
    (common tasks for several variables) '''
    yrSpan = [
            darrToi.time.dt.year.data[0],
            darrToi.time.dt.year.data[-1]]
    # ic(yrSpan)
    tGradDict = calc_temporal_grad(darrToi, years=yrSpan)
    tGrad = tGradDict["grad"].compute()
    darrTimeMn = darrToi.mean(dim='time')
    nsGrad, ewGrad = calc_spatial_grad_comp(darrTimeMn)
    
    return tGrad, nsGrad, ewGrad, yrSpan

def calc_climate_speed(darr, setDict):
    ''' Calculate the climate speed of a variable '''
    darrIntvlList = time_slice_darr(darr, setDict)
    climSpdList = list()
    for darr in darrIntvlList:
        tGrad, nsGrad, ewGrad, yrSpan = calc_spat_temp_grad(darr, setDict)
        totGrad = np.sqrt((nsGrad ** 2) + (ewGrad ** 2))
        sGrad = totGrad
        climSpd = tGrad / sGrad
        climSpd.attrs = darr.attrs
        climSpdList.append(climSpd)

    climSpdItvl = xr.concat(climSpdList, dim='interval')
    # TODO: make attributes flexible by variable (e.g., climate speed of X)
    climSpdItvl.attrs['long_name'] = 'Climate speed of 2m temperature' \
        + ' ' + str(yrSpan[0]) + '-' + str(yrSpan[1])
    climSpdItvl.attrs['units'] = 'km/yr'
    
    # Display for troubleshooting
    # ic(climSpdItvl.scenario)
    # ic(check_stats(tGrad.data))
    # ic(check_stats(sGrad.data))
    # ic(check_stats(climSpdItvl.data))

    return climSpdItvl
    
def calc_warming_rate(darr, setDict):
    ''' Calculate the warming rate '''
    darrIntvlList = time_slice_darr(darr, setDict)
    tGradList = list()
    for darr in darrIntvlList:
        tGrad, nsGrad, ewGrad, yrSpan = calc_spat_temp_grad(darr, setDict)
        tGradList.append(tGrad)
        tGrad.attrs['scenario'] = darr.scenario

    tGradItvl = xr.concat(tGradList, dim='interval')
    # TODO: make attributes flexible by variable (e.g., climate speed of X)
    tGradItvl.attrs['long_name'] = 'Temporal gradient of 2m temperature' \
        + ' ' + str(yrSpan[0]) + '-' + str(yrSpan[1])
    tGradItvl.attrs['units'] = 'degC/yr'

    return tGradItvl
    
def calc_spat_grad(darr, setDict):
    ''' Calculate the spatial gradient of a variable '''
    tGrad, nsGrad, ewGrad, yrSpan = calc_spat_temp_grad(darr, setDict)  
    totGrad = np.sqrt((nsGrad ** 2) + (ewGrad ** 2))
    sGradDs = xr.Dataset(
        data_vars = dict(
            ns=(["lat","lon"], nsGrad),
            ew=(["lat","lon"], ewGrad),
        ),
        coords={
            "lat": tGrad['lat'],
            "lon": tGrad['lon']
        }
    )
    sGradDs.attrs = darr.attrs
    sGradDs['ns'].attrs['long_name'] = 'Spatial grad of 2m temp (N-S)'
    sGradDs['ns'].attrs['scenario'] = darr.scenario
    sGradDs['ns'].attrs['units'] = 'C/km'
    sGradDs['ew'].attrs['long_name'] = 'Spatial grad of 2m temp (E-W)'
    sGradDs['ew'].attrs['scenario'] = darr.scenario
    sGradDs['ew'].attrs['units'] = 'C/km'
    sGradDs.attrs['scenario'] = darr.scenario

    sGradDsToPlot = sGradDs['ew']
    
    return sGradDsToPlot

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

def calc_climate_velocity(darr, setDict):
    ''' Calculate the climate velocity (vector form) of a variable '''
    tGrad, nsGrad, ewGrad, yrsSpan = calc_spat_temp_grad(darr, setDict)
    nsClimVel = -1 * tGrad / nsGrad
    ewClimVel = -1 * tGrad / ewGrad
    # ewGrad, nsGrad
    totGrad = np.sqrt((nsGrad ** 2) + (ewGrad ** 2))
    sGrad = totGrad
    climSpd = tGrad / sGrad
    climVelDs = xr.Dataset(
        data_vars = dict(
            ns=(["realization","lat","lon"], nsClimVel.data),
            ew=(["realization","lat","lon"], ewClimVel.data),
            cspd=(["realization","lat","lon"], climSpd.data)
        ),
        coords={
            "realization": darr['realization'].data,
            "lat": darr['lat'],
            "lon": darr['lon']
        }
    )
       
    # climVelDs = xr.Dataset(
    #     data_vars = dict(
    #         ns=(["lat","lon"], nsClimVel.data),
    #         ew=(["lat","lon"], ewClimVel.data),
    #         cspd=(["lat","lon"], climSpd.data)
    #     ),
    #     coords={
    #         "lat": tGrad['lat'].data,
    #         "lon": tGrad['lon'].data
    #     }
    # )
    climVelDs['ns'].attrs['long_name'] = 'Climate velocity of 2m temp (N-S)'
    climVelDs['ns'].attrs['scenario'] = darr.scenario
    climVelDs['ns'].attrs['units'] = 'km/yr'
    climVelDs['ew'].attrs['long_name'] = 'Climate velocity of 2m temp (E-W)'
    climVelDs['ew'].attrs['scenario'] = darr.scenario
    climVelDs['ew'].attrs['units'] = 'km/yr'
    climVelDs['cspd'].attrs['long_name'] = 'Climate speed of 2m temp'
    climVelDs['cspd'].attrs['scenario'] = darr.scenario
    climVelDs['cspd'].attrs['units'] = 'km/yr'
    climVelDs.attrs['scenario'] = darr.scenario

    return climVelDs

                                                                                                            
# def calc_climate_velocity(darr, setDict):
#     ''' Calculate the climate velocity (vector form) of a variable '''
#     tGrad, nsGrad, ewGrad, yrsSpan = calc_spat_temp_grad(darr, setDict)
#     nsClimVel = tGrad / nsGrad
#     ewClimVel = tGrad / ewGrad
#     # ewGrad, nsGrad
#     totGrad = np.sqrt((nsGrad ** 2) + (ewGrad ** 2))
#     sGrad = totGrad   
#     climSpd = tGrad / sGrad
    
#     climVelDs = xr.Dataset(
#         data_vars = dict(
#             ns=(["realization","lat","lon"], nsClimVel.data),
#             ew=(["realization","lat","lon"], ewClimVel.data),
#             cspd=(["realization","lat","lon"], climSpd.data)
#         ),
#         coords={
#             "realization": tGrad['realization'].data,
#             "lat": tGrad['lat'],
#             "lon": tGrad['lon']
#         }
#     )
#     climVelDs['ns'].attrs['long_name'] = 'Climate velocity of 2m temp (N-S)'
#     climVelDs['ns'].attrs['scenario'] = darr.scenario
#     climVelDs['ns'].attrs['units'] = 'km/yr'
#     climVelDs['ew'].attrs['long_name'] = 'Climate velocity of 2m temp (E-W)'
#     climVelDs['ew'].attrs['scenario'] = darr.scenario
#     climVelDs['ew'].attrs['units'] = 'km/yr'
#     climVelDs['cspd'].attrs['long_name'] = 'Climate speed of 2m temp'
#     climVelDs['cspd'].attrs['scenario'] = darr.scenario
#     climVelDs['cspd'].attrs['units'] = 'km/yr'
#     climVelDs.attrs['scenario'] = darr.scenario

#     return climVelDs
                                    
def calc_seasonal_shift(darr, setDict):
    ''' Calculate seasonal shift of a variable '''
    mnthsForMn = 3
    
    if 'CESM2' in darr.scenario:
        setYrs = setDict["calcIntvl"]["CESM2-ARISE"]
    elif 'UKESM-ARISE' in darr.scenario:
        setYrs = setDict["calcIntvl"]["UKESM-ARISE"]
    elif 'GLENS' in darr.scenario:
        setYrs = setDict["calcIntvl"]["GLENS"]
    timeSliceCesm = slice(
        cftime.DatetimeNoLeap(setYrs[0], 1, 1, 0, 0, 0, 0),
        cftime.DatetimeNoLeap(setYrs[1], 12, 31, 0, 0, 0, 0))
    timeSliceUkesm = slice(
        cftime.Datetime360Day(setYrs[0], 1, 1, 0, 0, 0, 0),
        cftime.Datetime360Day(setYrs[1], 12, 25, 0, 0, 0, 0))
    try:
        darrToi = darr.sel(time=timeSliceCesm)
    except: # This will get messy when another model runs ARISE :)
        darrToi = darr.sel(time=timeSliceUkesm)
    
    val = 3 # Month for which to calculate seasonal shift
    monthsDarr = darr.groupby("time.month").mean()
    preDarr = monthsDarr.sel(month=val-1)
    postDarr = monthsDarr.sel(month=val+1)
    rateOfChange = (postDarr - preDarr) / 2
    marchDarr = darr.sel(time=is_val(darr['time.month'], val))
    tGradDict = calc_temporal_grad(marchDarr, years=setYrs)
    tGrad = tGradDict["grad"].compute() 
    
    ssnShiftMnYr = tGrad / rateOfChange # months/year
    ssnShift = (ssnShiftMnYr * 10 * 365.25) / 12 # days/decade
    ssnShift.attrs = darr.attrs
    ssnShift.attrs['long_name'] = 'Seasonal shift of 2m temperature ' + str(val) \
        + ' ' + str(setYrs[0]) + '-' + str(setYrs[1])
    ssnShift.attrs['units'] = 'days/decade'
    ic(check_stats(ssnShift))
    return ssnShift
    
def calc_area_exposed(cSpdDarr, setDict, dataDict, threshold):
    ''' Calculate area exposed to certain climate velocity '''
    cSpdData = cSpdDarr.data # Must be dimension lat x lon, no realization or interval
    cellAreaDset = fdd.generate_gridcellarea(saveFlag=False)
    cellAreaDarr = cellAreaDset['grid_cell_area']
    cesmMask = fpd.get_mask(dataDict, cSpdDarr.scenario)
    if dataDict["landmaskFlag"] == 'land':
        loMask = cesmMask > 0
    elif dataDict["landmaskFlag"] == 'ocean':
        loMask = cesmMask == 0
    else:
        loMask = cesmMask > -1
    cellAreaLandmask = cellAreaDarr.copy()
    cellAreaLandmask.data[~loMask] = np.nan
    cellAreaData = cellAreaLandmask.data
    cSpdDataAbs = abs(cSpdData)
    # ic(check_stats(cSpdDataAbs))
    velMask = cSpdDataAbs >= threshold
    # ic(velMask, loMask)
    # ic(np.array_equal(velMask, loMask))
    # sys.exit('STOP')
    areaExposed = cellAreaLandmask.copy()
    areaExposed.data[~velMask] = np.nan
    
    totalLandmaskArea = cellAreaLandmask.sum(skipna=True)
    totalAreaExposed = np.nansum(areaExposed)
    percAreaExposed = (totalAreaExposed / totalLandmaskArea) * 100

    return percAreaExposed

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

def is_val(month, val):
    ''' Thanks to stackoverflow.com/a/40272331 for this idea! '''
    mask = month == val
    return mask
    
def calc_weighted_med(data):
    ''' Calculated latitude-weighted global median of an array '''
    latWeights = np.cos(np.deg2rad(data['lat']))
    darrWght = data.weighted(latWeights)
    ensMed = darrWght.quantile(0.5, dim=('lat','lon'), skipna=True)
    return ensMed