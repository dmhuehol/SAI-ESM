''' process_glens_fun
Contains functions to process GLENS data, taking care of operations like
extracting metadata from the filename, choosing and summing levels, or
managing area inputs.

Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import xarray as xr
import numpy as np
from cdo import *

def prep_raw_data(inPath,outPath,inCard,outFile):
    ''' Conduct basic preparation of raw data to make it suitable to run through other
    functions. Current list of tasks: merge decadal files, shift time -1 days,
    annual average. Requires Python CDO bindings, and thus cannot be run on the
    NCAR system.
    '''

    cdo = Cdo()
    inFile = inPath + inCard
    mergeFile = outPath + 'merge_' + outFile + '.nc'
    shiftFile = outPath + 'shift_' + outFile + '.nc'
    yearFile = outPath + 'annual_' + outFile + '.nc'

    cdo.mergetime(input=inFile, output=mergeFile)
    print('Merged')
    cdo.shifttime('-1days', input=mergeFile, output=shiftFile)
    print('Shifted')
    cdo.yearmonmean(input=shiftFile, output=yearFile)
    print('Annual averaged')
    print('Completed')

    return


def extract_meta_from_fname(fname):
    ''' Extract useful metadata from GLENS output filename '''

    # baselineFname = fname[len(fname):len(fname)] #only filename
    fnamePcs = fname.split('.') #period is delimiter in filename
    glensMeta = {"unkn1": fnamePcs[0],
        "unkn2": fnamePcs[1],
        "unkn3": fnamePcs[2],
        "unkn4": fnamePcs[3],
        "runType": fnamePcs[4],
        "ensMemNum": int(fnamePcs[5]),
        "unkn5": fnamePcs[6],
        "unkn6": fnamePcs[7],
        "dataVar": fnamePcs[8],
        "rawDates": fnamePcs[9],
        "strtYr": int(fnamePcs[9][0:4]),
        "strtMnth": int(fnamePcs[9][4:6]),
        "endYr": int(fnamePcs[9][7:11]),
        "endMnth": int(fnamePcs[9][11:13]),
        "extension": fnamePcs[10]
    }

    return glensMeta

def discover_data_var(glensDsetCntrl):
    ''' GLENS files contain lots of variables--this finds the one with actual data! '''

    fileKeys = list(glensDsetCntrl.keys())
    notDataKeys = ['time_bnds', 'date', 'datesec']
    dataKey = "empty"
    for cKey in fileKeys:
        if cKey in notDataKeys:
            print('not this one: ' + cKey)
        else:
            dataKey = cKey
            print('Data key is: ' + dataKey)
            if dataKey == "empty":
                sys.exit("Data discovery failed! Ending run now...")

    return dataKey

def find_matching_year_bounds(glensCntrl,glensFdbck):
    ''' Find indices that bound matching time periods in control and feedback data '''

    cntrlYrs = glensCntrl['time'].dt.year.data
    fdbckYrs = glensFdbck['time'].dt.year.data
    bothYrs,cntrlInd,fdbckInd = np.intersect1d(cntrlYrs, fdbckYrs, return_indices=True)
    firstYrInBoth = bothYrs[0]
    lastYrInBoth = bothYrs[len(bothYrs)-1]
    bndDct = {
        "strtYrMtch": firstYrInBoth,
        "endYrMtch": lastYrInBoth,
        "cntrlStrtMtch": cntrlInd[0],
        "fdbckStrtMtch": fdbckInd[0],
        "cntrlEndMtch": cntrlInd[len(cntrlInd)-1],
        "fdbckEndMtch": fdbckInd[len(fdbckInd)-1],
        "mtchYrs": bothYrs}

    return bndDct


def extract_doi(intervalsToPlot, years, timePeriod, darr, handlesToPlot):
    ''' Extract time intervals (usually decades) from GLENS data '''

    for cdc,cdv in enumerate(intervalsToPlot):
        startInd = np.where(years==cdv)[0][0]
        if cdv+timePeriod == 2100:
            endInd = len(years)
        else:
            endInd = np.where(years==cdv+timePeriod)[0][0]
        # handlesToPlot[cdv]["data"] = darr[startInd:endInd]
        handlesToPlot.append(darr[startInd:endInd])

    return handlesToPlot

def find_closest_level(darr, levOfInt, levName='lev'):
    ''' Find the index of the level which is closest to an input value '''

    levs = darr[levName].data
    levDiffs = np.abs(levs - levOfInt)
    closestVal = np.min(levDiffs)
    closestTup = np.where(levDiffs == closestVal)

    if len(closestTup) > 1:
        sys.exit('Conflicting levels! Input new level of interest and try again.')

    closestInd = closestTup[0][0] #np.where returns tuple, which can't be indexed

    return closestInd

def obtain_levels(darr, levOfInt, levName='lev'):
    ''' Obtain relevant level(s) from data based on levOfInt input '''

    levs = darr[levName].data
    if levOfInt == 'total':
        darr = darr.sum(dim=levName)
    elif levOfInt == 'troposphere':
        indTpause = find_closest_level(darr, 200, levName=levName) #simple split on 200hPa for now
        levMask = levs > levs[indTpause]
        darr = darr[:,levMask,:,:]
        darr = darr.sum(dim=levName)
    elif levOfInt == 'stratosphere':
        indTpause = find_closest_level(darr, 200, levName=levName) #simple split on 200hPa for now
        levMask = levs <= levs[indTpause]
        darr = darr[:,levMask,:,:]
        darr = darr.sum(dim=levName)
    elif np.size(levOfInt)==2:
        indHghrPres = find_closest_level(darr, max(levOfInt), levName=levName)
        indLwrPres = find_closest_level(darr, min(levOfInt), levName=levName)
        levOfInt = [levs[indLwrPres],levs[indHghrPres]]
        darr = darr[:,indLwrPres:indHghrPres+1,:,:]
        darr = darr.sum(dim=levName)
    else:
        indClosest = find_closest_level(darr, levOfInt, levName=levName)
        darr.attrs[levName] = darr[levName].data[indClosest]
        darr = darr.sel(lev=darr.attrs[levName])

    return darr

def make_level_string(darr, levOfInt):
    ''' Create strings describing level of data used in title and filenames. '''

    if isinstance(levOfInt,str):
        levStr = levOfInt
    elif np.size(levOfInt)==2:
        if (np.round_(levOfInt[0],decimals=1)==0) | (np.round_(levOfInt[1],decimals=1)==0):
            levStr = str(np.round_(levOfInt,decimals=6)) + ' mb'
        else:
            levStr = str(np.round_(levOfInt,decimals=1)) + ' mb'
    elif np.round_(darr.attrs['lev'],decimals=1) == 0:
        levStr = str(np.round_(darr.attrs['lev'],decimals=6)) + ' mb'
    else:
        levStr = str(np.round_(darr.attrs['lev'],decimals=1)) + ' mb'

    return levStr

def manage_area(darr, regionToPlot, areaAvgBool=True):
    ''' Manage area operations: obtain global, regional, or pointal output '''

    if regionToPlot == 'global':
        ic('global')
        locStr = 'global'
        locTitleStr = 'global'

        if areaAvgBool:
            latWeights = np.cos(np.deg2rad(darr['lat']))
            darrWght = darr.weighted(latWeights)
            darr = darrWght.mean(dim=['lat','lon'])

    elif isinstance(regionToPlot,dict):
        ic('regional')
        locStr = regionToPlot['regSaveStr']
        locTitleStr = regionToPlot['regSaveStr']

        lats = darr['lat'] #feedback and control are on same grid, fortunately
        lons = darr['lon']
        latMask = (lats>regionToPlot['regLats'][0]) & (lats<regionToPlot['regLats'][1])
        lonMask = (lons>regionToPlot['regLons'][0]) & (lons<regionToPlot['regLons'][1])
        darrBoxMask = darr[:,latMask,lonMask]

        if areaAvgBool:
            darr = darrBoxMask.mean(dim=['lat','lon'])

    elif isinstance(regionToPlot,list):
        ic('point')
        darr = darr.sel(lat=regionToPlot[0], lon=regionToPlot[1], method="nearest")

        latStr = str(np.round_(darr.lat.data,decimals=2))
        lonStr = str(np.round_(darr.lon.data,decimals=2))
        locStr = latStr + '_' + lonStr
        locTitleStr = '(' + latStr + ',' + lonStr + ')'
    else:
        sys.exit('Invalid region! Check value for regionToPlot.')

    return darr, locStr, locTitleStr

def isolate_change_quantile(darr, quantileOfInt):
    ''' Isolate the largest change by a quantile cutoff '''
    darrAbs = np.abs(darr)
    normValue = np.max(darrAbs)
    darrAbsNorm = darrAbs / normValue
    quantCut = darrAbsNorm.quantile(quantileOfInt)
    darrAbsNormQ = darrAbsNorm
    darrAbsNormQ.data[darrAbsNormQ.data < quantCut.data] = np.nan

    return darrAbsNormQ
