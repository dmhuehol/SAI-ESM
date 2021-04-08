import xarray as xr
import numpy as np
import sys
from cdo import *

# baselineFname = baselineStr[len(baselinePath):len(baselineStr)] #only filename

def prep_raw_data(inPath,inCard,outFile):
# Conduct basic preparation of raw data to make it suitable to run through other
    # functions. Current list of tasks: merge decadal files, shift time -1 days,
    # annual average
    cdo = Cdo()
    inFile = inPath + inCard
    mergeFile = inPath + 'merge_' + outFile + '.nc'
    shiftFile = inPath + 'shift_' + outFile + '.nc'
    yearFile = inPath + 'annual_' + outFile + '.nc'

    cdo.mergetime(input=inFile, output=mergeFile)
    print('Merged')
    cdo.shifttime('-1days', input=mergeFile, output=shiftFile)
    print('Shifted')
    cdo.yearmonmean(input=shiftFile, output=yearFile)
    print('Annual averaged')
    print('Completed')

    return


def extract_meta_from_fname(fname):
# Extract useful metadata from GLENS output filename
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
# Extract time intervals (usually decades) from GLENS data
    for cdc,cdv in enumerate(intervalsToPlot):
        startInd = np.where(years==cdv)[0][0]
        if cdv+timePeriod == 2100:
            endInd = len(years)
        else:
            endInd = np.where(years==cdv+timePeriod)[0][0]
        # handlesToPlot[cdv]["data"] = darr[startInd:endInd]
        handlesToPlot.append(darr[startInd:endInd])

    return handlesToPlot
