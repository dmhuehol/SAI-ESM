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
xr.set_options(keep_attrs=True)
import dask
dask.config.set(**{'array.slicing.split_large_chunks': True})
import numpy as np
import matplotlib.path as mpth
import glob

def discover_data_var(glensDsetCntrl):
    ''' GLENS files contain lots of variables--this finds the one with actual data! '''

    fileKeys = list(glensDsetCntrl.keys())
    notDataKeys = ['time_bnds', 'date', 'datesec', 'lev_bnds', 'gw', 'ch4vmr',
                   'co2vmr', 'ndcur', 'nscur', 'sol_tsi', 'nsteph', 'f11vmr',
                   'n2ovmr', 'f12vmr']
    notDataInDset = list()
    dataKey = "empty"

    for cKey in fileKeys:
        if cKey in notDataKeys:
            notDataInDset.append(cKey)
        else:
            dataKey = cKey
            print('Data key could be: ' + dataKey)

    if dataKey == "empty":
        sys.exit("Data discovery failed! Ending run now...")

    return dataKey

def find_matching_year_bounds(glensCntrl,glensFdbck):
    ''' Find indices that bound matching time periods in control and feedback output '''

    cntrlYrs = glensCntrl['time'].dt.year.data
    fdbckYrs = glensFdbck['time'].dt.year.data
    nanYrInd = np.where(np.isnan(glensCntrl[:,1,1,1])) #Years after the model completes may be present with NaN values
    bothYrs,cntrlInd,fdbckInd = np.intersect1d(cntrlYrs, fdbckYrs, return_indices=True) #Or they may be not present at all
    firstYrInBoth = bothYrs[0]
    try:
        lastYrInBoth = np.min(np.array([bothYrs[len(bothYrs)-1],cntrlYrs[nanYrInd[0][0]-1]])) #Either way, the last year is whichever method with the earliest end date
    except:
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
        if np.ndim(darr) == 3:
            darr = darr[levMask,:,:]
            darr = darr.sum(dim=levName)
        else:
            darr = darr[:,levMask,:,:]
            darr = darr.sum(dim=levName)
    elif levOfInt == 'stratosphere':
        indTpause = find_closest_level(darr, 200, levName=levName) #simple split on 200hPa for now
        levMask = levs <= levs[indTpause]
        if np.ndim(darr) == 3:
            darr = darr[levMask,:,:]
            darr = darr.sum(dim=levName)
        else:
            darr = darr[:,levMask,:,:]
            darr = darr.sum(dim=levName)
    elif np.size(levOfInt)==2:
        indHghrPres = find_closest_level(darr, max(levOfInt), levName=levName)
        indLwrPres = find_closest_level(darr, min(levOfInt), levName=levName)
        levOfInt = [levs[indLwrPres],levs[indHghrPres]]
        if np.ndim(darr) == 3:
                    darr = darr[indLwrPres:indHghrPres+1,:,:]
                    darr = darr.sum(dim=levName)
        else:
            darr = darr[:,indLwrPres:indHghrPres+1,:,:]
            darr = darr.sum(dim=levName)
    else:
        indClosest = find_closest_level(darr, levOfInt, levName=levName)
        darr.attrs[levName] = darr[levName].data[indClosest]
        darr = darr.sel(lev=darr.attrs[levName])

    return darr

def make_level_string(darr, levOfInt):
    ''' Create strings describing level of data used in title and filenames. '''

    try:
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
    except:
        levStr = ''

    return levStr

def manage_area(darr, regionToPlot, areaAvgBool=True):
    ''' Manage area operations: obtain global, regional, or pointal output '''

    if regionToPlot == 'global':
        locStr = 'global'
        locTitleStr = 'global'

        if areaAvgBool:
            latWeights = np.cos(np.deg2rad(darr['lat']))
            darrWght = darr.weighted(latWeights)
            darr = darrWght.mean(dim=['lat','lon'], skipna=True)

    elif isinstance(regionToPlot,dict):
        locStr = regionToPlot['regSaveStr']
        locTitleStr = regionToPlot['regSaveStr']

        lats = darr['lat'] #feedback and control are on same grid, fortunately
        lons = darr['lon']
        if len(regionToPlot['regLons'])>2: #non-rectangular region that does not cross Prime Meridian
            gridMask = make_polygon_mask(lats, lons, regionToPlot['regLats'], regionToPlot['regLons'])
            darrBoxMask = darr.copy()
            darrBoxMask.data[:,~gridMask] = np.nan
        elif isinstance(regionToPlot['regLons'], tuple): #non-rectangular region that crosses Prime Meridian
            sGridMaskList = list()
            for sc in np.arange(0,len(regionToPlot['regLons'])):
                sGridMask = make_polygon_mask(lats, lons, regionToPlot['regLats'][sc], regionToPlot['regLons'][sc])
                sGridMaskList.append(sGridMask)
            gridMask = np.logical_or.reduce(sGridMaskList)
            darrBoxMask = darr.copy()
            darrBoxMask.data[:,~gridMask] = np.nan
        else: #rectangular region
            latMask = (lats>regionToPlot['regLats'][0]) & (lats<regionToPlot['regLats'][1])
            if regionToPlot['regLons'][0] < regionToPlot['regLons'][1]: #rectangle does not cross Prime Meridian
                lonMask = (lons>regionToPlot['regLons'][0]) & (lons<regionToPlot['regLons'][1])
            else: #rectangle crosses Prime Meridian
                lonMask = (lons>regionToPlot['regLons'][0]) | (lons<regionToPlot['regLons'][1])
            darrBoxMask = darr[:,latMask,lonMask]

        if areaAvgBool:
            darr = darrBoxMask.mean(dim=['lat','lon'], skipna=True)
        else:
            darr = darrBoxMask

    elif isinstance(regionToPlot,list):
        darr = darr.sel(lat=regionToPlot[0], lon=regionToPlot[1], method="nearest")

        latStr = str(np.round_(darr.lat.data,decimals=2))
        lonStr = str(np.round_(darr.lon.data,decimals=2))
        locStr = latStr + '_' + lonStr
        locTitleStr = '(' + latStr + ',' + lonStr + ')'
    else:
        sys.exit('Invalid region! Check value for regionToPlot.')

    return darr, locStr, locTitleStr

def isolate_change_quantile(darr, quantileOfInt):
    ''' Isolate the largest values by a quantile cutoff '''
    darrAbs = np.abs(darr)
    normValue = np.max(darrAbs)
    darrAbsNorm = darrAbs / normValue
    quantCut = darrAbsNorm.quantile(quantileOfInt)

    darrNorm = darr / np.nanmax(darr)
    # quantCut = darrNorm.quantile(quantileOfInt)
    # ic(quantCut.data) #troubleshooting
    quantMask = np.logical_and(darrNorm>-quantCut, darrNorm<quantCut)
    darrNorm.data[quantMask] = np.nan

    darrNorm.attrs['units'] = 'dimless'

    return darrNorm

def open_data(dataDict):
    ''' Opens data and select data variable '''
    if '*' in dataDict["fnameCntrl"]:
        cntrlPath = dataDict["dataPath"] + dataDict["fnameCntrl"]
        fdbckPath = dataDict["dataPath"] + dataDict["fnameFdbck"]
        glensDsetCntrl = xr.open_mfdataset(cntrlPath,concat_dim='realization',combine='nested')
        glensDsetFdbck = xr.open_mfdataset(fdbckPath,concat_dim='realization',combine='nested')
        dataKey = discover_data_var(glensDsetCntrl)
        glensDarrCntrl = glensDsetCntrl[dataKey]
        glensDarrFdbck = glensDsetFdbck[dataKey]
    else:
        cntrlPath = dataDict["dataPath"] + dataDict["fnameCntrl"]
        fdbckPath = dataDict["dataPath"] + dataDict["fnameFdbck"]
        glensDsetCntrl = xr.open_dataset(cntrlPath)
        glensDsetFdbck = xr.open_dataset(fdbckPath)
        dataKey = discover_data_var(glensDsetCntrl)
        glensDarrCntrl = glensDsetCntrl[dataKey]
        glensDarrFdbck = glensDsetFdbck[dataKey]

    return glensDarrCntrl, glensDarrFdbck, dataKey

def get_ens_mem(files):
    ''' Find ensemble members from a list of files '''
    emem = list()
    for pth in files:
        pthPcs = pth.split('/')
        fName = pthPcs[len(pthPcs)-1]
        fNamePcs = fName.split('_')
        emem.append(fNamePcs[1])

    return emem

def manage_realizations(setDict, cntrlDarr, fdbckDarr, ememCntrl, ememFdbck):
    ''' Either obtain realization of interest or calculate ensemble mean, and
    create relevant filename '''
    if setDict['realization'] == 'mean':
        cntrlDarrMn = cntrlDarr.mean(dim='realization')
        fdbckDarrMn = fdbckDarr.mean(dim='realization')
        cntrlDarrOut = cntrlDarrMn.compute()
        fdbckDarrOut = fdbckDarrMn.compute()
        ememSaveCntrl = 'mnc' + 'r'.join(ememCntrl)
        ememSaveFdbck = 'mnf' + 'r'.join(ememFdbck)
        ememSave = ememSaveCntrl + '-' + ememSaveFdbck
    else:
        ememCntrlNum = list(map(int, ememCntrl))
        rCntrl = ememCntrlNum.index(setDict['realization'])
        cntrlDarrOut = cntrlDarr[rCntrl,:,:,:].compute()
        ememFdbckNum = list(map(int, ememFdbck))
        rFdbck = ememFdbckNum.index(setDict['realization'])
        fdbckDarrOut = fdbckDarr[rFdbck,:,:,:].compute()

        activeCntrlEmem = ememCntrl[rCntrl]
        activeFdbckEmem = ememFdbck[rFdbck]
        ememSave = 'rc' + activeCntrlEmem + '-' + 'rf' + activeFdbckEmem

    return cntrlDarrOut, fdbckDarrOut, ememSave

def call_to_open(dataDict, setDict):
    ''' Common data tasks for all basic plots '''
    glensDarrCntrl, glensDarrFdbck, dataKey = open_data(dataDict)

    cntrlFiles = sorted(glob.glob(dataDict['dataPath'] + dataDict['fnameCntrl']))
    fdbckFiles = sorted(glob.glob(dataDict['dataPath'] + dataDict['fnameFdbck']))
    ememCntrl = get_ens_mem(cntrlFiles)
    ememFdbck = get_ens_mem(fdbckFiles)
    glensCntrlRlz, glensFdbckRlz, ememSave = manage_realizations(setDict, glensDarrCntrl, glensDarrFdbck, ememCntrl, ememFdbck)
    cmnDict = {'dataKey': dataKey, 'ememSave': ememSave}

    if setDict["convert"] is not None:
        glensCntrlRlz = setDict["convert"](glensCntrlRlz)
        glensFdbckRlz = setDict["convert"](glensFdbckRlz)

    return glensCntrlRlz, glensFdbckRlz, cmnDict

def meta_book(setDict, dataDict, cntrlToPlot, labelsToPlot=None):
    ''' Compile the bits and pieces used in filenames and titles '''
    metaDict = {
        "cntrlStr": 'RCP8.5',
        "fdbckStr": 'SAI',
        "varStr": cntrlToPlot.long_name,
        "varSve": cntrlToPlot.long_name.replace(" ",""),
        "strtStr": str(cntrlToPlot['time'].data[0].year),
        "endStr": str(cntrlToPlot['time'].data[len(cntrlToPlot)-1].year),
        "frstDcd": str(setDict["startIntvl"][0]) + '-' + str(setDict["startIntvl"][1]),
        "lstDcd": str(setDict["endIntvl"][0]) + '-' + str(setDict["endIntvl"][1]),
        "tmStr": bcf_parser(labelsToPlot),
        "levStr": make_level_string(cntrlToPlot, setDict["levOfInt"]) if "levOfInt" in setDict.keys() else '',
        "levSve": make_level_string(cntrlToPlot, setDict["levOfInt"]).replace(" ","") if "levOfInt" in setDict.keys() else '',
        "ensStr": dataDict["ememSave"],
        "qntlStr": str(setDict["quantileOfInt"]),
        "yStr": cntrlToPlot.units,
        "unit": cntrlToPlot.attrs['units'],
        "pdfStyle": setDict["plotStyle"],
        "spcStr": 'spcavg' if setDict["areaAvgBool"] else 'nospcavg',
        "pid": {'g1p': 'globe_1p', 'g4p': 'globe_4p', 'ts': 'timeseries', 'pdf': 'pdf'},
        "glbType": {'vGl': 'vertical', 'bGl': 'baseline', 'fcStr': 'FdbckCntrl'}
    }

    return metaDict

def bcf_parser(labelsToPlot):
    ''' Parse labels to make filenames '''

    timeStr = list()
    if labelsToPlot != None:
        for lab in labelsToPlot:
            if 'Baseline' in lab:
                timeStr.append('b' + lab[0:9].replace("-",""))
            elif 'RCP8.5' in lab:
                timeStr.append('c' + lab[0:9].replace("-",""))
            elif 'SAI' in lab:
                timeStr.append('f' + lab[0:9].replace("-",""))

    timeStr = "_".join(timeStr)

    return timeStr

def make_polygon_mask(lats, lons, regionLats, regionLons):
    ''' Make mask for a non-rectangular polygonal region '''
    gridLon,gridLat = np.meshgrid(lons,lats) #make 2D lonxlat grids of lon/lat
    flatLon = np.ravel(gridLon)
    flatLat = np.ravel(gridLat)
    flatLatLon = np.transpose(np.vstack((flatLat,flatLon))) #Nx2
    regionPoly = np.transpose(np.vstack((regionLats,regionLons))) #Polyx2

    regionPath = mpth.Path(regionPoly) #defines the path of the region
    flatMask = regionPath.contains_points(flatLatLon) #the heavyweight--calculates whether points from the grid are exterior to the path or not
    gridMask = flatMask.reshape((len(lats),len(lons)))

    return gridMask

def norm_by_absmax(darr):
    ''' Normalize data to within [-1,1] wrt its max magnitude '''
    darrAbs = np.abs(darr)
    normValue = np.max(darrAbs)
    darrNorm = darr / normValue

    darrNorm.attrs['units'] = 'dimless'

    return darrNorm
