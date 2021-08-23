''' process_glens_fun
Contains functions to process GLENS data, taking care of operations like
extracting metadata from the filename, choosing and summing levels, or
managing area inputs.

Written by Daniel Hueholt | July 2021
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

import matplotlib.pyplot as plt

def discover_data_var(glensDsetCntrl):
    ''' GLENS files contain lots of variables--this finds the one with actual data! '''

    fileKeys = list(glensDsetCntrl.keys())
    notDataKeys = ['time_bnds', 'date', 'datesec', 'lev_bnds', 'gw', 'ch4vmr',
                   'co2vmr', 'ndcur', 'nscur', 'sol_tsi', 'nsteph', 'f11vmr',
                   'n2ovmr', 'f12vmr', 'lon_bnds', 'lat_bnds']
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

# def find_matching_year_bounds(rlzList):
#     ''' Find indices that bound matching time periods in control and feedback output '''
# # NEED TO FIX
#     rlzYearList = list()
#     for rlz in rlzList:
#     rlzYears = rlz['time'].dt.year.data
#     rlzYearList.append(rlz)
#
#     if rlzList[0].ndim == 4: #rlzList[0]=it doesn't matter which scenario, and index 0 will always be present
#         nanYrInd = np.where(np.isnan(rlzList[0][:,1,1,1])) #Years after the model completes may be present with NaN values
#     elif rlzList[0].ndim == 3: #time,lat,lon
#         nanYrInd = np.where(np.isnan(rlzList[0][:,1,1])) #Years after the model completes may be present with NaN values
#     else:
#         ic(rlzList[0].ndim)
#         # Not giving a sys.exit() here, want whatever the error is to show up
#
#
#     bothYrs,cntrlInd,fdbckInd = np.intersect1d(cntrlYrs, fdbckYrs, return_indices=True) #Or they may be not present at all
#     firstYrInBoth = bothYrs[0]
#     try:
#         lastYrInBoth = np.min(np.array([bothYrs[len(bothYrs)-1],cntrlYrs[nanYrInd[0][0]-1]])) #Either way, the last year is whichever method with the earliest end date
#     except:
#         lastYrInBoth = bothYrs[len(bothYrs)-1]
#     bndDct = {
#         "strtYrMtch": firstYrInBoth,
#         "endYrMtch": lastYrInBoth,
#         "cntrlStrtMtch": cntrlInd[0],
#         "fdbckStrtMtch": fdbckInd[0],
#         "cntrlEndMtch": cntrlInd[len(cntrlInd)-1],
#         "fdbckEndMtch": fdbckInd[len(fdbckInd)-1],
#         "mtchYrs": bothYrs}
#
#     return bndDct


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
    if levOfInt == None:
        return darr

    if 'plev' in darr.dims: #CMIP6
        levName = 'plev'
        levOfInt = levOfInt * 100

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
        try:
            darr = darr.sel(lev=darr.attrs[levName])
        except: #This is horrible
            darr = darr.sel(plev=darr.attrs[levName])

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

def open_data(dataDict, setDict):
    ''' Opens data and select data variable '''
    if '*' in dataDict["fnameCntrl"]:
        cntrlPath = dataDict["dataPath"] + dataDict["fnameCntrl"]
        fdbckPath = dataDict["dataPath"] + dataDict["fnameFdbck"]
        glens2Path = dataDict["dataPath"] + dataDict["fnameGlens2"]
        cmip6Path = dataDict["dataPath"] + dataDict["fnameCmip6"]
        glensDsetCntrl = xr.open_mfdataset(cntrlPath, concat_dim='realization', combine='nested')
        glensDsetFdbck = xr.open_mfdataset(fdbckPath, concat_dim='realization', combine='nested')
        try:
            glens2Dset = xr.open_mfdataset(glens2Path, concat_dim='realization', combine='nested')
        except:
            ic('No GLENS2 files located')
            glens2Dset = None
        try:
            cmip6Dset = xr.open_mfdataset(cmip6Path, concat_dim='realization', combine='nested')
        except:
            ic('No CMIP6 files located')
            cmip6Dset = None
        if setDict['landmaskFlag'] == 'land':
            glensDsetCntrl = glensDsetCntrl.where(glensDsetCntrl.landmask > 0)
            glensDsetFdbck = glensDsetFdbck.where(glensDsetFdbck.landmask > 0)
            glens2Dset = glens2Dset.where(glens2Dset.landmask > 0)
            cesmMask = xr.open_dataset('/Users/dhueholt/Documents/Summery_Summary/daniel_mask.nc')
            cmip6Dset = cmip6Dset.where(cesmMask.imask>0)
        dataKey = discover_data_var(glensDsetCntrl)
        glensDarrCntrl = glensDsetCntrl[dataKey]
        glensDarrCntrl.attrs['scenario'] = 'GLENS1:Control/RCP8.5'
        glensDarrFdbck = glensDsetFdbck[dataKey]
        glensDarrFdbck.attrs['scenario'] = 'GLENS1:Feedback/SAI/G1.2[8.5]'
        if glens2Dset is not None:
            glens2Darr = glens2Dset[dataKey]
            glens2Darr.attrs['scenario'] = 'GLENS2:Feedback/SAI/G1.5[4.5]'
        else:
            glens2Darr = None
        if cmip6Dset is not None:
            dataKeyCmip6 = discover_data_var(cmip6Dset) #CMIP6 has unique variable names
            cmip6Darr = cmip6Dset[dataKeyCmip6]
            cmip6Darr.attrs['scenario'] = 'CMIP6:Control/SSP2-4.5'
        else:
            cmip6Darr = None
    else:
        sys.exit("Check input! Token should have a wildcard (i.e. match multiple files).")

    darrCheckList = list([glensDarrCntrl,glensDarrFdbck,glens2Darr,cmip6Darr])
    darrList = list()
    for d in darrCheckList: #Surely there's a better way
        if d is None:
            pass
        else:
            darrList.append(d)

    return darrList, dataKey

def get_ens_mem(files):
    ''' Find ensemble members from a list of files '''
    emem = list()
    for pth in files:
        pthPcs = pth.split('/')
        fName = pthPcs[len(pthPcs)-1]
        fNamePcs = fName.split('_')
        emem.append(fNamePcs[1])

    return emem

def manage_realizations(setDict, darr, emem):
    ''' Either obtain realization of interest or calculate ensemble mean, and
    create relevant filename '''
    try:
        if 'GLENS1:Control' in darr.scenario:
            scnStr = 'g1c' #glens1control
        elif 'GLENS1:Feedback' in darr.scenario:
            scnStr = 'g1f' #glens1feedback
        elif 'GLENS2:Feedback' in darr.scenario:
            scnStr = 'g2f' #glens2feedback
        elif 'CMIP6:Control' in darr.scenario:
            scnStr = 'c6c' #cmip6control

        if setDict['realization'] == 'mean': #Output DataArray of ensemble mean
            darrMn = darr.mean(dim='realization')
            darrOut = darrMn.compute()
            ememSave = 'mn' + scnStr
        elif setDict['realization'] == 'ensplot': #Output DataArray of all members and ensemble mean
            darrMn = darr.mean(dim='realization')
            darrOut = xr.concat([darr,darrMn],dim='realization') #Add ensemble mean as another "realization"
            ememSave = 'ens' + scnStr
        else: #Output DataArray of single ensemble member
            ememNum = list(map(int, emem))
            rInd = ememNum.index(setDict['realization'])
            darrOut = darr[rInd,:,:,:].compute()

            activeEmem = emem[rInd]
            ememSave = scnStr + activeEmem
    except:
        darrOut = None
        ememSave = ''

    return darrOut, ememSave

def call_to_open(dataDict, setDict):
    ''' Common data tasks for all basic plots '''
    darrList, dataKey = open_data(dataDict, setDict)

    cntrlFiles = sorted(glob.glob(dataDict['dataPath'] + dataDict['fnameCntrl']))
    fdbckFiles = sorted(glob.glob(dataDict['dataPath'] + dataDict['fnameFdbck']))
    glens2Files = sorted(glob.glob(dataDict['dataPath'] + dataDict['fnameGlens2']))
    cmip6Files = sorted(glob.glob(dataDict['dataPath'] + dataDict['fnameCmip6']))
    ememList = list()
    for ec in (cntrlFiles, fdbckFiles, glens2Files, cmip6Files):
        ememList.append(get_ens_mem(ec))

    rlzList = list()
    ememStrList = list()
    for dc,dv in enumerate(darrList):
        rlzArr, ememStr = manage_realizations(setDict, dv, ememList[dc])
        rlzList.append(rlzArr)
        ememStrList.append(ememStr)

    ememStrList = list(filter(None,ememStrList))
    ememSave = '-'.join(ememStrList)
    cmnDict = {'dataKey': dataKey, 'ememSave': ememSave}

    if setDict["convert"] is not None:
        for cnvrtr in setDict["convert"]:
            for rc,rv in enumerate(rlzList):
                rlzList[rc] = cnvrtr(rv)

    return rlzList, cmnDict

def meta_book(setDict, dataDict, cntrlToPlot, labelsToPlot=None):
    ''' Compile the bits and pieces used in filenames and titles '''
    metaDict = {
        "cntrlStr": 'RCP8.5',
        "ssp245Str": 'SSP2-4.5',
        "fdbckStr": 'G1.2(8.5)',
        "fdbckStrG2": 'G1.5(2-4.5)',
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
        "yStr": cntrlToPlot.units,
        "unit": cntrlToPlot.attrs['units'],
        "pdfStyle": setDict["plotStyle"],
        "spcStr": make_spc_string(setDict),
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

def make_spc_string(setDict):
    ''' Make string reflecting area average Boolean '''
    try:
        if setDict["areaAvgBool"]:
            spcStr = 'spcavg'
        else:
            spcStr = 'nospcavg'
    except:
        spcStr = ''

    return spcStr

def apply_landmask(dset):
    ''' Apply embedded landmask to GLENS or SCIRIS output '''
    dset['landmask'] = dset['landmask'].fillna(0) #Change nan values to 0 for consistency


    return dsetMasked
