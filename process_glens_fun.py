''' process_glens_fun
Contains functions for operations like extracting metadata from filenames,
choosing and summing levels, or managing area inputs.

Written by Daniel Hueholt | August 2021
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

def discover_data_var(dset):
    ''' Find the data variable from among the many variables in a dataset '''
    fileKeys = list(dset.keys())
    notDataKeys = ['time_bnds', 'date', 'datesec', 'lev_bnds', 'gw', 'ch4vmr',
                   'co2vmr', 'ndcur', 'nscur', 'sol_tsi', 'nsteph', 'f11vmr',
                   'n2ovmr', 'f12vmr', 'lon_bnds', 'lat_bnds', 'ZSOI', 'BSW',
                   'WATSAT', 'landmask', 'ZLAKE', 'DZLAKE', 'SUCSAT',
                   'landfrac', 'topo', 'DZSOI', 'area', 'pftmask', 'HKSAT',
                   'nstep', 'mdcur', 'mscur', 'mcdate', 'mcsec']
    notDataInDset = list()
    dataKey = None

    for cKey in fileKeys:
        if cKey in notDataKeys:
            notDataInDset.append(cKey)
        else:
            dataKey = cKey
            print('Data key could be: ' + dataKey)

    if dataKey is None:
        sys.exit("Data discovery failed! Ending run now...")

    return dataKey

def extract_intvl(intervalsToPlot, years, timePeriod, darr, handlesToPlot):
    ''' Extract intervals of interest '''
    for intvl in intervalsToPlot:
        startInd = np.where(years==intvl)[0][0]
        if intvl+timePeriod == 2100:
            endInd = len(years)
        else:
            endInd = np.where(years==intvl+timePeriod)[0][0]
        handlesToPlot.append(darr.isel(time=np.arange(startInd,endInd)))

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
    elif levOfInt == 'troposphere': #simple split on 200hPa for now
        if levName == 'plev':
            indTpause = find_closest_level(darr, 20000, levName=levName)
        else:
            indTpause = find_closest_level(darr, 200, levName='lev')
        levMask = levs > levs[indTpause]
        if np.ndim(darr) == 3:
            darr = darr[levMask,:,:]
            darr = darr.sum(dim=levName)
        else:
            darr = darr[:,levMask,:,:]
            darr = darr.sum(dim=levName)
    elif levOfInt == 'stratosphere':
        if levName == 'plev':
            indTpause = find_closest_level(darr, 20000, levName=levName)
        else:
            indTpause = find_closest_level(darr, 200, levName='lev')
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
    ''' Create pressure level strings for use in title and filenames. '''
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
        levStr = '' #e.g. 2m variables that have only one level

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

    elif isinstance(regionToPlot,dict): #region_library objects
        locStr = regionToPlot['regSaveStr']
        locTitleStr = regionToPlot['regSaveStr']

        lats = darr['lat'] #GLENS, SCIRIS, and SSP2-4.5 Control are all on the same grid to within 10^-6
        lons = darr['lon']
        if len(regionToPlot['regLons'])>2: #non-rectangular region that does not cross Prime Meridian
            gridMask = make_polygon_mask(lats, lons, regionToPlot['regLats'], regionToPlot['regLons'])
            darrMask = darr.copy()
            try:
                darrMask.data[:,~gridMask] = np.nan
            except:
                darrMask.data[:,:,~gridMask] = np.nan #Can't use .sel since gridMask is two-dimensional
        elif isinstance(regionToPlot['regLons'], tuple): #non-rectangular region that crosses Prime Meridian
            sGridMaskList = list()
            for sc in np.arange(0,len(regionToPlot['regLons'])):
                sGridMask = make_polygon_mask(lats, lons, regionToPlot['regLats'][sc], regionToPlot['regLons'][sc])
                sGridMaskList.append(sGridMask)
            gridMask = np.logical_or.reduce(sGridMaskList)
            darrMask = darr.copy()
            try:
                darrMask.data[:,~gridMask] = np.nan
            except:
                darrMask.data[:,:,~gridMask] = np.nan #Can't use .sel since gridMask is two-dimensional
        else: #rectangular region
            latMask = (lats>regionToPlot['regLats'][0]) & (lats<regionToPlot['regLats'][1])
            if regionToPlot['regLons'][0] < regionToPlot['regLons'][1]: #rectangle does not cross Prime Meridian
                lonMask = (lons>regionToPlot['regLons'][0]) & (lons<regionToPlot['regLons'][1])
            else: #rectangle crosses Prime Meridian
                lonMask = (lons>regionToPlot['regLons'][0]) | (lons<regionToPlot['regLons'][1])
            latsOfInt = lats[latMask]
            lonsOfInt = lons[lonMask]
            darrMask = darr.sel(lat=latsOfInt,lon=lonsOfInt)

        if areaAvgBool:
            latWeights = np.cos(np.deg2rad(darrMask['lat']))
            darrWght = darrMask.weighted(latWeights)
            darr = darrWght.mean(dim=['lat','lon'], skipna=True)
            # darr = darrMask.mean(dim=['lat','lon'], skipna=True)
        else:
            darr = darrMask

    elif isinstance(regionToPlot,list): #List of lat/lon (must be rectangular and not crossing the Prime Meridian)
        darr = darr.sel(lat=regionToPlot[0], lon=regionToPlot[1], method="nearest")

        latStr = str(np.round_(darr.lat.data,decimals=2))
        lonStr = str(np.round_(darr.lon.data,decimals=2))
        locStr = latStr + '_' + lonStr
        locTitleStr = '(' + latStr + ',' + lonStr + ')'

    else:
        sys.exit('Invalid region! Check value for regionToPlot.')

    return darr, locStr, locTitleStr

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
    ''' Obtain realization of interest or calculate ensemble mean, then create
    relevant filename. This, meta_book, and open_data have substantial
    hard-coding and require reworking when a new model run is added.
    '''
    try:
        if 'GLENS:Control' in darr.scenario:
            scnStr = 'gc' #glenscontrol
        elif 'GLENS:Feedback' in darr.scenario:
            scnStr = 'gf' #glensfeedback
        elif 'SCIRIS:Feedback' in darr.scenario:
            scnStr = 'sci' #sciris
        elif 'SCIRIS:Control' in darr.scenario:
            scnStr = 's245c' #ssp245control
        else:
            ic('Unknown scenario!')
            #No sys.exit(); want to know what the error is

        if setDict['realization'] == 'mean': #Output DataArray of ensemble mean
            darrMn = darr.mean(dim='realization')
            darrOut = darrMn.compute()
            ememSave = 'mn' + scnStr
        elif setDict['realization'] == 'ensplot': #Output DataArray of all members and ensemble mean
            darrMn = darr.mean(dim='realization')
            darrOut = xr.concat([darr,darrMn],dim='realization').compute() #Add ensemble mean as another "realization"
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

def open_data(dataDict, setDict):
    ''' Opens data and select data variable. This, manage_realizations, and
    meta_book have substantial hard-coding and require reworking when a new
    model run is added.
    '''
    if '*' in dataDict["idGlensCntrl"]:
        cntrlPath = dataDict["dataPath"] + dataDict["idGlensCntrl"]
        fdbckPath = dataDict["dataPath"] + dataDict["idGlensFdbck"]
        scirisPath = dataDict["dataPath"] + dataDict["idSciris"]
        s245CntrlPath = dataDict["dataPath"] + dataDict["idS245Cntrl"]

        glensCntrlDset = xr.open_mfdataset(cntrlPath, concat_dim='realization', combine='nested')
        glensFdbckDset = xr.open_mfdataset(fdbckPath, concat_dim='realization', combine='nested')
        try:
            scirisDset = xr.open_mfdataset(scirisPath, concat_dim='realization', combine='nested')
        except:
            ic('No SCIRIS files located')
            scirisDset = None
        try:
            s245CntrlDset = xr.open_mfdataset(s245CntrlPath, concat_dim='realization', combine='nested')
        except:
            ic('No SSP2-4.5 Control files located')
            s245CntrlDset = None

        if setDict['landmaskFlag'] == 'land':
            glensCntrlDset = glensCntrlDset.where(glensCntrlDset.landmask > 0)
            glensFdbckDset = glensFdbckDset.where(glensFdbckDset.landmask > 0)
            if scirisDset is not None:
                scirisDset = scirisDset.where(scirisDset.landmask > 0)
            if s245CntrlDset is not None:
                cesmMask = xr.open_dataset(dataDict["idCesmMask"])
                s245CntrlDset = s245CntrlDset.where(cesmMask.imask > 0)

        dataKey = discover_data_var(glensCntrlDset)
        glensCntrlDarr = glensCntrlDset[dataKey]
        glensCntrlDarr.attrs['scenario'] = 'GLENS:Control/RCP8.5'
        glensFdbckDarr = glensFdbckDset[dataKey]
        glensFdbckDarr.attrs['scenario'] = 'GLENS:Feedback/SAI/G1.2(8.5)'
        if scirisDset is not None:
            scirisDarr = scirisDset[dataKey]
            scirisDarr.attrs['scenario'] = 'SCIRIS:Feedback/SAI/G1.5(4.5)'
        else:
            scirisDarr = None
        if s245CntrlDset is not None:
            dataKeyCmip6 = discover_data_var(s245CntrlDset) #CMIP6 has unique variable names tied to CMIP6 conventions specifically, not its status as the control for SCIRIS
            s245CntrlDarr = s245CntrlDset[dataKeyCmip6]
            s245CntrlDarr.attrs['scenario'] = 'CMIP6_CESM2WACCM/SCIRIS:Control/SSP2-4.5'
        else:
            s245CntrlDarr = None
    else:
        sys.exit("Check input! Token should have a wildcard (i.e. match multiple files).")

    darrCheckList = list([glensCntrlDarr,glensFdbckDarr,scirisDarr,s245CntrlDarr])
    darrList = [d for d in darrCheckList if d is not None]

    return darrList, dataKey

def call_to_open(dataDict, setDict):
    ''' Common data tasks for all basic plots '''
    darrList, dataKey = open_data(dataDict, setDict)

    glensCntrlFiles = sorted(glob.glob(dataDict['dataPath'] + dataDict['idGlensCntrl']))
    glensFdbckFiles = sorted(glob.glob(dataDict['dataPath'] + dataDict['idGlensFdbck']))
    scirisFiles = sorted(glob.glob(dataDict['dataPath'] + dataDict['idSciris']))
    s245CntrlFiles = sorted(glob.glob(dataDict['dataPath'] + dataDict['idS245Cntrl']))

    ememList = list()
    for ec in (glensCntrlFiles, glensFdbckFiles, scirisFiles, s245CntrlFiles):
        ememList.append(get_ens_mem(ec))
    rlzList = list()
    ememStrList = list()
    for dc,darr in enumerate(darrList):
        rlzArr, ememStr = manage_realizations(setDict, darr, ememList[dc])
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
    ''' Compile bits and pieces for filenames and titles. As many as possible
    are derived automatically, but some are manual and require reworking when
    a new model run is added.
    '''
    metaDict = {
        "cntrlStr": 'RCP8.5',
        "fdbckStr": 'G1.2(8.5)',
        "scirisStr": 'G1.5(2-4.5)',
        "s245Cntrl": 'SSP2-4.5',
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
        "ensPid": {'spg': 'spghtti', 'sprd': 'spread'},
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

def average_over_years(darr, startYear, endYear):
    ''' Take average over a period of time '''
    datasetYears = darr['time'].dt.year.data
    startInd = int(np.where(datasetYears == startYear)[0])
    endInd = int(np.where(datasetYears == endYear)[0])
    intrvlOfInt = darr[startInd:endInd]
    darrMeanToi = intrvlOfInt.mean(dim='time')

    return darrMeanToi

def norm_by_absmax(darr):
    ''' Normalize data to within [-1,1] with respect to its max magnitude '''
    darrAbs = np.abs(darr)
    normValue = np.max(darrAbs)
    darrNorm = darr / normValue
    darrNorm.attrs['units'] = 'dimless'

    return darrNorm

def make_spc_string(setDict):
    ''' Describe area average Boolean or dimensions of variability '''
    if "dimOfVrblty" in setDict.keys():
        spcStr = 'rlz' + str(setDict["dimOfVrblty"]['rlzBool'])[0] + 'tm' + str(setDict["dimOfVrblty"]['timeBool'])[0] + 'spc' + str(setDict["dimOfVrblty"]['spcBool'])[0]
    elif "areaAvgBool" in setDict.keys():
        if setDict["areaAvgBool"]:
            spcStr = 'spcavg'
        else:
            spcStr = 'nospcavg'
    else:
        spcStr = ''

    return spcStr

def make_dov_title(dovDict):
    ''' Describe dimensions of variability for title '''
    dovStr = ''
    comma = ','
    for key in dovDict.keys():
        if dovDict[key]:
            dovStr = dovStr + key[:-4] + comma
    if dovStr[len(dovStr)-1] == ',':
        dovStr = dovStr[:-1]

    return dovStr

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
