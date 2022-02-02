''' fun_process_data
Contains functions for operations like extracting metadata from filenames,
choosing and summing levels, or managing area inputs.

Written by Daniel Hueholt
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

import CustomExceptions
import fun_convert_unit as fcu

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
    originalPoi = timePeriod
    # ic(originalPoi)
    for intvl in intervalsToPlot:
        # ic(intvl)
        if intvl == 2015:
            timePeriod = 15
        else:
            timePeriod = originalPoi
        # ic(timePeriod)
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

        lats = darr['lat'] #GLENS, ARISE, and SSP2-4.5 Control are all on the same grid to within 10^-6
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
        elif 'ARISE:Feedback' in darr.scenario:
            scnStr = 'ari' #arise
        elif 'ARISE:Control' in darr.scenario:
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

def apply_mask(dset, dataDict, setDict):
    if setDict['landmaskFlag'] is not None:
        activeMaskDset = xr.open_dataset(dataDict["mask"])
        try:
            activeMask = activeMaskDset.landmask
        except:
            activeMask = activeMaskDset.imask

        if setDict['landmaskFlag'] == 'land':
            maskDset = dset.where(activeMask > 0)
        elif setDict['landmaskFlag'] == 'ocean':
            maskDset = dset.where(activeMask == 0)
        else:
            ic('Invalid landmaskFlag! Continuing with no mask applied')
            maskDset = dset
    else:
        ic('landmaskFlag is None, no mask applied')
        maskDset = dset

    return maskDset

def bind_scenario(darr, inID):
    if 'control_' in inID:
        darr.attrs['scenario'] = 'GLENS:Control/RCP8.5'
    elif 'feedback_' in inID:
        darr.attrs['scenario'] = 'GLENS:Feedback/SAI/G1.2(8.5)'
    elif 'SSP245-TSMLT-GAUSS' in inID:
        darr.attrs['scenario'] = 'ARISE:Feedback/SAI/G1.5(4.5)'
    elif 'BWSSP245' in inID:
        darr.attrs['scenario'] = 'CESM2-WACCM/ARISE:Control/SSP2-4.5'
    elif 'BWHIST' in inID:
        darr.attrs['scenario'] = 'CESM2-WACCM/ARISE:Control/Historical'
    elif 'SSP245cmip6' in inID:
        darr.attrs['scenario'] = 'CESM2-WACCM/ARISE:Control/SSP2-4.5/ProcessedForCMIP6'
    else:
        ic('Unable to match scenario, binding empty string to array')
        darr.attrs['scenario'] = ''

    return darr

def call_to_open(dataDict, setDict):
    ''' Common data tasks for all basic plots '''
    # Open datasets
    darrList = list()
    cmbnHistFutList = list()
    globsList = list()
    for dky in dataDict.keys():
        if 'id' in dky:
            try:
                inPath = dataDict["dataPath"] + dataDict[dky]
                globsList.append(sorted(glob.glob(inPath)))
                rawDset = xr.open_mfdataset(inPath, concat_dim='realization', combine='nested')
                maskDset = apply_mask(rawDset, dataDict, setDict)
                dataKey = discover_data_var(maskDset)
                maskDarr = maskDset[dataKey]
                scnDarr = bind_scenario(maskDarr, dataDict[dky])
                if 'ARISE:Control' in scnDarr.scenario: #Two parts to ARISE Control: historical and future
                    cmbnHistFutList.append(scnDarr) #These need to be kept separate
                else:
                    darrList.append(scnDarr) #Others can be put in darrList directly
            except: #Usually reached if input is None
                pass
    if len(cmbnHistFutList) == 2: #If both historical and future are input
        acntrlDarr = combine_hist_fut(cmbnHistFutList[0],cmbnHistFutList[1]) #Combine ARISE Control here
        darrList.append(acntrlDarr) #Append ARISE Control to darrList
        acntrlDarr.attrs['scenario'] = 'CESM2-WACCM/ARISE:Control/SSP2-4.5' #THIS IS CLUMSY, TODO: # FIX
    else:
        try:
            darrList.append(cmbnHistFutList[0]) #Append the one that's present
        except:
            pass #If there is no data, just move on
    if len(darrList) == 0:
        raise CustomExceptions.NoDataError('No data! Check input and try again.')

    # Manage ensemble members
    ememList = list()
    for ec in globsList:
        ememList.append(get_ens_mem(ec))
    scnList = list()
    ememStrList = list()
    for dc,darr in enumerate(darrList):
        scnArr, ememStr = manage_realizations(setDict, darr, ememList[dc])
        scnList.append(scnArr)
        ememStrList.append(ememStr)
    ememStrList = list(filter(None,ememStrList))
    ememSave = '-'.join(ememStrList)
    cmnDict = {'dataKey': dataKey, 'ememSave': ememSave}

    # Convert units (if necessary)
    if setDict["convert"] is not None:
        for cnvrtr in setDict["convert"]:
            for rc,rv in enumerate(scnList):
                scnList[rc] = cnvrtr(rv) #Use input converter function(s)

    return scnList, cmnDict

def meta_book(setDict, dataDict, cntrlToPlot, labelsToPlot=None):
    ''' Compile bits and pieces for filenames and titles. As many as possible
    are derived automatically, but some are manual and require reworking when
    a new model run is added.
    '''
    metaDict = {
        "cntrlStr": 'RCP8.5',
        "fdbckStr": 'G1.2(8.5)',
        "ariseStr": 'G1.5(2-4.5)',
        "s245Cntrl": 'SSP2-4.5',
        "varStr": var_str_lookup(cntrlToPlot.long_name, setDict, strType='title'),
        "varSve": var_str_lookup(cntrlToPlot.long_name, setDict, strType='save'),
        "strtStr": str(cntrlToPlot['time'].data[0].year),
        "endStr": str(cntrlToPlot['time'].data[len(cntrlToPlot)-1].year),
        "frstDcd": str(setDict["startIntvl"][0]) + '-' + str(setDict["startIntvl"][1]-1),
        "aFrstDcd": str(setDict["startIntvl"][2]) + '-' + str(setDict["startIntvl"][3]-1) if len(setDict["startIntvl"]) > 2 else '',
        "lstDcd": str(setDict["endIntvl"][0]) + '-' + str(setDict["endIntvl"][1]-1),
        "aLstDcd": str(setDict["endIntvl"][2]) + '-' + str(setDict["endIntvl"][3]-1) if len(setDict["startIntvl"]) > 2 else '',
        "tmStr": rcf_parser(labelsToPlot),
        "levStr": make_level_string(cntrlToPlot, setDict["levOfInt"]) if "levOfInt" in setDict.keys() else '',
        "levSve": make_level_string(cntrlToPlot, setDict["levOfInt"]).replace(" ","") if "levOfInt" in setDict.keys() else '',
        "ensStr": dataDict["ememSave"],
        "yStr": cntrlToPlot.units,
        "unit": cntrlToPlot.attrs['units'],
        "pdfStyle": setDict["plotStyle"],
        "spcStr": make_spc_string(setDict),
        "pid": {'g1p': 'globe_1p', 'g2p': 'globe_2p', 'g4p': 'globe_4p', 'g6p': 'globe_6p', 'ts': 'timeseries', 'pdf': 'pdf'},
        "ensPid": {'spg': 'spghtti', 'sprd': 'spread'},
        "glbType": {'vGl': 'vertical', 'bGl': 'baseline', 'fcStr': 'FdbckCntrl', 'gfcStr': 'glensFdbckCntrl'}
    }

    return metaDict

def rcf_parser(labelsToPlot):
    ''' Parse labels to make filenames '''
    timeStr = list()
    if labelsToPlot != None:
        for lab in labelsToPlot:
            if 'Reference' in lab:
                timeStr.append('r' + lab[0:9].replace("-","")) #Reference
            elif 'RCP8.5' in lab:
                timeStr.append('gc' + lab[0:9].replace("-","")) #GLENS control
            elif 'SSP2-4.5' in lab:
                timeStr.append('ac' + lab[0:9].replace("-","")) #ARISE control
            elif 'G1.2(8.5)' in lab:
                timeStr.append('gf' + lab[0:9].replace("-","")) #GLENS feedback
            elif 'G1.5(2-4.5)' in lab:
                timeStr.append('af' + lab[0:9].replace("-","")) #ARISE feedback
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
    ic(startYear, endYear)
    startInd = int(np.where(datasetYears == startYear)[0])
    endInd = int(np.where(datasetYears == endYear)[0])
    intrvlOfInt = darr[startInd:endInd]
    ic(intrvlOfInt.time)
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

def combine_hist_fut(darrCntrl, darrHist):
    ''' Combine historical and future output into a single DataArray '''
    darrHistForFormat = darrHist.sel(realization=0)
    darrHistNan = copy_blank_darr(darrHistForFormat)
    darrCombine = darrHistNan.copy()

    cntrlRlz = darrCntrl['realization']
    for rc in cntrlRlz:
        try:
            activeHist = darrHist.sel(realization=rc)
            activeCntrl = darrCntrl.sel(realization=rc)
            activeDarr = xr.concat((activeHist,activeCntrl), dim='time')
        except:
            activeCntrl = darrCntrl.sel(realization=rc)
            activeDarr = xr.concat((darrHistNan,activeCntrl), dim='time')
        darrCombine = xr.concat((darrCombine,activeDarr), dim='realization')
    darrCombine = darrCombine.sel(realization=np.arange(1,len(cntrlRlz)))

    return darrCombine

def copy_blank_darr(darr):
    ''' Make a copy of a DataArray with blank data but the same structure '''
    darrShape = np.shape(darr.data)
    darrNan = np.full(darrShape, np.nan)
    darrOut = darr.copy(data=darrNan)

    return darrOut

def var_str_lookup(longName, setDict, strType='title'):
    ''' For known variables, give a better name than the default '''
    if longName == 'Reference height temperature':
        if strType == 'title':
            outStr = '2m air temperature'
        elif strType == 'save':
            outStr = '2mtemp'
        else:
            outStr = None
    elif ((longName == 'Total precipitation (liq+ice)') or (longName == 'Total (convective and large-scale) precipitation  (liq + ice)')):
        if strType == 'title':
            outStr = 'Total precipitation (liq+ice)'
        elif strType == 'save':
            outStr = 'TotalPrecip(liq+ice)'
    elif longName == 'Ice Fraction from Coupler':
        if strType == 'title':
            outStr = 'Ice Fraction from Coupler'
        elif strType == 'save':
            outStr = 'IceFracCoupler'
    elif longName == 'Annual tropical nights':
        if strType == 'title':
            outStr = 'Annual tropical nights'
        elif strType == 'save':
            outStr = 'clxTR'
    else:
        outStr = longName

    if setDict["landmaskFlag"] == 'land':
        if strType == 'title':
            outStr = outStr + ' ' + 'land'
        elif strType == 'save':
            outStr = outStr + 'land'
        else:
            outStr = None
    elif setDict["landmaskFlag"] == 'ocean':
        if strType == 'title':
            outStr = outStr + ' ' + 'ocean'
        elif strType == 'save':
            outStr = outStr + 'ocean'
        else:
            outStr = None

    return outStr

def check_last_time(inDlyDarr):
    ''' Removes bonus timestep from ARISE daily data when it crosses year '''
    inTimes = inDlyDarr.time
    inYears = inTimes.dt.year
    bonusInd = len(inTimes)-1
    bonusTimeYr = inYears[bonusInd]
    lastNominalTimeYr = inYears[bonusInd-1]
    if bonusTimeYr == lastNominalTimeYr:
        outDlyDarr = inDlyDarr.copy()
    else:
        outDlyDarr = inDlyDarr.where(inYears<bonusTimeYr, drop=True)

    return outDlyDarr

def convert_for_consistency(inDarr):
    ''' CMIP6 has different unit standards, but the CESM2-WACCM CMIP6 run is
        used as the ARISE control. It's therefore necessary to convert this
        output so it's consistent. These conversions are not automatic and are
        purely determined by trial and error. '''
    if inDarr.standard_name == 'precipitation_flux':
        consistentDarr = fcu.flux_to_prect(inDarr)
    elif inDarr.standard_name == 'sea_ice_area_fraction':
        consistentDarr = fcu.perc_to_frac(inDarr)
    else:
        ic('Unknown CMIP6 variable! No unit conversion applied.')
        consistentDarr = inDarr

    return consistentDarr

def period_month_avg(darrList):
    darrPerAvgList = list()
    for darr in darrList:
        darrPerAvg = darr.groupby("time.year").mean()
        darrPerAvg = darrPerAvg.rename({'year':'time'})
        darrPerAvgList.append(darrPerAvg)

    return darrPerAvgList
