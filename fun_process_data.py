''' fun_process_data
Contains functions for operations like extracting metadata from filenames,
choosing and summing levels, or managing area inputs. These are designed to be
as flexible as possible, but manage_realizations, meta_book, and call_to_open
have unavoidable hard-coding and require reworking to add new model runs.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import collections as col
import dask
dask.config.set(**{'array.slicing.split_large_chunks': True})
from datetime import date
import glob
import numpy as np
import matplotlib.path as mpth
import xarray as xr
xr.set_options(keep_attrs=True)

import CustomExceptions
import fun_convert_unit as fcu

def call_to_open(dataDict, setDict):
    ''' Open data and carry out tasks common to all types '''
    lf = col.defaultdict(list) # List factory stackoverflow.com/a/2402678
    # Open datasets for each input experiment
    for dky in dataDict.keys():
        if 'id' in dky:
            try:
                inPath = dataDict["dataPath"] + dataDict[dky]
                inGlobs = sorted(glob.glob(inPath))
                lf['globs'].append(inGlobs) # Used to track realizations
                rawDset = xr.open_mfdataset(
                    inPath, concat_dim='realization', 
                    combine='nested', coords='minimal')
                try:
                    sourceId = rawDset.attrs['source_id']
                except:
                    sourceId = None
                dataKey = discover_data_var(rawDset)
                rawDarr = rawDset[dataKey]
                scnDarr = bind_scenario(rawDarr, dataDict[dky])
                if scnDarr.scenario == 'PreindustrialControl': 
                # TODO: make this elegant for all
                    scnDarr.attrs['scenario'] = sourceId + ':' + scnDarr.scenario
                maskDarr = apply_mask(scnDarr, dataDict, setDict)
                if 'CESM2-ARISE:Control' in scnDarr.scenario: # Two parts: historical and future
                    lf['chf'].append(maskDarr) # Keep separate for now
                else:
                    lf['darr'].append(maskDarr) # Put other data in darrList
            except Exception as fileOpenErr: # Reached for None input
                ic(fileOpenErr) # Display reason for error
                pass # Move on if possible

    # Handle the two parts of the CESM2-ARISE Control                
    if len(lf['chf']) == 2: # Hist & future input
        acntrlDarr = combine_hist_fut(lf['chf'][0], lf['chf'][1])
        lf['darr'].append(acntrlDarr)
        acntrlDarr.attrs['scenario'] = 'CESM2-WACCM/ARISE:Control/SSP2-4.5'
    else: # Either hist OR future input
        try:
            lf['darr'].append(lf['chf'][0]) # Use whichever is present
        except:
            pass # If there is no data, move on if possible
            
    if len(lf['darr']) == 0:
        raise CustomExceptions.NoDataError(
            'No data! Check input and try again.')

    # Manage ensemble members
    for ec in lf['globs']:
        lf['emem'].append(get_ens_mem(ec))
    lf['emem'] = list(filter(None, lf['emem']))
    for dc,darr in enumerate(lf['darr']):
        scnArr, ememStr = manage_realizations(
            setDict, darr, lf['emem'][dc])
        lf['scn'].append(scnArr)
        lf['ememStr'].append(ememStr)
    lf['ememStr'] = list(filter(None,lf['ememStr']))
    ememSave = '-'.join(lf['ememStr'])
    cmnDict = {'dataKey': dataKey, 'ememSave': ememSave}
    
    # Convert units and calculate variables
    if setDict["convert"] is not None:
        for fcufcv in setDict["convert"]:
            for rc,rv in enumerate(lf['scn']):
                try:
                    lf['scn'][rc] = fcufcv(rv) # Apply converter or calculator function(s)
                except:
                    lf['scn'][rc] = fcufcv(rv, setDict) # Sometimes they need setDict
    
    return lf['scn'], cmnDict

def apply_mask(darr, dataDict, setDict):
    ''' Apply land/ocean mask, returning masked dataset '''
    if setDict['landmaskFlag'] is not None:
        if "UKESM" in darr.scenario:
            activeMaskDset = xr.open_dataset(dataDict["maskUkesm"])
            activeMask = activeMaskDset.mask
        elif "CESM" in darr.scenario:
            activeMaskDset = xr.open_dataset(dataDict["mask"])
            try:
                activeMask = activeMaskDset.landmask
            except:
                activeMask = activeMaskDset.imask

        if setDict['landmaskFlag'] == 'land':
            maskDarr = darr.where(activeMask > 0)
        elif setDict['landmaskFlag'] == 'ocean':
            maskDarr = darr.where(activeMask == 0)
        else:
            ic('Invalid landmaskFlag! Continuing with no mask applied')
            maskDarr = darr
    else: # When no mask is applied
        maskDarr = darr

    return maskDarr

def discover_data_var(dset):
    ''' Find data var among the many keys in an ESM file '''
    fileKeys = list(dset.keys())
    notDataKeys = [ # Manually add keys that don't represent data here
        'time_bnds', 'date', 'datesec', 'lev_bnds', 'gw', 'ch4vmr',
        'co2vmr', 'ndcur', 'nscur', 'sol_tsi', 'nsteph', 'f11vmr',
        'n2ovmr', 'f12vmr', 'lon_bnds', 'lat_bnds', 'ZSOI', 'BSW',
        'WATSAT', 'landmask', 'ZLAKE', 'DZLAKE', 'SUCSAT', 'area',
        'landfrac', 'topo', 'DZSOI', 'pftmask', 'HKSAT', 'nstep',
        'mdcur', 'mscur', 'mcdate', 'mcsec', 'nbedrock',
        'binary_mhw_start']
    notDataInDset = list()
    dataKey = None

    for cKey in fileKeys: # Check all keys in file
        if cKey in notDataKeys:
            notDataInDset.append(cKey) # Useful for troubleshooting
        else:
            dataKey = cKey
            print('Data key could be: ' + dataKey)

    if dataKey is None:
        sys.exit("Data discovery failed! Ending run now...")

    return dataKey

def bind_scenario(darr, inID):
    ''' Bind scenario identifiers to xarray attributes '''
    if 'control' in inID:
        darr.attrs['scenario'] = 'CESM1-GLENS:Control/RCP8.5/No-SAI/RCP8.5'
    elif 'feedback' in inID:
        darr.attrs['scenario'] = 'CESM1-GLENS:Feedback/SAI/GLENS-SAI'
    elif 'DEFAULT' in inID:
        darr.attrs['scenario'] = 'CESM2-ARISE:Feedback/SAI/ARISE-SAI-1.5'
    elif 'DELAYED' in inID:
        darr.attrs['scenario'] = 'CESM2-ARISE-DelayedStart:Feedback/SAI/ARISE-SAI-1.5-DelayedStart'
    elif 'ARISE1P0' in inID:
        darr.attrs['scenario'] = 'CESM2-ARISE-1.0:Feedback/SAI/ARISE-SAI-1.0'
    elif 'BWSSP245' in inID:
        darr.attrs['scenario'] = 'CESM2-WACCM/CESM2-ARISE:Control/No-SAI/SSP2-4.5'
    elif 'BWHIST' in inID:
        darr.attrs['scenario'] = 'CESM2-WACCM/CESM2-ARISE:Control/Historical/No-SAI'
    elif 'SSP245cmip6' in inID:
        darr.attrs['scenario'] = 'CESM2-WACCM/CESM2-ARISE:Control/SSP2-4.5/CMIP6/No-SAI'
    elif 'ssp245' in inID:
        darr.attrs['scenario'] = 'UKESM-ARISE:Control/No-SAI/UKESM-SSP2-4.5'
    elif 'arise-sai-1p5' in inID:
        darr.attrs['scenario'] = 'UKESM-ARISE:Feedback/SAI/UKESM-ARISE-SAI-1.5'
    elif 'piControl' in inID:        
        darr.attrs['scenario'] = 'PreindustrialControl'
    else:
        ic('Unable to match scenario, binding empty string to array')
        darr.attrs['scenario'] = ''

    return darr

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
            activeDarr = xr.concat(
                (activeHist, activeCntrl), dim='time')
        except:
            activeCntrl = darrCntrl.sel(realization=rc)
            activeDarr = xr.concat(
                (darrHistNan,activeCntrl), dim='time')
        darrCombine = xr.concat(
            (darrCombine, activeDarr), dim='realization')
    darrCombine = darrCombine.sel( # Ensure correct num of realizations
        realization=np.arange(0, len(cntrlRlz)))

    return darrCombine

def copy_blank_darr(darr):
    ''' Make a copy of a DataArray with blank data but the same structure '''
    darrShape = np.shape(darr.data)
    darrNan = np.full(darrShape, np.nan)
    darrOut = darr.copy(data=darrNan)

    return darrOut

def manage_realizations(setDict, darr, emem):
    ''' Obtain single realization or ensemble mean and make useful filename '''
    try:
        if 'GLENS:Control' in darr.scenario:
            scnStr = 'gc' #glenscontrol
        elif 'GLENS:Feedback' in darr.scenario:
            scnStr = 'gf' #glensfeedback
        elif 'ARISE:Feedback' in darr.scenario:
            scnStr = 'arif' #arisefeedback
        elif 'ARISE-DelayedStart:Feedback' in darr.scenario:
            scnStr = 'aridsf' #arisedelayedstartfeedback
        elif 'ARISE-1.0:Feedback' in darr.scenario:
            scnStr = 'ari1p0f' #arise1p0feedback
        elif 'ARISE:Control' in darr.scenario:
            scnStr = 'aric' #arisecontrol
        elif 'UKESM-ARISE:No-SAI' in darr.scenario:
            scnStr = 'ukan' #ukarisenosai
        elif 'UKESM-ARISE:SAI' in darr.scenario:
            scnStr = 'ukas' #ukarisesai
        elif 'PreindustrialControl' in darr.scenario:
            scnStr = 'pi'
        else:
            ic('Unknown scenario!')
            #No sys.exit(); want to know what the error is

        if setDict['realization'] == 'mean': #Output DataArray of ens mean
            darrMn = darr.mean(dim='realization')
            darrOut = darrMn.compute()
            ememSave = 'mn' + scnStr
        elif setDict['realization'] == 'ensplot': #Output DataArray of members AND ens mean
            if setDict["areaAvgBool"] == 'sum':
                darrOut = darr.compute()
            else:
                darrMn = darr.mean(dim='realization')
                darrOut = xr.concat([darr,darrMn], dim='realization').compute() #Add ens mean as another "realization"
            ememSave = 'ens' + scnStr
        elif setDict['realization'] == 'allplot':
            #TODO: this should be only behavior, and could maybe
            # even result in removing this function completely.
            # Taking ensemble statistics should be the very LAST step
            # before plotting, rather than coming as part of the basic
            # data opening tasks in call_to_open.
            # This change is easy for the slice globes which already
            # don't use any of this, but it will break the basicplots
            # and ensplots functions and require considerable 
            # refactoring to repair. This is not the highest priority
            # right now but will be dealt with eventually!
            darrOut = darr.compute()
            ememSave = 'ens' + scnStr
        else: #Output DataArray of single ens member
            try:
                ememNum = list(map(int, emem))
            except:
                ememNum  = list(map(str, emem))
            rInd = ememNum.index(setDict['realization'])
            darrOut = darr.isel(realization=rInd).compute()
            activeEmem = emem[rInd]
            ememSave = scnStr + activeEmem
    except:
        darrOut = None
        ememSave = ''

    return darrOut, ememSave
    
def get_ens_mem(files):
    ''' Find ensemble members from a list of files '''
    emem = list()
    for pth in files:
        pthPcs = pth.split('/')
        fName = pthPcs[len(pthPcs)-1]
        fNamePcs = fName.split('_')
        # ic(fName, fNamePcs)
        emem.append(fNamePcs[1])
        
    return emem

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
        if areaAvgBool == True:
            latWeights = np.cos(np.deg2rad(darr['lat']))
            darrWght = darr.weighted(latWeights)
            darr = darrWght.mean(dim=['lat','lon'], skipna=True)
        elif areaAvgBool == 'sum':
            latWeights = np.cos(np.deg2rad(darr['lat']))
            darrWght = darr.weighted(latWeights)
            darr = darrWght.sum(dim=['lat','lon'], skipna=True)
    elif isinstance(regionToPlot,dict): #region_library objects
        if len(regionToPlot["regLats"]) == 1: # Point region_library object
            if 'lat' in darr.dims: #If selecting data as opposed to making strings
                darr = darr.sel(
                    lat=regionToPlot['regLats'],
                    lon=regionToPlot['regLons'],
                    method="nearest")
                ic(darr['lat'].data, darr['lon'].data) #Display true lat/lon as 'nearest' is used
                darr = np.squeeze(darr) #Drop length 1 lat/lon dimensions
                areaAvgBool = None
            locStr = regionToPlot["regSaveStr"]
            locTitleStr = regionToPlot["regStr"]
            return darr, locStr, locTitleStr # Shortcut the rest of the function
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

        if areaAvgBool == True:
            latWeights = np.cos(np.deg2rad(darrMask['lat']))
            darrWght = darrMask.weighted(latWeights)
            darr = darrWght.mean(dim=['lat','lon'], skipna=True)
            # darr = darrMask.mean(dim=['lat','lon'], skipna=True) #No latitude weighting
        elif areaAvgBool == 'sum':
            darr = darrMask.sum(dim=['lat','lon'], skipna=True)
        elif areaAvgBool == False:
            darr = darrMask

    elif isinstance(regionToPlot,list): #List of lat/lon (must be rectangular and not crossing the Prime Meridian)
        darr = darr.sel(lat=regionToPlot[0], lon=regionToPlot[1], method="nearest")

        latStr = str(np.round_(darr.lat.data,decimals=2))
        lonStr = str(np.round_(darr.lon.data,decimals=2))
        locStr = latStr + '_' + lonStr
        locTitleStr = '(' + latStr + ',' + lonStr + ')'

    else:
        sys.exit('Invalid region! Check value for regionToPlot.')

    if areaAvgBool == 'sum': #NaNs freak out summations because it results in a 0
        try: #spaghetti
            cursedZeros = (darr == 0)
            darr[cursedZeros] = np.nan
        except: #spread
            for rc in np.arange(0,len(darr['realization'])):
                cursedZeros = (darr.isel(realization=rc) == 0)
                darr.isel(realization=rc)[cursedZeros] = np.nan

    return darr, locStr, locTitleStr

def meta_book(setDict, dataDict, darr):
    ''' Compile bits and pieces of metadata for filenames and titles '''

    metaDict = {
        "cntrlStr": 'RCP8.5',
        "fdbckStr": 'GLENS-SAI',
        "ariseStr": 'ARISE-SAI-1.5',
        "s245Cntrl": 'SSP2-4.5',
        "varStr": var_str_lookup(darr.long_name, setDict, strType='title'),
        "varSve": var_str_lookup(darr.long_name, setDict, strType='save'),
        "strtStr": str(darr['time'].data[0].year) if 'time' in darr.dims else '',
        "endStr": str(darr['time'].data[len(darr)-1].year)  if 'time' in darr.dims else '',
        "frstDcd": {
            "GLENS": '',
            "CESM2-ARISE": '',
            "UKESM-ARISE": ''
            },
        "lstDcd": {
            "GLENS": '',
            "CESM2-ARISE": '',
            "UKESM-ARISE": ''
            },
        "levStr": make_level_string(darr, setDict["levOfInt"]) if "levOfInt" in setDict.keys() else '',
        "levSve": make_level_string(darr, setDict["levOfInt"]).replace(" ","") if "levOfInt" in setDict.keys() else '',
        "ensStr": dataDict["ememSave"],
        "yStr": darr.units,
        "unit": darr.attrs['units'],
        "spcStr": make_spc_string(setDict),
        "pid": {'g1p': 'globe_1p', 'g2p': 'globe_2p', 'g4p': 'globe_4p', 'g6p': 'globe_6p', 'ts': 'timeseries'},
        "ensPid": {'spg': 'spghtti', 'sprd': 'spread'},
        "glbType": {'vGl': 'vertical', 'bGl': 'baseline', 'fcStr': 'FdbckCntrl', 'gfcStr': 'glensFdbckCntrl'}
    }

    # Fill out start and end interval dictionaries
    if "strtIntvl" in setDict.keys():
        if setDict["strtIntvl"]["GLENS"] is not None:
            metaDict["frstDcd"]["GLENS"] = str(setDict["strtIntvl"]["GLENS"][0]) \
                + '-' + str(setDict["strtIntvl"]["GLENS"][1]-1)
        else:
            metaDict["frstDcd"]["GLENS"] = ''
        if setDict["strtIntvl"]["CESM2-ARISE"] is not None:
            metaDict["frstDcd"]["CESM2-ARISE"] = str(setDict["strtIntvl"]["CESM2-ARISE"][0]) \
                + '-' + str(setDict["strtIntvl"]["CESM2-ARISE"][1]-1)
        else:
            metaDict["frstDcd"]["CESM2-ARISE"] = ''
        if setDict["strtIntvl"]["UKESM-ARISE"] is not None:
            metaDict["frstDcd"]["UKESM-ARISE"] = str(setDict["strtIntvl"]["UKESM-ARISE"][0]) \
                + '-' + str(setDict["strtIntvl"]["UKESM-ARISE"][1]-1)
        else:
            metaDict["frstDcd"]["UKESM-ARISE"] = ''

    if "endIntvl" in setDict.keys():
        if setDict["endIntvl"]["GLENS"] is not None:
            metaDict["lstDcd"]["GLENS"] = str(setDict["endIntvl"]["GLENS"][0]) \
                + '-' + str(setDict["endIntvl"]["GLENS"][1]-1)
        else:
            metaDict["lstDcd"]["GLENS"] = ''
        if setDict["endIntvl"]["CESM2-ARISE"] is not None:
            metaDict["lstDcd"]["CESM2-ARISE"] = str(setDict["endIntvl"]["CESM2-ARISE"][0]) \
                + '-' + str(setDict["endIntvl"]["CESM2-ARISE"][1]-1)
        else:
            metaDict["lstDcd"]["CESM2-ARISE"] = ''
        if setDict["endIntvl"]["UKESM-ARISE"] is not None:
            metaDict["lstDcd"]["UKESM-ARISE"] = str(setDict["endIntvl"]["UKESM-ARISE"][0]) \
                + '-' + str(setDict["endIntvl"]["UKESM-ARISE"][1]-1)
        else:
            metaDict["lstDcd"]["UKESM-ARISE"] = ''

    return metaDict

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
    elif longName == 'Near-Surface Air Temperature':
        if strType == 'title':
            outStr = '2m air temperature'
        elif strType == 'save':
            outStr = '2mtemp'
    elif 'Decadal climate distance of 2m temperature' in longName:
        if strType == 'title':
            outStr = longName #Default is fine here
        elif strType == 'save':
            yrs = longName[-9:].replace('-','')
            outStr = 'dcd-2mtemp-' + yrs
    elif 'Climate speed of 2m temperature' in longName:
        if strType == 'title':
            outStr = longName #Default is fine here
        elif strType == 'save':
            yrs = longName[-9:].replace('-','')
            outStr = 'cspd-2mtemp-' + yrs
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

def period_month_avg(darrList):
    darrPerAvgList = list()
    for darr in darrList:
        darrPerAvg = darr.groupby("time.year").mean()
        darrPerAvg = darrPerAvg.rename({'year':'time'})
        darrPerAvgList.append(darrPerAvg)

    return darrPerAvgList

def make_ord_array(inTimes):
    ''' Make ordinal time array required by e.g. mhws package '''
    tCftime = xr.cftime_range(inTimes[0], inTimes[len(inTimes)-1])
    checkMissing = ic(list(set(inTimes).symmetric_difference(set(tCftime)))) #Print missing timesteps between input and cfrange
    ordList = list()
    for tcf in tCftime:
        dtAc = tcf
        dtAcOrd = date(dtAc.year, dtAc.month, dtAc.day).toordinal()
        ordList.append(dtAcOrd)
    ordArr = np.array(ordList)

    return ordArr
