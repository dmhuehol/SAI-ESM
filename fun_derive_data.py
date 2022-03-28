''' fun_derive_data
Functions to derive and save data from a base dataset as a new file.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

from cftime import DatetimeNoLeap as dtnl
import cftime
from datetime import date
import numpy as np
import xarray as xr
import climdex.temperature as tdex
import climdex.precipitation as pdex
import marineHeatWaves as mhws

import fun_convert_unit as fcu
import fun_process_data as fpd
import region_library as rlib

### Climdex extremes
def derive_annual_tropical_nights(inFileTrefht, outPath):
    ''' Obtain annual tropical nights from daily min temperature data '''
    tIndices = tdex.indices(time_dim='time')
    inKey = 'TREFHTMN' #TREFHTMN, TSMN
    outKey = 'clxTR'

    trefhtDset = xr.open_dataset(inFileTrefht)
    trefhtDarr = trefhtDset[inKey]
    trefhtCelDarr = fcu.kel_to_cel(trefhtDarr)
    ic('Beginning clxTR calculation!')
    trefhtCelDarrChkLstTm = fpd.check_last_time(trefhtCelDarr)
    atnDarr = tIndices.annual_tropical_nights(trefhtCelDarrChkLstTm)
    ic('clxTR calculation complete!')

    timeList = list()
    for yr in atnDarr.year:
        timeList.append(dtnl(yr,7,15,12,0,0,0)) #Year with standard fill values
    newDset = xr.Dataset(
        {outKey: (("time","lat","lon"), atnDarr.data)},
        coords={
            "time": (('time'), timeList),
            "lat": atnDarr['lat'],
            "lon": atnDarr['lon']
        }
    )
    newDset[outKey].attrs = trefhtCelDarr.attrs
    newDset[outKey].attrs['long_name'] = 'Annual tropical nights'
    newDset[outKey].attrs['units'] = 'd/yr'

    inPcs = inFileTrefht.split('/') #inFileTrefht is the entire path to file
    inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
    strOut = inFn.replace(inKey, outKey) #Replace var name with extreme key
    outFile = outPath + strOut
    newDset.to_netcdf(outFile) #Save data
    ic(strOut, outFile, newDset) #icecream all useful parts of the output

def derive_annual_count_precip_gte20(inFilePrect, outPath):
    ''' Count annual days with precip >= 20mm from daily precip data '''
    pIndices = pdex.indices(time_dim='time')
    inKey = 'PRECT'
    outKey = 'r20mm'

    prectDset = xr.open_dataset(inFilePrect)
    prectDarr = prectDset[inKey]
    prectMmDarr = fcu.m_to_mm(prectDarr)
    prectMmDayDarr = fcu.persec_perday(prectMmDarr)
    ic('Beginning r20mm calculation!')
    prectMmDayDarrChkLstTm = fpd.check_last_time(prectMmDayDarr)
    r20mmDarr = pIndices.annual_r20mm(prectMmDayDarrChkLstTm)
    ic('r20mm calculation complete!')

    timeList = list()
    for yr in r20mmDarr.year:
        timeList.append(dtnl(yr,7,15,12,0,0,0)) #Year with standard fill values
    newDset = xr.Dataset(
        {outKey: (("time","lat","lon"), r20mmDarr.data)},
        coords={
            "time": (('time'), timeList),
            "lat": r20mmDarr['lat'],
            "lon": r20mmDarr['lon']
        }
    )
    newDset[outKey].attrs = prectMmDayDarr.attrs
    newDset[outKey].attrs['long_name'] = 'Annual count of days with precip gte 20mm'
    newDset[outKey].attrs['units'] = 'd/yr'

    inPcs = inFilePrect.split('/') #inFilePrect is the entire path to file
    inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
    strOut = inFn.replace(inKey, outKey) #Replace var name with extreme key
    outFile = outPath + strOut
    newDset.to_netcdf(outFile) #Save data
    ic(strOut, outFile, newDset) #icecream all useful parts of the output

def derive_simple_intensity_index(inFilePrect, outPath):
    ''' Derive simple intensity index from daily precipitation data '''
    pIndices = pdex.indices(time_dim='time')
    inKey = 'PRECT'
    outKey = 'sdii'

    prectDset = xr.open_dataset(inFilePrect)
    prectDarr = prectDset[inKey]
    prectMmDarr = fcu.m_to_mm(prectDarr)
    prectMmDayDarr = fcu.persec_perday(prectMmDarr)
    ic('Beginning sdii calculation!')
    prectMmDayDarrChkLstTm = fpd.check_last_time(prectMmDayDarr)
    sdiiDarr = pIndices.sdii(prectMmDayDarrChkLstTm, period='1Y')
    ic('sdii calculation complete!')

    newDset = xr.Dataset(
        {outKey: (("time","lat","lon"), sdiiDarr.data)},
        coords={
            "time": (('time'), sdiiDarr.time), #Times are datenums at "end of (fake) year"
            "lat": sdiiDarr['lat'],
            "lon": sdiiDarr['lon']
        }
    )
    newDset[outKey].attrs = prectMmDayDarr.attrs
    newDset[outKey].attrs['long_name'] = 'Simple intensity index'
    newDset[outKey].attrs['units'] = 'mm/day'

    inPcs = inFilePrect.split('/') #inFilePrect is the entire path to file
    inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
    strOut = inFn.replace(inKey, outKey) #Replace var name with extreme key
    outFile = outPath + strOut
    newDset.to_netcdf(outFile) #Save data
    ic(strOut, outFile, newDset) #icecream all useful parts of the output

def derive_prcptot(inFilePrect, outPath):
    ''' Derive total precip over period '''
    pIndices = pdex.indices(time_dim='time')
    inKey = 'PRECT'
    outKey = 'PRCPTOT'

    prectDset = xr.open_dataset(inFilePrect)
    prectDarr = prectDset[inKey]
    prectMmDarr = fcu.m_to_mm(prectDarr)
    prectMmDayDarr = fcu.persec_perday(prectMmDarr)
    ic('Beginning PRCPTOT calculation!')
    prectMmDayDarrChkLstTm = fpd.check_last_time(prectMmDayDarr)
    prcptotDarr = pIndices.prcptot(prectMmDayDarrChkLstTm, period='1Y', varname='PRECT')
    ic('PRCPTOT calculation complete!')

    newDset = xr.Dataset(
        {outKey: (("time","lat","lon"), prcptotDarr.data)},
        coords={
            "time": (('time'), prcptotDarr.time), #Times are datenums at "end of (fake) year"
            "lat": prcptotDarr['lat'],
            "lon": prcptotDarr['lon']
        }
    )
    newDset[outKey].attrs = prectMmDayDarr.attrs
    newDset[outKey].attrs['long_name'] = 'Total precip over period'
    newDset[outKey].attrs['units'] = 'mm'

    inPcs = inFilePrect.split('/') #inFilePrect is the entire path to file
    inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
    strOut = inFn.replace(inKey, outKey) #Replace var name with extreme key
    outFile = outPath + strOut
    newDset.to_netcdf(outFile) #Save data
    ic(strOut, outFile, newDset) #icecream all useful parts of the output

### Ocean extremes
def derive_mhw_presence_fromDefFile(inFileSst, outPath):
    ''' Calculate marine heatwave presence from daily SST data using the Hobday
    et al. 2016 definition. Returns binary classification (1=present 0=absent).
    Intensity data is not saved. Calculation relative to FIXED BASELINE:
    2010-2020 for GLENS, 2025-2035 for ARISE-SAI. Baseline definition file MUST
    be calculated first.'''
    mhwDefDict = {
        "defPath": '/Users/dhueholt/Documents/GLENS_data/extreme_MHW/definitionFiles/',
        "defPathCasper": '/glade/work/dhueholt/definitionFiles/',
        "defFile": 'baseDefFile_rlzMn_GLENS.nc',
        "defKey": 'mn_SST'
    }
    annSum = True
    inKey = 'SST'
    outKey = 'binary_mhw_pres'
    outKeyEv = 'binary_mhw_start'
    regOfInt = rlib.WesternAustraliaMHW_point() #Use point location in almost all circumstances

    sstDset = xr.open_dataset(inFileSst)
    sstDarr = sstDset[inKey]
    sstReg, locStr, _ = fpd.manage_area(sstDarr, regOfInt, areaAvgBool=True)
    sstRegFullTimes = sstReg.resample(time='1D').asfreq() #Add missing timesteps with NaN value
    sstRegDat = sstRegFullTimes.data.squeeze()

    # Set up ordinal time array (required by marineHeatWaves package)
    times = sstRegFullTimes.time.data
    tCftime = xr.cftime_range(times[0], times[len(times)-1])
    ordList = list()
    for tcf in tCftime:
        dtAc = tcf
        dtAcOrd = date(dtAc.year, dtAc.month, dtAc.day).toordinal()
        ordList.append(dtAcOrd)
    ordArr = np.array(ordList)

    # Load baseline definition file
    mhwDef = xr.open_dataset(mhwDefDict["defPath"] + mhwDefDict["defFile"])
    mhwDefTimes = mhwDef.time.data #POSSIBLY BLEHC
    mhwDefSst = mhwDef[mhwDefDict["defKey"]].data
    altClim = list([ordArr, mhwDefSst]) #Format required by mhws alternateClimatology feature

    mhwsDict, climDict = mhws.detect(ordArr, sstRegDat, alternateClimatology=altClim)

    # Make binary MHW presence/absence array
    binMhwPres = np.zeros(np.shape(ordArr)) #Initiate blank array of the proper size
    mhwStrtInd = mhwsDict['index_start']
    ic(mhwStrtInd)
    mhwDurInd = mhwsDict['duration']
    # ic(np.mean(mhwsDict['duration']))
    for actInd,startInd in enumerate(mhwStrtInd):
        actDur = mhwDurInd[actInd] #Duration of active MHW
        binMhwPres[startInd:startInd+actDur] = 1 #Times with active MHW get a 1
    ic(np.size(np.nonzero(binMhwPres))) # Number of MHW indices

    # Make binary MHW number of events array
    binMhwEv = np.zeros(np.shape(ordArr))
    ic(mhwStrtInd)
    binMhwEv[mhwStrtInd] = 1

    # Save xarray Dataset with all MHW data
    newDset = xr.Dataset(
        data_vars=dict(
            a=(["time"], binMhwPres),
            aa=(["time"], binMhwEv)
        ),
        coords=dict(
            time=times
        ),
        attrs=dict(description='MHW data')
    )
    newDset = newDset.rename_vars({'a': outKey, 'aa': outKeyEv})
    newDset[outKey].attrs['long_name'] = 'Binary presence/absence of MHWs'
    newDset[outKeyEv].attrs['long_name'] = 'Binary starting indices of independent MHWs'

    inPcs = inFileSst.split('/') #inFilePrect is the entire path to file
    inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
    outKeyReg = 'binary_mhw' + locStr
    strOut = inFn.replace(inKey, outKeyReg) #Replace var name with extreme key
    outFile = outPath + strOut

    if annSum == True:
        newDset = newDset.groupby("time.year").sum()
        outFile = outFile.replace('.nc', 'ann.nc')
        timeList = list()
        for yr in newDset.year:
            timeList.append(dtnl(yr,7,15,12,0,0,0)) #Year with standard fill values
        outDset = xr.Dataset(
            data_vars=dict(
                a=(["time"], newDset[outKey].data),
                aa=(["time"], newDset[outKeyEv].data)
            ),
            coords=dict(
                time=timeList
            ),
            attrs=dict(description='MHW data')
        )
        outDset = outDset.rename_vars({'a': outKey, 'aa': outKeyEv})
        outDset[outKey].attrs['long_name'] = 'Binary presence-absence of MHWs'
        outDset[outKeyEv].attrs['long_name'] = 'Binary starting indices of independent MHWs'
        outDset.to_netcdf(outFile) #Save data
    else:
        outDset = newDset.copy()
        outDset.to_netcdf(outFile) #Save data
    ic(strOut, outFile, outDset) #icecream all useful parts of the output

def derive_mhw_presence_reg(inFileSst, outPath):
    ''' Calculate marine heatwave presence from daily SST data using the Hobday
    et al. 2016 definition. Returns binary classification (1=present 0=absent).
    Intensity data is not saved.'''
    inKey = 'SST'
    outKey = 'binary_mhw_pres'
    outKeyEv = 'binary_mhw_start'
    regOfInt = rlib.WesternAustraliaMHW()

    sstDset = xr.open_dataset(inFileSst)
    sstDarr = sstDset[inKey]
    sstRegMn, locStr, _ = fpd.manage_area(sstDarr, regOfInt, areaAvgBool=True)
    sstRegMnDat = sstRegMn.data

    # Set up ordinal time array (required by marineHeatWaves package)
    times = sstDarr.time.data
    tCftime = xr.cftime_range(times[0], times[len(times)-1])
    ordList = list()
    for tcf in tCftime:
        dtAc = tcf
        dtAcOrd = date(dtAc.year, dtAc.month, dtAc.day).toordinal()
        ordList.append(dtAcOrd)
    ordArr = np.array(ordList)

    mhwsDict, climDict = mhws.detect(ordArr, sstRegMnDat)

    # Make binary MHW presence/absence array
    binMhwPres = np.zeros(np.shape(ordArr)) #Initiate blank array of the proper size
    mhwStrtInd = mhwsDict['index_start']
    ic(mhwStrtInd)
    mhwDurInd = mhwsDict['duration']
    # ic(np.mean(mhwsDict['duration']))
    for actInd,startInd in enumerate(mhwStrtInd):
        actDur = mhwDurInd[actInd] #Duration of active MHW
        binMhwPres[startInd:startInd+actDur] = 1 #Times with active MHW get a 1
    ic(np.size(np.nonzero(binMhwPres))) # Number of MHW indices

    # Make binary MHW number of events array
    binMhwEv = np.zeros(np.shape(ordArr))
    ic(mhwStrtInd)
    binMhwEv[mhwStrtInd] = 1

    # Save xarray Dataset with all MHW data
    newDset = xr.Dataset(
        data_vars=dict(
            a=(["time"], binMhwPres),
            aa=(["time"], binMhwEv)
        ),
        coords=dict(
            time=times
        ),
        attrs=dict(description='MHW data')
    )
    newDset = newDset.rename_vars({'a': outKey, 'aa': outKeyEv})
    newDset[outKey].attrs['long_name'] = 'Binary presence/absence of MHWs'
    newDset[outKeyEv].attrs['long_name'] = 'Binary starting indices of independent MHWs'

    inPcs = inFileSst.split('/') #inFilePrect is the entire path to file
    inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
    outKeyReg = outKey + outKeyEv + locStr
    strOut = inFn.replace(inKey, outKeyReg) #Replace var name with extreme key
    outFile = outPath + strOut
    # newDset.to_netcdf(outFile) #Save data
    ic(strOut, outFile, newDset) #icecream all useful parts of the output


### Ice metrics
def derive_icearea(inFileIcefrac, outPath):
    ''' Derive ice area from ice fraction. Ice area is the ice fraction of each
    grid cell multiplied by its area. This is not to be confused with ice
    EXTENT.'''
    inKey = 'ICEFRAC'
    outKey = 'ICEAREA'

    icefracDset = xr.open_dataset(inFileIcefrac)
    icefracDarr = icefracDset[inKey]
    icefracData = icefracDarr.data
    cellAreaDset = generate_gridcellarea(saveFlag=False)
    cellAreaDarr = cellAreaDset['grid_cell_area']
    cellAreaData = cellAreaDarr.data
    iceArea = np.multiply(icefracData,cellAreaData)

    newDset = xr.Dataset(
        {outKey: (("time","lat","lon"), iceArea)},
        coords={
            "time": (('time'), icefracDarr.time), #Times are datenums at "end of (fake) year"
            "lat": icefracDarr['lat'],
            "lon": icefracDarr['lon']
        }
    )
    newDset[outKey].attrs = icefracDarr.attrs
    newDset[outKey].attrs['long_name'] = 'Ice area'
    newDset[outKey].attrs['units'] = 'km^2'

    inPcs = inFileIcefrac.split('/') #inFilePrect is the entire path to file
    inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
    strOut = inFn.replace(inKey, outKey) #Replace var name with extreme key
    outFile = outPath + strOut
    newDset.to_netcdf(outFile) #Save data
    ic(strOut, outFile, newDset) #icecream all useful parts of the output

def derive_iceextent(inFileIcefrac, outPath):
    ''' Derive ice extent from ice fraction. Ice extent is the area of all grid
    cells with >15% ice fraction, and is frequently used due to its
    correspondence to observations. Not to be confused with ice AREA. '''
    inKey = 'ICEFRAC'
    outKey = 'ICEEXTENT'
    threshold = 0.15 #15% is convention used by e.g. NSIDC, EGU

    icefracDset = xr.open_dataset(inFileIcefrac)
    icefracDarr = icefracDset[inKey]
    icefracData = icefracDarr.data
    icefracTimes = icefracDarr.time
    cellAreaDset = generate_gridcellarea(saveFlag=False)
    cellAreaDarr = cellAreaDset['grid_cell_area']
    cellAreaData = cellAreaDarr.data
    # mask cellarea using icefrac > 0.15
    # need to build for each timestep! icefrac is monthly
    # then sum cellarea
    cellAreaMaskedList = list()
    for tc in np.arange(0,len(icefracTimes)):
        iceCoveredMask = icefracData[tc,:,:] >= threshold
        cellAreaDarrMasked = cellAreaDarr.copy()
        cellAreaDarrMasked.data[~iceCoveredMask] = np.nan
        cellAreaMaskedList.append(cellAreaDarrMasked)
    iceExtent = xr.concat(cellAreaMaskedList, icefracTimes)

    newDset = xr.Dataset(
        {outKey: (("time","lat","lon"), iceExtent)},
        coords={
            "time": (('time'), icefracDarr.time), #Times are datenums at "end of (fake) year"
            "lat": icefracDarr['lat'],
            "lon": icefracDarr['lon']
        }
    )
    newDset[outKey].attrs = icefracDarr.attrs
    newDset[outKey].attrs['long_name'] = 'Ice extent'
    newDset[outKey].attrs['units'] = 'km^2'

    inPcs = inFileIcefrac.split('/') #inFilePrect is the entire path to file
    inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
    strOut = inFn.replace(inKey, outKey) #Replace var name with extreme key
    outFile = outPath + strOut
    newDset.to_netcdf(outFile) #Save data
    ic(strOut, outFile, newDset) #icecream all useful parts of the output

### Misc
def generate_gridcellarea(saveFlag=False):
    ''' Generate a file with grid cell area data '''
    outFile = '/Users/dhueholt/Documents/GLENS_data/Misc/gridCellArea.nc'
    outKey = 'grid_cell_area'

    latNew = np.arange(-90,91,0.94240838)
    latNew = latNew[:-1]
    lonNew = np.arange(0,360,1.25)
    earthRad = 6371000 / 1000 #Earth's radius in km

    lonGrid,latGrid = np.meshgrid(lonNew, latNew)
    latGridRad = np.deg2rad(latGrid)
    dLat = np.gradient(latGrid)
    dLon = np.gradient(lonGrid)
    dDistY = dLat[0] * earthRad * np.pi / 180
    dDistX = (dLon[1]/180) * np.pi * earthRad * np.cos(latGridRad)

    gridCellArea = np.abs(dDistX * dDistY)

    newDset = xr.Dataset(
        {outKey: (("lat","lon"), gridCellArea)},
        coords={
            "lat": latNew,
            "lon": lonNew
        }
    )
    newDset[outKey].attrs['long_name'] = 'Grid cell area'
    newDset[outKey].attrs['units'] = 'km^2'
    if saveFlag:
        newDset.to_netcdf(outFile) #Save data

    return newDset
