'''wrap_derive_mhw_definitionFile
Defines MHW baseline for reference period at a given POINT location (typical) or
REGION (of dubious use). Clumsy for now as a first attempt, will be refactored
if project moves more in this direction.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import xarray as xr
import numpy as np
from datetime import date
import marineHeatWaves as mhws

import fun_process_data as fpd
import region_library as rlib

# Note that data ID keys are slightly different than most wrap_ scripts because
# these have raw-type filenames!
dataDict = {
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/daily_SST/',
    "idGlensCntrl": '*control*', #'*control*' or None
    "idGlensFdbck": None, #'*feedback*' or None
    "idArise": None, #'*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": None, # d'*BWSSP245*' or None
    "idS245Hist": None, #'*BWHIST*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc' #cesm_component_mask.nc
}
setDict = {
    "landmaskFlag": None, #None or 'land'
    "startIntvl": [2011,2030], #dg
    "endIntvl": [2041,2060], #dg
    "areaAvgBool": True, #True=mean, 'sum'=sum (e.g. ice extent), False=functionality that hasn't been tested in two years, and yes I know that's not a boolean
    "convert": None, #TUPLE of converter(s), or None if using default units
    "realization": 'ensplot',
    "regOfInt": (rlib.WesternAustraliaMHW_point(),),
}
outDict = {
    "outKeyRlzSeq": "rlzSeq_SST",
    "outKeyMn": "mn_SST",
    "savePath": '/Users/dhueholt/Documents/GLENS_data/extreme_MHW/definitionFiles/',
}

darrPointList = list()
for reg in setDict["regOfInt"]:
    darrList, cmnDict = fpd.call_to_open(dataDict, setDict)
    for darr in darrList:
        darrPoint, locStr, locTitleStr = fpd.manage_area(darr, reg, areaAvgBool=setDict["areaAvgBool"])
        darrPointList.append(darrPoint)
    ic(darrPointList)

    r0 = darrPointList[0].isel(realization=0).data.squeeze()
    r1 = darrPointList[0].isel(realization=1).data.squeeze()
    r2 = darrPointList[0].isel(realization=2).data.squeeze()
    r3 = darrPointList[0].isel(realization=3).data.squeeze()
    rmn = darrPointList[0].isel(realization=4).data.squeeze()

    times = darrPointList[0].time.data
    tCftime = xr.cftime_range(times[0], times[len(times)-1])
    ordList = list()
    for tcf in tCftime:
        dtAc = tcf
        dtAcOrd = date(dtAc.year, dtAc.month, dtAc.day).toordinal()
        ordList.append(dtAcOrd)
    ordArr = np.array(ordList)
    mhwsDictMn, climDictMn = mhws.detect(ordArr, rmn)

    rRlzSeq = np.concatenate((r0,r1,r2,r3))
    rlzSeqPer = len(times) * 4
    tCftimeRlzSeq = xr.cftime_range(start=times[0], periods=rlzSeqPer, freq='D')
    ordListRlzSeq = list()
    for tcfrs in tCftimeRlzSeq:
        dtAcRlzSeq = tcfrs
        dtAcOrdRlzSeq = date(dtAcRlzSeq.year, dtAcRlzSeq.month, dtAcRlzSeq.day).toordinal()
        ordListRlzSeq.append(dtAcOrdRlzSeq)
    ordArrRlzSeq = np.array(ordListRlzSeq)

    mhwsDictRlzSeq, climDictRlzSeq = mhws.detect(ordArrRlzSeq, rRlzSeq)
    ic(mhwsDictRlzSeq, climDictRlzSeq)

    rlzSeqClimDset = xr.Dataset(
        {outDict["outKeyRlzSeq"]: (("time"), rRlzSeq)},
        coords={
            "time": (('time'), ordArrRlzSeq), #Times are FAKE and DON'T EXIST
        }
    )
    rlzSeqClimDset[outDict["outKeyRlzSeq"]].attrs = darrPointList[0].attrs
    rlzSeqClimDset[outDict["outKeyRlzSeq"]].attrs['long_name'] = 'LENGTHENED SST'
    rlzSeqClimDset[outDict["outKeyRlzSeq"]].attrs['note'] = 'Realizations placed end to end with FAKE TIMES supplied, use for MHW calc'
    rlzSeqClimDset[outDict["outKeyRlzSeq"]].attrs['lat'] = reg["regLats"]
    rlzSeqClimDset[outDict["outKeyRlzSeq"]].attrs['lon'] = reg["regLons"]

    strOutRlzSeq = 'baseDefFile_rlzSeq_GLENS'
    outFileRlzSeq = outDict["savePath"] + strOutRlzSeq + '.nc'
    rlzSeqClimDset.to_netcdf(outFileRlzSeq) #Save data
    ic(strOutRlzSeq, outFileRlzSeq, rlzSeqClimDset) #icecream all useful parts of the output

    rlzMnClimDset = xr.Dataset(
        {outDict["outKeyMn"]: (("time"), rmn)},
        coords={
            "time": (('time'), ordArr),
        }
    )
    rlzMnClimDset[outDict["outKeyMn"]].attrs = darrPointList[0].attrs
    rlzMnClimDset[outDict["outKeyMn"]].attrs['long_name'] = 'Rlz mn SST'
    rlzMnClimDset[outDict["outKeyMn"]].attrs['note'] = 'Ensemble mean SST'
    rlzMnClimDset[outDict["outKeyMn"]].attrs['lat'] = reg["regLats"]
    rlzMnClimDset[outDict["outKeyMn"]].attrs['lon'] = reg["regLons"]

    strOutRlzMn = 'baseDefFile_rlzMn_GLENS'
    outFileRlzMn = outDict["savePath"] + strOutRlzMn + '.nc'
    rlzMnClimDset.to_netcdf(outFileRlzMn) #Save data
    ic(strOutRlzMn, outFileRlzMn, rlzMnClimDset) #icecream all useful parts of the output


# sys.exit('STOP')
