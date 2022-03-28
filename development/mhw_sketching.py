'''
Attempting different threshold practices for defining marine heatwaves
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
    "dataPath": '/Users/dhueholt/Documents/GLENS_data/daily_SST/GLENS_reference/',
    "idGlensCntrl": '*control*', #'*control*' or None
    "idGlensFdbck": None, #'*feedback*' or None
    "idArise": None, #'*SSP245-TSMLT-GAUSS*' or None
    "idS245Cntrl": None, #'*BWSSP245*' or None
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
    "savePath": '/Users/dhueholt/Documents/GLENS_data/daily_SST/',
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

mhwsDict0, climDict0 = mhws.detect(ordArr, r0)
mhwsDict1, climDict1 = mhws.detect(ordArr, r1)
mhwsDict2, climDict2 = mhws.detect(ordArr, r2)
mhwsDict3, climDict3 = mhws.detect(ordArr, r3)
mhwsDictMn, climDictMn = mhws.detect(ordArr, rmn)

ic(mhwsDict0, climDict0)
ic(mhwsDict1, climDict1)
ic(mhwsDict2, climDict2)
ic(mhwsDict3, climDict3)
ic(mhwsDictMn, climDictMn)

ic(climDictMn['thresh'])
ic(len(climDictMn['thresh']))

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

# sys.exit('STOP')
