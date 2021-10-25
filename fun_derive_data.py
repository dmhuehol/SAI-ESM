''' fun_derive_data
Functions to derive and save data from a base dataset as a new file.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import numpy as np
import xarray as xr
from cftime import DatetimeNoLeap as dtnl
import climdex.temperature as tdex
import climdex.precipitation as pdex

import fun_convert_unit as fcu

### Climdex extremes
def derive_annual_tropical_nights(inFileTrefht, outPath):
    tIndices = tdex.indices(time_dim='time')
    inKey = 'TREFHTMN'
    outKey = 'clxTR'

    trefhtDset = xr.open_dataset(inFileTrefht)
    trefhtDarr = trefhtDset[inKey]
    trefhtCelDarr = fcu.kel_to_cel(trefhtDarr)
    ic('Beginning clxTR calculation!')
    atnDarr = tIndices.annual_tropical_nights(trefhtCelDarr)
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
