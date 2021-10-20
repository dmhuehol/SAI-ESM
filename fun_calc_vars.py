''' fun_calc_vars
Functions to derive variables from root datasets.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import numpy as np
import xarray as xr

def sst_from_ptmp(darr):
    levName = 'z_t'
    levs = darr[levName].data
    levMask = levs == 500
    darrSurf = darr[:,levMask,:,:]

    return darrSurf

def uohc_from_ptmp(darr):
    ''' Sum upper 2000 meters of ocean potential temperature to obtain upper
        ocean heat content '''
    levs = darr['z_t'].data
    levMask = levs <= 2 * 10**5
    darrU2000 = darr[:,levMask,:,:]
    darrUohc = darrU2000.sum(dim='z_t')
    rhoSeawater = 1025 #kg/m3 OVERSIMPLIFIED, use data in future
    cpSeawater = 3850 #J/kgK OVERSIMPLIFIED
    intDistance = -2000 #m integration distance
    darrUohc = darrUohc * rhoSeawater * cpSeawater * intDistance
    vStr = list(['TEMP', 'UOHC'])
    darrUohc.attrs = darr.attrs
    darrUohc.attrs['units'] = 'J'
    darrUohc.attrs['long_name'] = 'Upper ocean heat content'
    darrUohc.attrs['outFile'] = darrUohc.attrs['inFile'].replace(vStr[0],vStr[1])
    # ic(darrUohc)
    dsetUohc = darrUohc.to_dataset(name=vStr[1]) #xarray doesn't like renaming DataArrays, but allows naming a new dataset
    darrUohc = dsetUohc[vStr[1]] #Then just re-extract the DataArray! It's an annoying but necessary hack

    return darrUohc
