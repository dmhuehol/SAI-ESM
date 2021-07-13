''' fun_convert_unit
Contains functions for unit conversions. Module is sectioned by type of variable,
i.e. temperature, chemistry, etc.

Written by Daniel Hueholt | July 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import xarray as xr
import numpy as np

### Chemistry

def molmol_to_ppm(darrMolmol):
    ''' Convert mol/mol to parts per million
    acmg.seas.harvard.edu/people/faculty/djj/book/bookchap1.html '''
    darrPpm = darrMolmol * 10**6
    darrPpm.attrs = darrMolmol.attrs
    darrPpm.attrs['units'] = 'ppm'

    return darrPpm

def molmol_to_ppb(darrMolmol):
    ''' Convert mol/mol to parts per billion
    acmg.seas.harvard.edu/people/faculty/djj/book/bookchap1.html '''
    darrPpb = darrMolmol * 10**9
    darrPpb.attrs = darrMolmol.attrs
    darrPpb.attrs['units'] = 'ppb'

    return darrPpb

def molmol_to_pptr(darrMolmol):
    ''' Convert mol/mol to parts per trillion
    acmg.seas.harvard.edu/people/faculty/djj/book/bookchap1.html '''
    darrPptr = darrMolmol * 10**12
    darrPptr.attrs = darrMolmol.attrs
    darrPptr.attrs['units'] = 'parts per trillion'

    return darrPptr

### Moisture
def kgkg_to_gkg(darrKgkg):
    ''' Convert kg/kg to g/kg '''
    darrGkg = darrKgkg * 1000
    darrGkg.attrs = darrKgkg.attrs
    darrGkg.attrs['units'] = 'g/kg'

    return darrGkg

### Precipitation
def m_to_cm(darrM):
    ''' Convert meters to centimeters '''
    darrCm = darrM * 100
    darrCm.attrs = darrM.attrs
    darrCm.attrs['units'] = darrCm.attrs['units'].replace("m/",'cm/')

    return darrCm

### Temperature
def kel_to_cel(darrKel):
    ''' Convert K to deg C '''
    darrCel = darrKel - 273.15
    darrCel.attrs = darrKel.attrs
    darrCel.attrs['units'] = 'deg C'

    return darrCel

### General
def persec_peryr(darrPerSec):
    ''' Convert per second to per year '''
    darrPerYr = darrPerSec * 3.154*10**7
    darrPerYr.attrs = darrPerSec.attrs
    darrPerYr.attrs['units'] = darrPerSec.attrs['units'].replace("/s",'/yr')
    if "rate" in darrPerYr.attrs['long_name']:
        darrPerYr.attrs['long_name'] = darrPerYr.attrs['long_name'].replace("rate",'')

    return darrPerYr
