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

### Temperature
def kelvin_to_celsius(darrKel):
    darrCel = darrKel - 273.15
    darrCel.attrs = darrKel.attrs
    darrCel.attrs['units'] = 'celsius'

    return darrCel


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

### Temperature
def kel_to_cel(darrKel):
    ''' Convert K to deg C '''
    darrCel = darrKel - 273.15
    darrCel.attrs = darrKel.attrs
    darrCel.attrs['units'] = 'deg C'

    return darrCel
