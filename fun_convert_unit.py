# Do unit conversions for various variables.
#
# Written by Daniel Hueholt | May 2021
# Graduate Research Assistant at Colorado State University

from icecream import ic
import sys

import xarray as xr
import numpy as np

def molmol_to_ppm(darrMolmol):
# Convert mol/mol to parts per million
# acmg.seas.harvard.edu/people/faculty/djj/book/bookchap1.html
    darrPpm = darrMolmol * 10**6
    darrPpm.attrs['units'] = 'ppm'

    return darrPpm

def molmol_to_ppb(darrMolmol):
# Convert mol/mol to parts per billion
# acmg.seas.harvard.edu/people/faculty/djj/book/bookchap1.html
    darrPpb = darrMolmol * 10**9
    darrPpb.attrs['units'] = 'ppb'

    return darrPpb

def molmol_to_pptr(darrMolmol):
# Convert mol/mol to parts per trillion
# acmg.seas.harvard.edu/people/faculty/djj/book/bookchap1.html
    darrPptr = darrMolmol * 10**12
    darrPptr.attrs['units'] = 'parts per trillion'

    return darrPptr