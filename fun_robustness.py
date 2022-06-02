''' fun_robustness
Contains functions to calculate the robustness of signals in data.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import numpy as np
import xarray as xr

def rbst_spread(cntrlDarr, fdbckDarr, spreadFlag='max'):
    ''' Evaluate robustness against ensemble spread by counting the number of
        Feedback members that fall outside of the Control max/min '''
    from numpy.random import default_rng
    rng = default_rng(13)
    vals = rng.integers(low=0, high=len(cntrlDarr.realization)-1, size=(192,288))
    ic(np.shape(vals))
    robustness = vals #number of ensemble members outside of spread

    return robustness

def mask_rbst(darr, robustness, nRlz):
    ''' Mask array based on robustness '''
    threshold = (nRlz)/2
    maskRobustness = robustness > threshold
    robustDarr = darr.copy()
    robustDarr = robustDarr.where(maskRobustness)

    return robustDarr
