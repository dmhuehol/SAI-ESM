''' fun_robustness
Contains functions to calculate the robustness of signals in data.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import numpy as np
import xarray as xr

def rbst_num_mn_ecev(cntrlDarr, fdbckDarr, spreadFlag='min', sprd=[15,20,5,10]):
    ''' "Each-Every" robustness. Evaluates robustness against ensemble spread
        by: for each Feedback time period, count the number of Control members
        that the given period is less/greater than. This is the only robustness
        function that fully accounts for NaN values. '''
    cntrlSprd = cntrlDarr.isel(time=np.arange(sprd[0],sprd[1]))
    cntrlSprdTimeMn = cntrlSprd.mean(dim='time')
    fdbckSprd = fdbckDarr.isel(time=np.arange(sprd[2],sprd[3]))
    fdbckSprdTimeMn = fdbckSprd.mean(dim='time')

    # Loop compares EACH feedback to EVERY control realization!
    countFdbckOutCntrl = np.full(np.shape(fdbckSprdTimeMn), np.nan)
    for rc,rv in enumerate(fdbckSprdTimeMn[:-1]): #Skip final ensemble member index, which is the ensemble mean
        if spreadFlag == 'max':
            fdbckOutCntrl = rv > cntrlSprdTimeMn
        elif spreadFlag == 'min':
            fdbckOutCntrl = rv < cntrlSprdTimeMn
        countFdbckOutCntrl[rc,:,:] = np.count_nonzero(fdbckOutCntrl.data, axis=2) #axis=2 is realization dimension
        nans = np.isnan(rv.data) #NaNs may be present, e.g. land area for ocean data
        countFdbckOutCntrl[rc,nans] = np.nan

    robustness = countFdbckOutCntrl #number of ensemble members outside of spread

    return robustness, nans

def rbst_num_mn_ecec(cntrlDarr, fdbckDarr, spreadFlag='min', sprd=[15,20,5,10]):
    ''' "Each-each" Evaluate robustness against ensemble spread by counting the
        number of Feedback members with a TIME MEAN less/greater than the
        corresponding Control member. Relies on there being a theoretically-
        motivated reason to pair a Feedback member of interest with a particular
        Control member. '''
    cntrlSprd = cntrlDarr.isel(time=np.arange(sprd[0],sprd[1]))
    cntrlSprdTimeMn = cntrlSprd.mean(dim='time')
    fdbckSprd = fdbckDarr.isel(time=np.arange(sprd[2],sprd[3]))
    fdbckSprdTimeMn = fdbckSprd.mean(dim='time')

    fdbckOutCntrl = fdbckSprdTimeMn < cntrlSprdTimeMn
    countOutCntrl = np.count_nonzero(fdbckOutCntrl.data, axis=0) #axis=0 is realization dimension
    ic(np.max(countOutCntrl),np.mean(countOutCntrl), np.median(countOutCntrl), np.mean(countOutCntrl))

    robustness = countOutCntrl

    return robustness

def rbst_spread(cntrlDarr, fdbckDarr, spreadFlag='min', sprd=[15,20,5,10]):
    ''' Evaluate robustness against ensemble spread by counting the number of
        Feedback members that fall outside of the Control max/min. This is an
        impractically stringent constraint, but is kept around for historical
        reasons. '''
    cntrlSprd = cntrlDarr.isel(time=np.arange(sprd[0],sprd[1]))

    if spreadFlag == 'max':
        cntrlSprdThresh = cntrlSprd.max()
    elif spreadFlag == 'min':
        cntrlSprdThresh = cntrlSprd.min()
    else:
        raise NotImplementedError('Check spread flag--must be "max" or "min"')

    fdbckSprd = fdbckDarr.isel(time=np.arange(sprd[2],sprd[3]))
    if spreadFlag == 'max':
        fdbckOutCntrl = fdbckSprd > cntrlSprdThresh
    elif spreadFlag == 'min':
        fdbckOutCntrl = fdbckSprd < cntrlSprdThresh
    fdbckOutCntrlYrs = np.count_nonzero(fdbckOutCntrl.data,axis=1)
    yrsThresh = 5
    fdbckOutCntrlYrsThresh = fdbckOutCntrlYrs >= yrsThresh
    fdbckOutCntrlYrsNum = np.count_nonzero(fdbckOutCntrlYrsThresh)
    # ic(len(fdbckOutCntrlYrs),fdbckOutCntrlYrsNum)

    robustness = fdbckOutCntrlYrsNum #number of ensemble members outside of spread
    ic(robustness)

    return robustness

def beat_rbst(rbstIn, beat=None):
    ''' Count how many members "beat" a given threshold value. '''
    if beat is None: #Default beat is greater than half of members
        numRlz = np.shape(rbstIn)[0]
        beat = numRlz / 2
        ic(beat)

    rlzAbvBeat = rbstIn > beat
    countGtBeat = np.count_nonzero(rlzAbvBeat, axis=0)
    rbstOut = countGtBeat
    ic(np.max(rbstOut),np.min(rbstOut),
       np.median(rbstOut), np.mean(rbstOut))

    return rbstOut


def mask_rbst(darr, robustness, nRlz):
    ''' Mask array based on robustness '''
    threshold = (nRlz)/2
    maskRobustness = robustness > threshold
    robustDarr = darr.copy()
    robustDarr = robustDarr.where(maskRobustness)

    return robustDarr

def get_quantiles(robustness):
    ''' Calculate and print quantiles for robustness array '''
    rbstQuant = {
        "0.99": np.nanquantile(robustness, 0.99),
        "0.9": np.nanquantile(robustness, 0.9),
        "0.7": np.nanquantile(robustness, 0.7),
        "0.5": np.nanquantile(robustness, 0.5),
        "0.3": np.nanquantile(robustness, 0.3),
        "0.2": np.nanquantile(robustness, 0.2),
        "0.1": np.nanquantile(robustness, 0.1),
        "0.05": np.nanquantile(robustness, 0.05),
        "0.01": np.nanquantile(robustness, 0.01)
    }
    ic(rbstQuant)
