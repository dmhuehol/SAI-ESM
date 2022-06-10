''' fun_robustness
Contains functions to calculate the robustness of signals in data.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import numpy as np
import xarray as xr

def rbst_spread(cntrlDarr, fdbckDarr, spreadFlag='min', sprd=[15,20,5,10]):
    ''' Evaluate robustness against ensemble spread by counting the number of
        Feedback members that fall outside of the Control max/min '''
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

def rbst_num_mn(cntrlDarr, fdbckDarr, spreadFlag='min', sprd=[15,20,5,10]):
    ''' Evaluate robustness against ensemble spread by counting the number of
        Feedback members with a TIME MEAN less/greater than time means in
        Control members.  '''
    cntrlSprd = cntrlDarr.isel(time=np.arange(sprd[0],sprd[1]))
    cntrlSprdTimeMn = cntrlSprd.mean(dim='time')
    fdbckSprd = fdbckDarr.isel(time=np.arange(sprd[2],sprd[3]))
    fdbckSprdTimeMn = fdbckSprd.mean(dim='time')

    fdbckOutCntrlRlz = np.full(np.shape(fdbckSprdTimeMn), np.nan)
    for rc,rv in enumerate(fdbckSprdTimeMn):
        if spreadFlag == 'max':
            fdbckOutCntrl = rv > cntrlSprdTimeMn
        elif spreadFlag == 'min':
            fdbckOutCntrl = rv < cntrlSprdTimeMn
        fdbckOutCntrlRlz[rc] = np.count_nonzero(fdbckOutCntrl.data)

    numThresh = 11 #n(Control members) to be less/greater than
    fdbckOutCntrlNumThresh = fdbckOutCntrlRlz >= numThresh
    fdbckOutCntrlNum = np.count_nonzero(fdbckOutCntrlNumThresh)

    robustness = fdbckOutCntrlNum #number of ensemble members outside of spread

    return robustness

def rbst_num_mn_ecec(cntrlDarr, fdbckDarr, spreadFlag='min', sprd=[15,20,5,10]):
    ''' Evaluate robustness against ensemble spread by counting the number of
        Feedback members with a TIME MEAN less/greater than Control.
        "Each-Each"'''
    cntrlSprd = cntrlDarr.isel(time=np.arange(sprd[0],sprd[1]))
    cntrlSprdTimeMn = cntrlSprd.mean(dim='time')
    fdbckSprd = fdbckDarr.isel(time=np.arange(sprd[2],sprd[3]))
    fdbckSprdTimeMn = fdbckSprd.mean(dim='time')

    fdbckOutCntrl = fdbckSprdTimeMn < cntrlSprdTimeMn
    countOutCntrl = np.count_nonzero(fdbckOutCntrl.data, axis=0)
    ic(np.max(countOutCntrl),np.mean(countOutCntrl), np.median(countOutCntrl), np.mean(countOutCntrl))

    # No thresholding for now
    # numThresh = 11 #n(Control members) to be less/greater than
    # fdbckOutCntrlNumThresh = fdbckOutCntrlRlz >= numThresh
    # fdbckOutCntrlNum = np.count_nonzero(fdbckOutCntrlNumThresh)
    #
    robustness = countOutCntrl

    return robustness

def rbst_num_mn_ecev(cntrlDarr, fdbckDarr, spreadFlag='min', sprd=[15,20,5,10]):
    ''' Evaluate robustness against ensemble spread by: for each Feedback
        time period, count the number of Control members that the Feedback
        is less/greater than.
        "Each-Every" '''
    cntrlSprd = cntrlDarr.isel(time=np.arange(sprd[0],sprd[1]))
    cntrlSprdTimeMn = cntrlSprd.mean(dim='time')
    fdbckSprd = fdbckDarr.isel(time=np.arange(sprd[2],sprd[3]))
    fdbckSprdTimeMn = fdbckSprd.mean(dim='time')

    # Loop compares EACH feedback to EVERY control realization!
    countFdbckOutCntrl = np.full(np.shape(fdbckSprdTimeMn), np.nan)
    for rc,rv in enumerate(fdbckSprdTimeMn[:-1]):
        if spreadFlag == 'max':
            fdbckOutCntrl = rv > cntrlSprdTimeMn
        elif spreadFlag == 'min':
            fdbckOutCntrl = rv < cntrlSprdTimeMn
        countFdbckOutCntrl[rc,:,:] = np.count_nonzero(fdbckOutCntrl.data, axis=2)

    robustness = countFdbckOutCntrl #number of ensemble members outside of spread

    return robustness

def thresh_count_rbst(rbstIn, thresh=None):
    ''' Threshold a robustness array on an input number of members, then count
        how many exceed. '''
    if thresh is None: #Default threshold is greater than half of members
        numRlz = np.shape(rbstIn)[0]
        thresh = numRlz / 2
        ic(thresh)

    rlzAbvThresh = rbstIn > thresh
    countGtThresh = np.count_nonzero(rlzAbvThresh, axis=0)
    rbstOut = countGtThresh
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
