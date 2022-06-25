''' fun_robustness
Contains functions to calculate the robustness of signals in data.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from icecream import ic
import sys

import cftime
import numpy as np
import xarray as xr

def handle_robustness(rlzList):
    ''' Handles robustness calculation '''
    rbd = { #Settings for robustness calculation
        "sprdFlag": 'outside', #calc based on above/below/outside of Control spread
        "beatNum": 11, #beat number is number of Control members to beat
        "muteQuThr": None, #threshold to image mute; None to disable
        "nRlz": None #Set automatically
    }

    # Select data of interest
    if len(rlzList) > 2:
        sys.exit('Robustness can only be run for GLENS OR ARISE, not both.')
    for s in rlzList:
        if 'Control' in s.scenario:
            actCntrlDarr = s
        elif 'Feedback' in s.scenario:
            actFdbckDarr = s
    rbd["nRlz"] = len(actCntrlDarr.realization)-1 #Skip ens mean at last index
    ic(rbd) #Show robustness dictionary for easy troubleshooting

    # Calculate robustness metric
    if 'GLENS' in actCntrlDarr.scenario:
        rbstEcEv, nans = rbst_num_mn_ecev(actCntrlDarr, actFdbckDarr, spreadFlag=rbd["sprdFlag"])
    elif 'ARISE' in actCntrlDarr.scenario:
        rbstEcEv, nans = rbst_num_mn_ecev(actCntrlDarr, actFdbckDarr, spreadFlag=rbd["sprdFlag"], sprd=[2040,2044])
    # Can use "outside" logic every time–no difference if "above" or "below" since
    # dictionary contains null arrays
    rbstnsAbv = beat_rbst(rbstEcEv["above"], beat=rbd["beatNum"])
    rbstnsBlw = beat_rbst(rbstEcEv["below"], beat=rbd["beatNum"])
    rbstns = np.maximum(rbstnsAbv, rbstnsBlw) #Composite of both above and below
    ic(rbstns) #Eyeball robustness array just for fun
    rbstns = rbstns.astype(np.float)
    rbstns[nans] = np.nan #NaNs from e.g. land area in ocean data

    return rbd, rbstns

def rbst_num_mn_ecev(cntrlDarr, fdbckDarr, spreadFlag='above', sprd=[2025,2029]):
    ''' "Each-Every" robustness. Evaluates robustness against ensemble spread
        by: for each Feedback time period, count the number of Control members
        that the given period is less/greater than. This is the only robustness
        function that fully accounts for NaN values. '''
    timeSlice = slice(cftime.DatetimeNoLeap(sprd[0], 7, 15, 12, 0, 0, 0),
                      cftime.DatetimeNoLeap(sprd[1], 7, 15, 12, 0, 0, 0))
    cntrlSprd = cntrlDarr.sel(time=timeSlice)
    cntrlSprdTimeMn = cntrlSprd.mean(dim='time')
    fdbckSprd = fdbckDarr.sel(time=timeSlice)
    fdbckSprdTimeMn = fdbckSprd.mean(dim='time')

    # Loop compares EACH feedback to EVERY control realization!
    countFdbckAbvCntrl = np.full(np.shape(fdbckSprdTimeMn), np.nan)
    countFdbckBlwCntrl = np.full(np.shape(fdbckSprdTimeMn), np.nan)
    for rc,rv in enumerate(fdbckSprdTimeMn[:-1]): #Skip final ensemble member index, which is the ensemble mean
        nans = np.isnan(rv.data) #NaNs may be present, e.g. land area for ocean data
        if spreadFlag == 'above':
            fdbckOutCntrl = rv > cntrlSprdTimeMn
            countFdbckAbvCntrl[rc,:,:] = np.count_nonzero(fdbckOutCntrl.data, axis=2) #axis=2 is realization dimension
            countFdbckAbvCntrl[rc,nans] = np.nan
        elif spreadFlag == 'below':
            fdbckOutCntrl = rv < cntrlSprdTimeMn
            countFdbckBlwCntrl[rc,:,:] = np.count_nonzero(fdbckOutCntrl.data, axis=2) #axis=2 is realization dimension
            countFdbckBlwCntrl[rc,nans] = np.nan
        elif spreadFlag == 'outside':
            fdbckAbvCntrl = rv > cntrlSprdTimeMn
            countFdbckAbvCntrl[rc,:,:] = np.count_nonzero(fdbckAbvCntrl.data, axis=2) #axis=2 is realization dimension
            countFdbckAbvCntrl[rc,nans] = np.nan
            fdbckBlwCntrl = rv < cntrlSprdTimeMn
            countFdbckBlwCntrl[rc,:,:] = np.count_nonzero(fdbckBlwCntrl.data, axis=2) #axis=2 is realization dimension
            countFdbckBlwCntrl[rc,nans] = np.nan

    robustness = {
        "above": countFdbckAbvCntrl,
        "below": countFdbckBlwCntrl,
    }

    return robustness, nans

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


def mask_rbst(darr, robustness, nRlz, threshold, spreadFlag):
    ''' Mask array based on robustness '''
    if (spreadFlag == 'below') or (spreadFlag == 'outside'):
        maskRobustness = robustness < threshold
    elif spreadFlag == 'above':
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

    return rbstQuant
