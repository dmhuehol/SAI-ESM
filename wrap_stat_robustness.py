''' wrap_stat_robustness
Calculate statistical metrics for robustness.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import numpy as np
from scipy import stats as st

import fun_robustness

sd = {
    "trials": 10000
}

rbd = { # Commented values used in Hueholt et al. 2023
    "beatNum": 6, #ARISE=6, GLENS=11
    "muteThresh": 7, #ARISE=7, GLENS=15
    "nRlz": 10 #ARISE=10, GLENS=21
}

holdRbst = list()
for t in np.arange(0, sd["trials"]):
    vec1 = np.random.uniform(low=0.0, high=1.0, size=rbd["nRlz"]) #mimic "SAI"
    vec2 = np.random.uniform(low=0.0, high=1.0, size=rbd["nRlz"]) #mimic "no-SAI"
    # ic(vec1, vec2) #troubleshooting

    # This is the robustness calculation implemented for our random vectors!
    # Compare to calc_robustness_ecev
    # This corresponds to a timeseries from a single point with nRlz time means
    countVec1AbvVec2 = np.full(np.shape(vec1), np.nan)
    countVec1BlwVec2 = np.full(np.shape(vec1), np.nan)
    for rc, rv in enumerate(vec1):
        vec1AbvVec2 = rv > vec2
        countVec1AbvVec2[rc] = np.count_nonzero(vec1AbvVec2)
        vec1BlwVec2 = rv < vec2
        countVec1BlwVec2[rc] = np.count_nonzero(vec1BlwVec2)

    rob = {
        "above": countVec1AbvVec2,
        "below": countVec1BlwVec2,
    }
    # ic(rob) #troubleshooting

    rbstAbv = fun_robustness.beat_rbst(rob["above"], beat=rbd["beatNum"])
    rbstBlw = fun_robustness.beat_rbst(rob["below"], beat=rbd["beatNum"])
    # ic(rbstAbv, rbstBlw) #troubleshooting
    rbst = np.maximum(rbstAbv, rbstBlw)
    holdRbst.append(rbst)
    # ic(rbst) #troubleshooting

# ic(holdRbst)
ic(fun_robustness.get_quantiles(holdRbst))
ic(st.mode(holdRbst))

# We prefer to discuss robustness in terms of the percent of 
# randomly-generated values that the threshold falls outside. We believe
# this is the most intuitive way to understand the information that
# robustness conveys about the consistency of a response. However, it is 
# possible to formally describe this in terms of statistical significance
# using a Sign Test as in the following code. The thresholds we recommend
# for GLENS (beatNum=11, muteThresh=15) and ARISE (beatNum=6, muteThresh=7)
# correspond to significance at the p<0.1 value.
# See Hueholt et al. 2023 supplementary for more details.
td = {
    "p": (rbd["nRlz"] - rbd["beatNum"]) / rbd["nRlz"], # inherit from beat number
    "n": rbd["nRlz"], # number of tests
    "x": rbd["muteThresh"] # number of successes
}
ic(td)
holdProb = list()
for xc in np.arange(td["x"], td["n"]+1):
    prob = fun_robustness.binomial_test(td["p"], td["n"], xc)
    holdProb.append(prob)
pValueBinomialTest = np.sum(holdProb)
ic(pValueBinomialTest)