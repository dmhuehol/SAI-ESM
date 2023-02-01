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

rbd = {
    "beatNum": 6,
    "muteThresh": 6,
    "nRlz": 10 #ARISE=10, GLENS=21
}

holdRbst = list()
for t in np.arange(0, sd["trials"]):
    vec1 = np.random.uniform(low=0.0, high=1.0, size=rbd["nRlz"]) #"Feedback"
    vec2 = np.random.uniform(low=0.0, high=1.0, size=rbd["nRlz"]) #"Control"
    # ic(vec1, vec2)

    # This is the robustness calculation implemented for our random vectors!
    # Compare to calc_robustness_ecev
    # This is essentially a timeseries from a single point with nRlz time means
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
    # ic(rob)

    rbstAbv = fun_robustness.beat_rbst(rob["above"], beat=rbd["beatNum"])
    rbstBlw = fun_robustness.beat_rbst(rob["below"], beat=rbd["beatNum"])
    # ic(rbstAbv, rbstBlw)
    rbst = np.maximum(rbstAbv, rbstBlw)
    holdRbst.append(rbst)
    # ic(rbst)

# ic(holdRbst)
ic(fun_robustness.get_quantiles(holdRbst))
ic(st.mode(holdRbst))

# Want to know: what is a statistically significant value of holdRbst at 0.1
# 0.05 0.01 level
# t,p = st.ttest_ind(7, holdRbst)
# ic(t,p)
z = (rbd["muteThresh"]-np.mean(holdRbst))/np.std(holdRbst)
ic(z)
p = st.norm.sf(abs(-z))
ic(p)

ic(':D')
