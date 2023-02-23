''' wrap_demo_gradient 
Demonstrate gradient calculation used for spatial gradients for climate
velocity and related variables.

Written by Daniel Hueholt '''

from icecream import ic
import sys

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as sndimg

from fun_calc_var import check_stats

# Make data of stacked vectors that increase by one for each element
# from left to right.
vec = np.array([0,1,2,3,4,5,6,7,8,9])
data = np.reshape(np.tile(vec, 9), (9,10))
ic(data, check_stats(data))

# 3x3 
nsSob = sndimg.sobel(data, axis=0, mode='reflect')
ewSob = sndimg.sobel(data, axis=1, mode='reflect')
ic(nsSob, check_stats(nsSob))
ic(ewSob, check_stats(ewSob))

# nsGrad = np.gradient(data, axis=0)
# ewGrad = np.gradient(data, axis=1)
# ic(nsGrad, check_stats(nsGrad))
# ic(ewGrad, check_stats(ewGrad))