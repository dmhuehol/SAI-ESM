# Define useful regions for plotting to be called from other functions.
# Latitudes are in deg N, longitudes are in 360-format deg E to match the GLENS
# format.

import numpy as np

def GulfOfMexico():
    regDict = {
        "regStr": 'Gulf of Mexico',
        "regLats": np.array([19.5,30]),
        "regLons": np.array([263,280])
    }

    return regDict
