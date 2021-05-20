''' region_library
Define useful regions for plotting to be called from other functions.
Latitudes are in deg N, longitudes are in 360-format deg E to match the GLENS
format.

Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
'''

import numpy as np

### Oceans

## North Atlantic Basin

def GulfOfMexico():
    regDict = {
        "regStr": 'Gulf of Mexico',
        "regSaveStr": 'GulfOfMexico',
        "regLats": np.array([19.5,30]),
        "regLons": np.array([263,280])
    }

    return regDict

## South Indian Ocean

def LeeuwinCurrent():
    regDict = {
        "regStr": 'Leeuwin Current',
        "regSaveStr": 'LeeuwinCurrent',
        "regLats": np.array([-35,-22]),
        "regLons": np.array([108,115])
    }

    return regDict

## Southern Ocean

def DrakePassage():
    regDict = {
        "regStr": 'Drake Passage',
        "regSaveStr": 'DrakePassage',
        "regLats": np.array([-73,-55]),
        "regLons": np.array([290,305])
    }

    return regDict

def NoLandLatitude():
    regDict = {
        "regStr": 'NoLandLat60S',
        "regSaveStr": 'NoLandLat60S',
        "regLats": np.array([-65,-55]),
        "regLons": np.array([0,360])
    }

    return regDict

## Equatorial special regions

def Nino34():
    regDict = {
        "regStr": 'Nino 3.4',
        "regSaveStr": 'Nino34',
        "regLats": np.array([-5,5]),
        "regLons": np.array([190,240])
    }

    return regDict

### Terrestrial

## Africa

def SouthernAfrica():
    regDict = {
        "regStr": 'SouthernAfrica',
        "regSaveStr": 'SrnAfrica',
        "regLats": np.array([-36,16]),
        "regLons": np.array([12,37])
    }

    return regDict

## Antarctica

def AntarcticCircle():
    regDict = {
        "regStr": 'AntarcticCircle',
        "regSaveStr": 'AntrctcCrcl',
        "regLats": np.array([-90,-66.5]),
        "regLons": np.array([0,360])
    }

    return regDict

def Below50S():
    regDict = {
        "regStr": 'Below50S',
        "regSaveStr": 'Below50S',
        "regLats": np.array([-90,-50]),
        "regLons": np.array([0,360])
    }

    return regDict

## Asia

def NortheastAsia():
    regDict = {
        "regStr": 'NortheastAsia',
        "regSaveStr": 'NortheastAsia',
        "regLats": np.array([40,70]),
        "regLons": np.array([106,115])
    }

    return regDict

## Australia

def AustralianContinent():
    regDict = {
        "regStr": 'Australia',
        "regSaveStr": 'Australia',
        "regLats": np.array([-43,-12]),
        "regLons": np.array([111,115])
    }

    return regDict

## Europe

def EasternEurope():
    regDict = {
        "regStr": 'EasternEurope',
        "regSaveStr": 'EasternEurope',
        "regLats": np.array([36,56]),
        "regLons": np.array([14,36])
    }

    return regDict

## Placeholder

def Globe():
    regDict = {
        "regStr": 'Global',
        "regSaveStr": 'global',
        "regLats": np.array([np.nan, np.nan]),
        "regLons": np.array([np.nan, np.nan])
    }

    return regDict
