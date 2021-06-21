''' region_library
Define useful regions for plotting to be called from other functions.
Latitudes are in deg N, longitudes are in 360-format deg E to match the GLENS
format.

Written by Daniel Hueholt | May 2021
Graduate Research Assistant at Colorado State University
'''

import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import plotting_tools as plt_tls


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

## Arctic

def ArcticCircle():
    regDict = {
        "regStr": 'ArcticCircle',
        "regSaveStr": 'ArctcCrcl',
        "regLats": np.array([66.5,90]),
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

### Hemispheres

def NorthernHemisphere():
    regDict = {
        "regStr": 'NorthernHemisphere',
        "regSaveStr": 'NHemi',
        "regLats": np.array([0,90]),
        "regLons": np.array([0,360])
    }

    return regDict

def SouthernHemisphere():
    regDict = {
        "regStr": 'SouthernHemisphere',
        "regSaveStr": 'SHemi',
        "regLats": np.array([-90,0]),
        "regLons": np.array([0,360])
    }

    return regDict

### Placeholder

def Globe():
    regDict = {
        "regStr": 'Global',
        "regSaveStr": 'global',
        "regLats": np.array([np.nan, np.nan]),
        "regLons": np.array([np.nan, np.nan])
    }

    return regDict

### Functions

def west180_to_360(west180):
    ''' Convert from deg 180 to deg 360 (as used in GLENS) '''
    east360 = west180 % 360 #wrap to 360 degrees

    return east360

def test_region(region):
    ''' Plots box on map to verify latitude/longitudes '''
    lats = np.linspace(region["regLats"][0], region["regLats"][1], 100)
    lons = np.linspace(region["regLons"][0], region["regLons"][1], 100)
    plotOnes = np.ones((len(lats), len(lons)))

    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    plt_tls.drawOnGlobe(ax, plotOnes, lats, lons, cmap='viridis', vmin=0, vmax=2, cbarBool=True, fastBool=True, extent='max')
    plt.title(region["regStr"])
    plt.show()
