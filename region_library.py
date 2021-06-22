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

### IPCC regions used in WG1-AR5

def AlaskaNorthwestCanada():
    regDict = {
        "regStr": 'Alaska/Northwest Canada',
        "regSaveStr": 'AKNWCan',
        "regLats": np.array([60, 72.6]),
        "regLons": np.array([192, 255])
    }

    return regDict

def CentralAsia():
    regDict = {
        "regStr": 'Central Asia',
        "regSaveStr": 'CentrlAsia',
        "regLats": np.array([30, 50]),
        "regLons": np.array([60, 75])
    }

    return regDict

def CanadaGreenlandIceland():
    regDict = {
        "regStr": 'Canada/Greenland/Iceland',
        "regSaveStr": 'CanGrnIce',
        "regLats": np.array([50, 85]),
        "regLons": np.array([255, 350])
    }

    return regDict

def EastAfrica():
    regDict = {
        "regStr": 'East Africa',
        "regSaveStr": 'EastAfrica',
        "regLats": np.array([-11.4, 15]),
        "regLons": np.array([25, 52])
    }

    return regDict

def EastAsia():
    regDict = {
        "regStr": 'East Asia',
        "regSaveStr": 'EastAsia',
        "regLats": np.array([20, 50]),
        "regLons": np.array([100, 145])
    }

    return regDict

def EastNorthAmerica():
    regDict = {
        "regStr": 'East N America',
        "regSaveStr": 'EstNAm',
        "regLats": np.array([25, 50]),
        "regLons": np.array([275, 300])
    }

    return regDict

def NorthAsia():
    regDict = {
        "regStr": 'North Asia',
        "regSaveStr": 'NrthAsia',
        "regLats": np.array([50, 70]),
        "regLons": np.array([40, 180])
    }

    return regDict

def NorthAustralia():
    regDict = {
        "regStr": 'North Australia',
        "regSaveStr": 'NrthAustrla',
        "regLats": np.array([-30, -10]),
        "regLons": np.array([110, 155])
    }

    return regDict

def NortheastBrazil():
    regDict = {
        "regStr": 'Northeast Brazil',
        "regSaveStr": 'NrthestBrazil',
        "regLats": np.array([-20, 0]),
        "regLons": np.array([310, 326])
    }

    return regDict

def SouthAustraliaNewZealand():
    regDict = {
        "regStr": 'South Australia/New Zealand',
        "regSaveStr": 'SAusNewZlnd',
        "regLats": np.array([-50, -30]),
        "regLons": np.array([110, 180])
    }

    return regDict

def SoutheastAsia():
    regDict = {
        "regStr": 'Southeast Asia',
        "regSaveStr": 'SEAsia',
        "regLats": np.array([-10, 20]),
        "regLons": np.array([95, 155])
    }

    return regDict

def TibetanPlateau():
    regDict = {
        "regStr": 'Tibetan Plateau',
        "regSaveStr": 'TibetPlat',
        "regLats": np.array([30, 50]),
        "regLons": np.array([75, 100])
    }

    return regDict

def WestAsia():
    regDict = {
        "regStr": 'West Asia',
        "regSaveStr": 'WestAsia',
        "regLats": np.array([15, 50]),
        "regLons": np.array([40, 60])
    }

    return regDict

def WestNorthAmerica():
    regDict = {
        "regStr": 'West N America',
        "regSaveStr": 'WstNAm',
        "regLats": np.array([28.6, 60]),
        "regLons": np.array([230, 255])
    }

    return regDict

def Antarctica():
    regDict = {
        "regStr": 'Antarctica',
        "regSaveStr": 'Antarctica',
        "regLats": np.array([-90, -50]),
        "regLons": np.array([0, 360])
    }

    return regDict

def Arctic():
    regDict = {
        "regStr": 'Arctic',
        "regSaveStr": 'Arctic',
        "regLats": np.array([67.5, 90]),
        "regLons": np.array([0, 360])
    }

    return regDict

def PacificIslandsRegion2():
    regDict = {
        "regStr": 'Pacific Islands Region[2]',
        "regSaveStr": 'PacIslReg2',
        "regLats": np.array([5, 25]),
        "regLons": np.array([155, 210])
    }

    return regDict

def PacificIslandsRegion3():
    regDict = {
        "regStr": 'Pacific Islands Region[3]',
        "regSaveStr": 'PacIslReg3',
        "regLats": np.array([-5, 5]),
        "regLons": np.array([155, 230])
    }

    return regDict

def SouthernTropicalPacific():
    regDict = {
        "regStr": 'Southern Tropical Pacific',
        "regSaveStr": 'STropPac',
        "regLats": np.array([-25,-5]),
        "regLons": np.array([155, 230])
    }

    return regDict

def WestIndianOcean():
    regDict = {
        "regStr": 'West Indian Ocean',
        "regSaveStr": 'WIndOcn',
        "regLats": np.array([-25, 5]),
        "regLons": np.array([52, 75])
    }

    return regDict

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
    # plt.show()
    plt.savefig('/Users/dhueholt/Documents/GLENS_fig/20210622_ipccRegions/' + region["regSaveStr"] + '.png', dpi=400)
