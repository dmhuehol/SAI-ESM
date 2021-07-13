''' region_library
Define useful regions for plotting to be called from other functions.
Latitudes are in deg N, longitudes are in 360-format deg E to match the GLENS
format.

Written by Daniel Hueholt | July 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import plotting_tools as plt_tls
import process_glens_fun as pgf

### IPCC regions used in WG1-AR5

def AlaskaNorthwestCanada():
    regDict = {
        "regStr": 'Alaska/Northwest Canada',
        "regSaveStr": 'AKNWCan',
        "regLats": np.array([60, 72.6]),
        "regLons": np.array([192, 255])
    }

    return regDict

def Amazon():
    regDict = {
        "regStr": 'Amazon',
        "regSaveStr": 'Amazon',
        "regLats": np.array([-20, -1.2, 11.4, 11.4, -20]),
        "regLons": np.array([293.6, 280.3, 291.2, 310, 310])
    }

    return regDict

def CentralAsia():
    regDict = {
        "regStr": 'Central Asia',
        "regSaveStr": 'CentralAsia',
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

def CentralAmericaMexico():
    regDict = {
        "regStr": 'Central America/Mexico',
        "regSaveStr": 'CAmMex',
        "regLats": np.array([11.4, -1.2, 28.6, 28.6]),
        "regLons": np.array([291.2, 280.3, 241.7, 269.7])
    }

    return regDict

def CentralEurope():
    regDict = {
        "regStr": 'Central Europe',
        "regSaveStr": 'CentralEurope',
        "regLats": (np.array([45, 48, 51, 45]), np.array([45, 51, 61.3, 45])),
        "regLons": (np.array([350, 350, 360, 360]), np.array([-1, -1, 40, 40]))
    }

    return regDict

def CentralNorthAmerica():
    regDict = {
        "regStr": 'Central N America',
        "regSaveStr": 'CenNAm',
        "regLats": np.array([28.6, 50]),
        "regLons": np.array([255, 275])
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
        "regSaveStr": 'NorthAsia',
        "regLats": np.array([50, 70]),
        "regLons": np.array([40, 180])
    }

    return regDict

def NorthAustralia():
    regDict = {
        "regStr": 'North Australia',
        "regSaveStr": 'NorthAustrla',
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

def NorthEurope():
    regDict = {
        "regStr": 'North Europe',
        "regSaveStr": 'NorthEurope',
        "regLats": (np.array([48, 75, 75, 51]), np.array([51, 75, 75, 61.3])),
        "regLons": (np.array([350, 350, 360, 360]), np.array([-1, -1, 40, 40]))
    }

    return regDict

def Sahara():
    regDict = {
        "regStr": 'Sahara',
        "regSaveStr": 'Sahara',
        "regLats": np.array([15, 30]),
        "regLons": np.array([340, 40])
    }

    return regDict

def SmallIslandsRegionsCaribbean():
    regDict = {
        "regStr": 'Small islands regions: Caribbean',
        "regSaveStr": 'SmIslCarbbn',
        "regLats": np.array([11.4, 25, 25, 11.4]),
        "regLons": np.array([291.2, 274.2, 300, 300])
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

def SouthAsia():
    regDict = {
        "regStr": 'South Asia',
        "regSaveStr": 'SouthAsia',
        "regLats": np.array([5, 30, 30, 20, 20, 5]),
        "regLons": np.array([60, 60, 100, 100, 95, 95])
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

def SoutheasternSouthAmerica():
    regDict = {
        "regStr": 'Southeastern South America',
        "regSaveStr": 'SESAm',
        "regLats": np.array([-20, -56.7, -56.7, -50, -20]),
        "regLons": np.array([320.6, 320.6, 292.7, 287.9, 293.6])
    }

    return regDict

def SouthernAfrica():
    regDict = {
        "regStr": 'Southern Africa',
        "regSaveStr": 'SthrnAfrica',
        "regLats": np.array([-35, -11.4]),
        "regLons": np.array([350, 52])
    }

    return regDict

def SouthEuropeMediterranean():
    regDict = {
        "regStr": 'South Europe/Mediterranean',
        "regSaveStr": 'SEurMed',
        "regLats": np.array([30, 45]),
        "regLons": np.array([350, 40])
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

def WestAfrica():
    regDict = {
        "regStr": 'West Africa',
        "regSaveStr": 'WestAfrica',
        "regLats": np.array([-11.4, 15]),
        "regLons": np.array([340, 25])
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

def WestCoastSouthAmerica():
    regDict = {
        "regStr": 'West Coast South America',
        "regSaveStr": 'WCstSAm',
        "regLats": np.array([-1.2, -20, -50, -56.7, -56.7, 0.5]),
        "regLons": np.array([280.3, 293.6, 287.9, 292.7, 278, 278])
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

def Tropics():
    regDict = {
        "regStr": 'Tropics',
        "regSaveStr": 'Tropics',
        "regLats": np.array([-23.45,23.45]),
        "regLons": np.array([0,360])
    }

    return regDict

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

def SouthernAfrica_hg():
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

### Atlases
def atlas_ipcc_wg1ar5():
    boxIPCC = (AlaskaNorthwestCanada(), CentralAsia(), CanadaGreenlandIceland(), CentralNorthAmerica(),
               EastAfrica(), EastAsia(), EastNorthAmerica(), NorthAsia(), NorthAustralia(), NortheastBrazil(),
               SouthAustraliaNewZealand(), SoutheastAsia(), TibetanPlateau(), WestAsia(), WestNorthAmerica(),
               Antarctica(), Arctic(), PacificIslandsRegion2(), PacificIslandsRegion3(), SouthernTropicalPacific(),
               WestIndianOcean())
    boxMerCrossIPCC = (SouthEuropeMediterranean(), SouthernAfrica(), Sahara(), WestAfrica())
    polyIPCC = (Amazon(), CentralAmericaMexico(), SmallIslandsRegionsCaribbean(), SouthAsia(),
                SoutheasternSouthAmerica(), WestCoastSouthAmerica())
    polyMerCrossIPCC = (CentralEurope(),NorthEurope())

    ipccWg1Ar5Dict = {
        "box": boxIPCC,
        "boxMerCross": boxMerCrossIPCC,
        "poly": polyIPCC,
        "polyMerCross": polyMerCrossIPCC,
        "allRegions": boxIPCC + boxMerCrossIPCC + polyIPCC + polyMerCrossIPCC
    }

    return ipccWg1Ar5Dict

### Functions

def west180_to_360(west180):
    ''' Convert from deg 180 to deg 360 (as used in GLENS) '''
    east360 = west180 % 360 #wrap to 360 degrees

    return east360

def test_region(region, colors, fig, ax):
    ''' Plots box on map to verify latitude/longitudes '''
    lats = np.arange(-90,91,1)
    lons = np.arange(0,360,1)

    if len(region['regLons'])>2: #non-rectangular region that does not cross Prime Meridian
        gridMask = pgf.make_polygon_mask(lats, lons, region['regLats'], region['regLons'])
        latsToPlot = lats
        lonsToPlot = lons
    elif isinstance(region['regLons'], tuple): #non-rectangular region that crosses Prime Meridian
        sgridMaskList = list()
        for sc in np.arange(0,len(region['regLons'])):
                sgridMask = pgf.make_polygon_mask(lats, lons, region['regLats'][sc], region['regLons'][sc])
                sgridMaskList.append(sgridMask)
        gridMask = np.logical_or.reduce(sgridMaskList)
        latsToPlot = lats
        lonsToPlot = lons
    else: #rectangular regions
        latMask = (lats>region['regLats'][0]) & (lats<region['regLats'][1])
        if region['regLons'][0] < region['regLons'][1]: #rectangle does not cross Prime Meridian
            lonMask = (lons>region['regLons'][0]) & (lons<region['regLons'][1])
        else: #rectangle crosses the Prime Meridian
            lonMask = (lons>region['regLons'][0]) | (lons<region['regLons'][1])
        latsToPlot = lats[latMask] #rectangular regions can just grab lat/lon of interest directly
        lonsToPlot = lons[lonMask]

    plotOnes = np.ones((len(latsToPlot),len(lonsToPlot)))
    try:
        plotOnes[~gridMask] = np.nan #non-rectangular regions must retain background grid and mask outside of region
    except:
        ic() #do nothing

    plt_tls.drawOnGlobe(ax, plotOnes, latsToPlot, lonsToPlot, cmap=colors, vmin=0, vmax=2, cbarBool=False, fastBool=True, extent='max')

    return fig, ax
