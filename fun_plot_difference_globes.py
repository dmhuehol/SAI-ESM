''' fun_plot_difference_globes
Contains the difference globe plotting functions used in Hueholt et al. 2023
"Assessing Outcomes in Stratospheric Aerosol Injection Scenarios Shortly After
Deployment". The same three dictionary inputs (defining the input files, plot
settings, and output, respectively) are used by each function.

plot_single_basic_difference_globe: 1 panel difference globe with colorbar and
    automatic title, see Fig. 1, 3-8 in Hueholt et al. 2023
plot_paper_robust_globe: 1 panel robustness globe
    See Fig. S2 in Hueholt et al. 2023

These functions are run from wrap_basicplots_script and
wrap_paperplots_basicplots_script. To make timeseries and plots that directly
display ensemble members, see wrap_ensplots_script and fun_ens_plot.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import numpy as np
import xarray as xr
xr.set_options(keep_attrs=True)
import cftime
import scipy.stats as stats
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn
import cmocean
import cmasher

import matplotlib.font_manager as fm
fontPath = '/Users/dhueholt/Library/Fonts/'  #Location of font files
for font in fm.findSystemFonts(fontPath):
    fm.fontManager.addfont(font)

import fun_convert_unit as fcu
import fun_process_data as fpd
import fun_plot_tools as fpt
import fun_robustness as fr
import fun_special_plot as fsp
import region_library as rlib

## DIFFERENCE GLOBES
def plot_single_basic_difference_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 1 panel difference globe '''
    # Set up panels
    toiStart, toiEnd = fpt.make_panels(rlzList, setDict)
    
    if setDict["plotPanel"] == 'snapR85':
        snapR85 = toiEnd['RCP8.5'] - toiStart['RCP8.5']
        panel = snapR85
    elif setDict["plotPanel"] == 'snapS245':
        snapS245 = toiEnd['SSP2-4.5'] - toiStart['SSP2-4.5']
        panel = snapS245
    elif setDict["plotPanel"] == 'snapGLENS':
        snapGLENS = toiEnd['GLENS-SAI'] - toiStart['RCP8.5']
        panel = snapGLENS
    elif setDict["plotPanel"] == 'snapARISE15':
        snapARISE15 = toiEnd['ARISE-SAI-1.5'] - toiStart['SSP2-4.5']
        panel = snapARISE15
    elif setDict["plotPanel"] == 'snapARISEDS':
        snapARISEDS = toiEnd['ARISE-SAI-1.5-DelayedStart'] - toiStart['SSP2-4.5']
        panel = snapARISEDS
    elif setDict["plotPanel"] == 'snapARISE10':
        snapARISE10 = toiEnd['ARISE-SAI-1.0'] - toiStart['SSP2-4.5']
        panel = snapARISE10
    elif setDict["plotPanel"] == 'intiGLENS':
        intiGLENS = toiEnd['GLENS-SAI'] - toiEnd['RCP8.5']
        panel = intiGLENS
    elif setDict["plotPanel"] == 'intiARISE15':
        intiARISE15 = toiEnd['ARISE-SAI-1.5'] - toiEnd['SSP2-4.5']
        panel = intiARISE15
    elif setDict["plotPanel"] == 'intiARISEDS':
        intiARISEDS = toiEnd['ARISE-SAI-1.5-DelayedStart'] - toiEnd['SSP2-4.5']
        panel = intiARISEDS
    elif setDict["plotPanel"] == 'intiARISE10':
        intiARISE10 = toiEnd['ARISE-SAI-1.0'] - toiEnd['SSP2-4.5']
        panel = intiARISE10
    elif setDict["plotPanel"] == 'snapUKS245':
        snapUKS245 = toiEnd['UKESM-SSP2-4.5'] - toiStart['UKESM-SSP2-4.5']
        panel = snapUKS245
    elif setDict["plotPanel"] == 'snapUKARISE15':
        snapUKS245 = toiEnd['UKESM-ARISE-SAI-1.5'] \
            - toiStart['UKESM-SSP2-4.5']
        panel = snapUKS245
    elif setDict["plotPanel"] == 'intiUKARISE15':
        intiUKARISE15 = toiEnd['UKESM-ARISE-SAI-1.5'] \
            - toiEnd['UKESM-SSP2-4.5']
        panel = intiUKARISE15
    else:
        blank = toiEnd['GLENS-SAI'].copy()
        blank.data = toiEnd['GLENS-SAI'] - toiEnd['GLENS-SAI']
        panel = blank
    # panel = snapS245 #Override setDict input
    panelStr = setDict["plotPanel"]

    # Plotting –– map
    CL = 0.
    # mapProj = cartopy.crs.Orthographic(180, -90)#N: (0,90) S: (180,-90) 
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    fig = plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance if setDict["cmap"] is None else setDict["cmap"]
    cbAuto = np.sort(
        [-np.nanquantile(panel.data,0.75), np.nanquantile(panel.data,0.75)])
    cbVals = cbAuto if setDict["cbVals"] is None else setDict["cbVals"]
    md = fpd.meta_book(setDict, dataDict, rlzList[0])
    lats = rlzList[0].lat
    lons = rlzList[0].lon
    plt.rcParams.update({'font.family': 'Open Sans'})
    plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
    fpt.drawOnGlobe(
        ax, panel, lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1],
        cbarBool=False, fastBool=True, extent='max',
        addCyclicPoint=setDict["addCyclicPoint"], alph=1)

    # Plotting –– image muting by adding separate layer of muted data
    if setDict["robustnessBool"]:
        if 'GLENS' in setDict["plotPanel"]:
            rbd = { #Settings for robustness calculation
                "beatNum": 11, #beat number is number of Control members to beat
                "muteThr": 15, #threshold to image mute; None to disable
                "sprd": [2025,2029],
                "nRlz": None, #Set automatically
                "exp": 'GLENS'
            }
        elif 'UK' in setDict["plotPanel"]: #TODO: improve scenario checking
            rbd = {
                "beatNum": 3,
                "muteThr": 4,
                "sprd": [2040,2044],
                "nRlz": None,
                "exp": 'UKARISE15'}
        elif 'ARISE15' in setDict["plotPanel"]:
            rbd = { #Settings for robustness calculation
                "beatNum": 6, #beat number is number of Control members to beat
                "muteThr": 7, #threshold to image mute; None to disable
                "sprd": [2040,2044],
                "nRlz": None, #Set automatically
                "exp": 'ARISE15'
            }
        else:
            rbd = None # Occurs for unknown scenario and won't work (by design)
        rbd, rbstns = fr.handle_robustness(rlzList, rbd)
        muteThr = rbd["muteThr"]
        ic(muteThr)
        robustDarr = fr.mask_rbst(panel, rbstns, rbd["nRlz"], muteThr)
        fpt.drawOnGlobe(
            ax, robustDarr, lats, lons, cmap='Greys', vmin=cbVals[0],
            vmax=cbVals[1], cbarBool=False, fastBool=True, extent='max', addCyclicPoint=setDict["addCyclicPoint"], alph=0.65)
    # plt.title(" ") #No automatic title, 1-panel is used for custom figures

    # Plotting –– settings for output file
    savePrfx = '' #Easy modification for unique filename
    if 'CESM2-ARISE' in panel.scenario:
        saveStr = panelStr + '_' + md['varSve'] + '_' + md['levSve'] \
            + '_' + md['lstDcd']["CESM2-ARISE"] + '_' + md['ensStr'] \
            + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    elif 'GLENS' in panel.scenario:
        saveStr = panelStr + '_' + md['varSve'] + '_' + md['levSve'] \
            + '_' + md['lstDcd']["GLENS"] + '_' + md['ensStr'] + '_' \
            + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    elif 'UKESM-ARISE' in panel.scenario:
        saveStr = panelStr + '_' + md['varSve'] + '_' + md['levSve'] \
            + '_' + md['lstDcd']["UKESM-ARISE"] + '_' + md['ensStr'] \
            + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    else:
        ic('Unable to generate saveStr for scenario')
        saveStr = panelStr

    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_single_robust_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 1 panel robustness globe '''
    # Calculate robustness
    if setDict["robustnessBool"]:
        if 'GLENS' in setDict["plotPanel"]:
            rbd = { #Settings for robustness calculation
                "beatNum": 11, #beat number is number of Control members to beat
                "muteThr": 15, #threshold to image mute; None to disable
                "sprd": [2025,2029],
                "nRlz": None, #Set automatically
                "exp": 'GLENS'
            }
        elif 'ARISE15' in setDict["plotPanel"]:
            rbd = { #Settings for robustness calculation
                "beatNum": 6, #beat number is number of Control members to beat
                "muteThr": 7, #threshold to image mute; None to disable
                "sprd": [2040,2044],
                "nRlz": None, #Set automatically
                "exp": 'ARISE15'   #Set automatically
            }
        else:
            rbd = None # Occurs for unknown scenario and won't work (by design)
        rbd, rbstns = fr.handle_robustness(rlzList, rbd)
    else:
        sys.exit('Cannot run robustness globe if robustness is False!')

    # Plotting –– map setup
    panel = rbstns
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    # mapProj = cartopy.crs.Orthographic(0, 90)#N: (0,90) S: (180,-90) #For polar variables
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    disc = cmasher.get_sub_cmap(cmasher.cm.ghostlight, 0, 1, N=rbd["nRlz"])

    # Plotting –– image muting by altering discrete colorbar
    # Image muting is implemented for cmasher ghostlight and custom pink
    discColors = cmasher.take_cmap_colors(cmasher.cm.ghostlight,
                                            N=rbd["nRlz"], cmap_range=(0,1))
    if rbd["exp"] == 'GLENS':
        muteList = fpt.mute_by_numbers(rbd["muteThr"])
    elif rbd["exp"] == 'ARISE15':
        muteList = fpt.mute_by_numbers_arise(rbd["muteThr"])
    disc = mpl.colors.ListedColormap(muteList)

    # Plotting –– map
    cmap = disc
    cbVals = [0, rbd["nRlz"]] #Robustness ranges up to size of the ensemble
    cbarBool = False
    md = fpd.meta_book(setDict, dataDict, rlzList[0])
    lats = rlzList[0].lat
    lons = rlzList[0].lon
    plt.rcParams.update({'font.family': 'Open Sans'})
    plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight

    cb,Im = fpt.drawOnGlobe(ax, panel, lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1],
                            cbarBool=cbarBool, fastBool=True, extent='max',
                            addCyclicPoint=setDict["addCyclicPoint"], alph=1)
    plt.title('') #Add manually in Keynote

    # Plotting –– settings for output file
    savePrfx = '' #Easy modification for unique filename
    saveStr = 'robust_' + setDict["plotPanel"] + '_' + md['varSve'] + '_' + \
              md['levSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + \
              md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)
