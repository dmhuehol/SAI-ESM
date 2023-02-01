''' fun_basic_plot
Contains the difference globe plotting functions used in Hueholt et al. 2023
"Assessing Outcomes in Stratospheric Aerosol Injection Scenarios Shortly After
Deployment". The same three dictionary inputs (defining the input files, plot
settings, and output, respectively) are used by each function.

plot_single_basic_difference_globe: 1 panel difference globe with colorbar and
    automatic title
plot_paper_robust_globe: 1 panel robustness globe
    See Fig. S2 in Hueholt et al. 2023
plot_paper_panels_globe: 1 panel difference globe with aesthetics used in
    Hueholt et al. 2023. Figures 1, 3-8 are made of these figures combined with
    titles and other annotations added manually in Keynote.

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
    # Calculate robustness
    if setDict["robustnessBool"]:
        rbd, rbstns = fr.handle_robustness(rlzList)

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
    elif setDict["plotPanel"] == 'intiGLENS':
        intiGLENS = toiEnd['GLENS-SAI'] - toiEnd['RCP8.5']
        panel = intiGLENS
    elif setDict["plotPanel"] == 'intiARISE15':
        intiARISE15 = toiEnd['ARISE-SAI-1.5'] - toiEnd['SSP2-4.5']
        panel = intiARISE15
    else:
        blank = toiEnd['GLENS-SAI'].copy()
        blank.data = toiEnd['GLENS-SAI'] - toiEnd['GLENS-SAI']
        panel = blank
    # panel = snapS245 #Override setDict input
    panelStr = setDict["plotPanel"]

    # Plotting –– map
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    fig = plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance if setDict["cmap"] is None else setDict["cmap"]
    cbAuto = np.sort(
        [-np.nanquantile(panel.data,0.75), np.nanquantile(panel.data,0.75)])
    cbVals = cbAuto if setDict["cbVals"] is None else setDict["cbVals"]
    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
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
        muteThr = rbd["muteThr"]
        ic(muteThr)
        robustDarr = fr.mask_rbst(panel, rbstns, rbd["nRlz"], muteThr)
        fpt.drawOnGlobe(
            ax, robustDarr, lats, lons, cmap='Greys', vmin=cbVals[0],
            vmax=cbVals[1], cbarBool=False, fastBool=True, extent='max', addCyclicPoint=setDict["addCyclicPoint"], alph=0.65)
    # plt.title(" ") #No automatic title, 1-panel is used for custom figures

    # Plotting –– settings for output file
    savePrfx = '' #Easy modification for unique filename
    if 'ARISE' in panel.scenario:
        saveStr = panelStr + '_' + md['varSve'] + '_' + md['levSve'] + '_' \
                   + md['aLstDcd'] + '_' + md['ensStr'] + '_' \
                   + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    else:
        saveStr = panelStr + '_' + md['varSve'] + '_' + md['levSve'] + '_' \
                   + md['lstDcd'] + '_' + md['ensStr'] + '_' \
                   + md['pid']['g1p'] + '_' + md['glbType']['fcStr']

    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

## SPECIAL CASE GLOBE FUNCTIONS
def plot_paper_robust_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 1 panel robustness globe '''
    # Calculate robustness
    if setDict["robustnessBool"] is False:
        sys.exit('Cannot run robustness globe if robustness is False!')
    rbd, rbstns = fr.handle_robustness(rlzList)

    # Plotting –– quantiles vs members
    # fsp.quantiles_vs_members(rbstns, rbd["nRlz"], savePath=outDict["savePath"])
    # sys.exit('STOP') #Uncomment to plot only q vs m

    # Plotting –– map setup
    panel = rbstns
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    # mapProj = cartopy.crs.Orthographic(0, 90)#N: (0,90) S: (180,-90) #For polar variables
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    disc = cmasher.get_sub_cmap(cmasher.cm.ghostlight, 0, 1, N=rbd["nRlz"])
    # cmasher ghostlight, eclipse, cmyt pixel green, custom pink are good
    # cmasher apple to emphasize low robustness

    # Plotting –– image muting by altering discrete colorbar
    # Image muting is implemented for cmasher ghostlight and custom pink
    if rbd["muteQuThr"] is not None:
        muteThresh = int(np.ceil(np.nanquantile(rbstns, rbd["muteQuThr"]))) #Threshold to mute below
        ic(muteThresh)
        discColors = cmasher.take_cmap_colors(cmasher.cm.ghostlight,
                                              N=rbd["nRlz"], cmap_range=(0,1))
        muteList = fpt.mute_by_numbers(muteThresh)
        disc = mpl.colors.ListedColormap(muteList)

    # Plotting –– map
    cmap = disc if setDict["cmap"] is None else setDict["cmap"]
    cbAuto = [0, rbd["nRlz"]]
    cbVals = cbAuto if setDict["cbVals"] is None else setDict["cbVals"]
    cbarBool = False
    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    lats = rlzList[0].lat
    lons = rlzList[0].lon
    plt.rcParams.update({'font.family': 'Palanquin'})
    plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight

    cb,Im = fpt.drawOnGlobe(ax, panel, lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1],
                            cbarBool=cbarBool, fastBool=True, extent='max',
                            addCyclicPoint=setDict["addCyclicPoint"], alph=1)
    plt.title('') #Add manually in Keynote

    # Plotting –– settings for output file
    savePrfx = 'robust_' #Easy modification for unique filename
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_paper_panels_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 1 panel difference globes in loop for paper. This is
    (for now at least) designed to be REPLICABLE, not flexible. '''
    # Set up panels
    toiStart, toiEnd = fpt.make_panels(rlzList, setDict)
    diffToiR85 = toiEnd['RCP8.5'] - toiStart['RCP8.5']
    diffToiR85.attrs['pnl'] = 'RedGLENS'
    diffToiS245 = toiEnd['SSP2-4.5'] - toiStart['SSP2-4.5']
    diffToiS245.attrs['pnl'] = 'RedARISE'
    diffToiG12R85 = toiEnd['GLENS-SAI'] - toiStart['RCP8.5']
    diffToiG12R85.attrs['pnl'] = 'SnapGLENS'
    diffToiG15S245 = toiEnd['ARISE-SAI-1.5'] - toiStart['SSP2-4.5']
    diffToiG15S245.attrs['pnl'] = 'SnapARISE'
    intiG12R85 = toiEnd['GLENS-SAI'] - toiEnd['RCP8.5']
    intiG12R85.attrs['pnl'] = 'IntImpGLENS'
    intiG15S245 = toiEnd['ARISE-SAI-1.5'] - toiEnd['SSP2-4.5']
    intiG15S245.attrs['pnl'] = 'IntImpARISE'

    if 'TREFHT' in dataDict["dataPath"]: #NOTE HOW HORRIBLY HARDCODED THIS IS!
        panelList = (diffToiR85, diffToiS245, diffToiG12R85, diffToiG15S245, intiG12R85, intiG15S245)
    else:
        panelList = (diffToiG12R85, diffToiG15S245, intiG12R85, intiG15S245)

    # Plotting
    for panel in panelList:
        CL = 0.
        mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
        plt.figure(figsize=(12, 2.73*2))
        ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
        cmap = cmocean.cm.balance if setDict["cmap"] is None else setDict["cmap"]
        cbAuto = [-panel.quantile(0.75).data, panel.quantile(0.75).data]
        cbVals = cbAuto if setDict["cbVals"] is None else setDict["cbVals"]
        md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
        lats = rlzList[0].lat
        lons = rlzList[0].lon
        plt.rcParams.update({'font.family': 'Fira Sans'})
        plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight

        fpt.drawOnGlobe(ax, panel, lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=False, fastBool=True, extent='max', addCyclicPoint=setDict["addCyclicPoint"])

        savePrfx = '' #Easy modification for unique filename
        pnlId = ic(panel.attrs['pnl'])
        saveStr = pnlId + '_' + md['varSve'] + '_' + md['levSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
        savename = outDict["savePath"] + savePrfx + saveStr + '.png'
        # savename = outDict["savePath"] + savePrfx + saveStr + '.eps'
        plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
        # plt.savefig(savename, format='eps', dpi=10)
        plt.close()
        ic(savename)

def plot_single_robust_globe(rlzList, dataDict, setDict, outDict):
    ''' Plot 1 panel robustness globe '''
    # Calculate robustness
    if setDict["robustnessBool"] is False:
        sys.exit('Cannot run robustness globe if robustness is False!')
    rbd, rbstns = fr.handle_robustness(rlzList)

    # Plotting –– map setup
    panel = rbstns
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    # mapProj = cartopy.crs.Orthographic(0, 90) # N: (0,90) S: (180,-90) #Polar
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    disc = cmasher.get_sub_cmap(cmasher.cm.ghostlight, 0, 1, N=rbd["nRlz"]+1)

    # Plotting –– image muting by altering discrete colorbar
    # Implemented for cmasher ghostlight and custom pink
    if rbd["muteThr"] is not None:
        discColors = cmasher.take_cmap_colors(
            cmasher.cm.ghostlight, N=rbd["nRlz"], cmap_range=(0,1))
        if 'GLENS' in rbd["exp"]:
            muteList = fpt.mute_by_numbers(rbd["muteThr"])
        elif 'ARISE' in rbd["exp"]:
            muteList = fpt.mute_by_numbers_arise(rbd["muteThr"])
        disc = mpl.colors.ListedColormap(muteList)

    # Plotting –– map
    cmap = disc if setDict["cmap"] is None else setDict["cmap"]
    cbAuto = [0, rbd["nRlz"]]
    cbVals = cbAuto if setDict["cbVals"] is None else setDict["cbVals"]
    cbarBool = True
    md = fpd.meta_book(setDict, dataDict, rlzList[0], labelsToPlot=None)
    lats = rlzList[0].lat
    lons = rlzList[0].lon
    plt.rcParams.update({'font.family': 'Palanquin'})
    plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
    cb, Im = fpt.drawOnGlobe(
        ax, panel, lats, lons, cmap, vmin=cbVals[0], vmax=cbVals[1],
        cbarBool=cbarBool, fastBool=True, extent='max',
        addCyclicPoint=setDict["addCyclicPoint"], alph=1)
    if rbd["muteThr"] is not None: #Muting
        cb.set_ticks([0, rbd["muteThr"], rbd["nRlz"]])
    else: #No muting
        cb.set_ticks(np.linspace(0,rbd["nRlz"],3).astype(int))
    cb.ax.tick_params(labelsize=11)
    cb.set_label('number of members', size='small', fontweight='light')
    cb.remove()
    # plt.title('Count of SAI members outside of ' + str(rbd["beatNum"]) + ' no-SAI members: ' + md['varStr'], fontsize=16, fontweight='light') #No automatic title, 1-panel is used for custom figures
    plt.title('')

    # Plotting –– settings for output file
    savePrfx = 'robust_' #Easy modification for unique filename
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    plt.close()
    ic(savename)

