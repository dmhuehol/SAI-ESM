''' fun_special_plot
Contains some special plotting functions for specific uses, e.g. plotting
colorbars alone or plotting quantiles vs. members for robustness.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys
import warnings

import cftime
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

import matplotlib.font_manager as fm
fontPath = '/Users/dhueholt/Library/Fonts/'  #Location of font files
for font in fm.findSystemFonts(fontPath):
    fm.fontManager.addfont(font)

import fun_calc_var as fcv
import fun_process_data as fpd
import fun_plot_tools as fpt

def save_colorbar(cbarDict, savePath, saveName, dpiVal=400):
    ''' Plots and saves a colorbar
        cbarDict should have a cmap, range, direction (horizontal vs. vertical),
        and label. Ex:
        cbarDict = {
            "cmap": cmocean.cm.delta,
            "range": [-15,15],
            "direction": 'vertical',
            "label": 'percent'
        }
    '''
    plt.rcParams.update({'font.size': 0})
    plt.rcParams.update({'font.family': 'Lato'})

    cbarRange = np.array([cbarDict["range"]])

    if cbarDict["direction"] == 'horizontal':
        plt.figure(figsize=(9,2.5))
        img = plt.imshow(cbarRange, cmap=cbarDict["cmap"])
        plt.gca().set_visible(False)
        colorAx = plt.axes([0.1,0.2,0.8,0.3])
        cb = plt.colorbar(orientation='horizontal', cax=colorAx, extend='both')
        for label in cb.ax.get_xticklabels():
            print(label)
            # label.set_fontproperties(FiraSansThin) #Set font
            label.set_fontsize(0) #Set font size
    elif cbarDict["direction"] == 'vertical':
        plt.figure(figsize=(2.5,9))
        img = plt.imshow(cbarRange, cmap=cbarDict["cmap"])
        plt.gca().set_visible(False)
        colorAx = plt.axes([0.1,0.2,0.2,0.6])
        cb = plt.colorbar(orientation='vertical', cax=colorAx)
        for label in cb.ax.get_yticklabels():
            print(label)
            # label.set_fontproperties(FiraSansThin)
            # label.set_fontsize(18)
    else:
        sys.exit('Direction must be either horizontal or vertical.')

    if cbarDict["label"] == '':
        cb.set_ticks([])
    # cb.set_label(cbarDict["label"], size='large')
    # cb.set_label(cbarDict["label"], size='large', fontproperties=FiraSansThin) # Set font
    plt.savefig(savePath + saveName + '.png', dpi=dpiVal)

def quantiles_vs_members(robustness, nRlz, savePath=None):
    ''' Plots quantiles vs number of members '''
    plt.rcParams.update({'font.family': 'Palanquin'})
    plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight

    qList = list()
    quantiles = np.linspace(0,1,500)
    for q in quantiles:
        qList.append(np.nanquantile(robustness,q))

    plt.figure(figsize=(9,6))
    plt.plot(quantiles, qList, color='#ff4d81', linewidth=2.5)

    plt.xlabel('Quantiles', fontsize=18, fontweight='light')
    plt.ylabel('Ensemble members', fontsize=18, fontweight='light')
    plt.xlim(-0.01,1.01)
    plt.xticks([0,0.2,0.4,0.6,0.8,1], fontsize=14)
    plt.ylim(0,nRlz+1)
    plt.yticks(np.arange(0,nRlz,3), fontsize=14)
    plt.title('Quantiles vs members: Annual mean 2m temperature', fontsize=22, fontweight='light')

    if savePath == None:
        plt.show()
    else:
        plt.savefig(savePath + 'QuantileVsMembers_AnnualMean2mTempARISE' + '.png')

def plot_rob_spaghetti_demo(darrList, dataDict, setDict, outDict):
    ''' Timeseries with ensemble members visualized as spaghetti plot. '''
    plotRlzMn = False # Plot ensemble mean in addition to every member
    setYear = [2030, 2045]
    timeSlice = slice(
        cftime.DatetimeNoLeap(setYear[0], 7, 15, 12, 0, 0, 0),
        cftime.DatetimeNoLeap(setYear[1], 7, 15, 12, 0, 0, 0))
    robTime = [2040,2044]
    robSlice = slice(
        cftime.DatetimeNoLeap(robTime[0], 7, 15, 12, 0, 0, 0),
        cftime.DatetimeNoLeap(robTime[1], 7, 15, 12, 0, 0, 0))

    # Plot timeseries
    plt.rcParams.update({'font.size': 14})
    plt.rcParams.update({'font.family': 'Lato'})
    plt.rcParams.update({'font.weight': 'normal'})
    fig,ax = plt.subplots()
    scnToPlot = list() # Make list of all scenarios to be plotted
    for scnDarr in darrList:
        rlzInScn = scnDarr['realization'].data # Number of members in scenario
        scnToPlot.append(scnDarr.scenario) # Add scenario to list
        if setDict["areaAvgBool"] == True: # Ens mean stored as last member
            rlzMn = scnDarr[len(rlzInScn)-1]
            rlzInScn = rlzInScn[:-1] # Don't double-plot ens mean
        elif setDict["areaAvgBool"] == 'sum': # Copy and take ens mean later
            rlzMn = scnDarr.copy()
        for rc in rlzInScn: # For each individual member
            rlzToi = scnDarr.sel(realization=rc, time=timeSlice)
            rlzLoi = fpd.obtain_levels(rlzToi, setDict["levOfInt"])
            rlzToPlot, locStr, locTitleStr = fpd.manage_area(
                rlzLoi, setDict["regOfInt"], areaAvgBool=setDict["areaAvgBool"])
            md = fpd.meta_book(
                setDict, dataDict, rlzToPlot, labelsToPlot=None) # Get metadata
            actCol, actLab = fpt.line_from_scenario(rlzToPlot.scenario, md)
            yrsToPlot = rlzToPlot['time'].dt.year.data
            plt.plot(
                yrsToPlot, rlzToPlot, color=actCol,
                linewidth=0.5, alpha=0.6) # Plot each member TIMESERIES
            rlzRobTm = scnDarr.sel(realization=rc, time=robSlice)
            rlzTmMn = rlzRobTm.mean(dim='time')
            rlzLoiTmMn = fpd.obtain_levels(rlzTmMn, setDict["levOfInt"])
            rlzTmMnToPlot, _, _ = fpd.manage_area(
                rlzLoiTmMn, setDict["regOfInt"],
                areaAvgBool=setDict["areaAvgBool"])
            plt.plot([robTime[0],robTime[1]], [rlzTmMnToPlot,rlzTmMnToPlot], color=actCol)

            if plotRlzMn:
                rlzToiMn = rlzMn.sel(time=timeSlice)
                rlzLoiMn = fpd.obtain_levels(rlzToiMn, setDict["levOfInt"])
                rlzAoiMn = rlzLoiMn.isel()
                if setDict["areaAvgBool"] == 'sum':
                    # Ens mean must be taken AFTER area average
                    rlzAoiMn = rlzAoiMn.mean(dim='realization')
                plt.plot(
                    yrsToPlot, rlzAoiMn, color=actCol, label=actLab,
                    linewidth=2.5) # Plot ensemble mean

    # Plot metadata and settings
    plt.yticks(ticks=setDict["yticks"],labels=None)
    b,t = plt.ylim(setDict["ylim"])
    fpt.plot_metaobjects(scnToPlot, fig, b, t, lw=0.4)
    plt.ylabel('')
    plt.xlim(setYear[0], setYear[1])
    plt.xticks(ticks=[setYear[0],setYear[0]+5,setYear[0]+10,setYear[0]+14],labels=None)

    # Save image
    savePrfx = 'ROB_'
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + str(setYear[0]) \
        + str(setYear[1]) + '_' + locStr + '_' + md['ensStr'] + '_' \
        + md['ensPid']['spg']
    savename = outDict["savePath"] + savePrfx + saveStr + '.pdf'
    plt.savefig(savename, bbox_inches='tight')
    plt.close()
    ic(savename)

def plot_rangeplot(loRlzList, loDataDictList, setDict, outDict):
    ''' Make rangeplot '''
    plt.rcParams.update({'font.size': 8})
    plt.rcParams.update({'font.family': 'Red Hat Display'})
    plt.rcParams.update({'font.weight': 'normal'})
    fig, ax = plt.subplots()
    scenarioList = list()
    dataList = list()
    
    for loc,rlzList in enumerate(loRlzList):
        if loc == 0: #Land, because there isn't a better way to track this
            disperseThrsh = 2
            disperseThrshMinus = -2
        elif loc == 1: #Ocean
            disperseThrsh = 7
            disperseThrshMinus = -7
        plotDict = fpt.make_scenario_dict(rlzList, setDict)
        for pscc,pscn in enumerate(plotDict.keys()):
            for itvl in plotDict[pscn].interval:
                ensMed = fcv.calc_weighted_med(plotDict[pscn])
                if setDict["magBool"] == True:
                    ensMedAbs = abs(ensMed)
                else:
                    ensMedAbs = ensMed.copy()
                ic(ensMedAbs)
                minEnsMed = np.min(ensMedAbs.data)
                maxEnsMed = np.max(ensMedAbs.data)
                mnEnsMed = np.mean(ensMedAbs.data)
                if mnEnsMed > 0: #Guess sign of response
                    gtThreshInd = ensMedAbs > disperseThrsh
                else:
                    gtThreshInd = ensMedAbs < disperseThrshMinus #greater than the MAGNITUDE, anyway
                ensMedAbsGtThr = ensMedAbs[gtThreshInd]
                ensMedAbsLteThr = ensMedAbs[~gtThreshInd]
                md = fpd.meta_book(
                        setDict, loDataDictList[loc], ensMed) # Get metadata
                actCol, actLab = fpt.line_from_scenario(ensMed.scenario, md)
                fcol = 'k'
                yVal, fcol = fpt.markers_from_scenario(
                    ensMed.scenario, loDataDictList[loc]["landmaskFlag"])
                if yVal is None:
                    yVal = pscc/2
                plt.plot(
                    [minEnsMed, maxEnsMed], [yVal, yVal], color=actCol, label=actLab,
                    linewidth=0.5)
                if (pscn == 'CESM2-WACCM:PreindustrialControl') | (pscn == 'LastMillennium'):
                    plt.scatter(
                        ensMedAbsGtThr.data, 
                        np.zeros(len(ensMedAbsGtThr.interval))+yVal,
                        color=actCol, facecolor=fcol)
                    plt.scatter(
                        ensMedAbsLteThr.data, 
                        np.zeros(len(ensMedAbsLteThr.interval))+yVal,
                        color=actCol, facecolor='none')    
                    # for imc, imv in enumerate(ensMedAbs.data):
                    #     plt.annotate(
                    #         str(imc), (imv, yVal+0.03), fontsize=8, color=actCol)
                    plt.scatter(
                        mnEnsMed, yVal, s=200, marker='|', color=actCol)
                else:
                    try:
                        plt.scatter(
                            ensMedAbsGtThr.sel(interval=itvl).data, 
                            np.zeros(len(ensMedAbsGtThr.realization))+yVal,
                            color=actCol, facecolor=fcol)
                        plt.scatter(
                            ensMedAbsLteThr.sel(interval=itvl).data, 
                            np.zeros(len(ensMedAbsLteThr.realization))+yVal,
                            color=actCol, facecolor='none')
                    except:
                        plt.scatter(
                            ensMedAbsGtThr.data, 
                            np.zeros(len(ensMedAbsGtThr.realization))+yVal,
                            color=actCol, facecolor=fcol)
                        plt.scatter(
                            ensMedAbsLteThr.data, 
                            np.zeros(len(ensMedAbsLteThr.realization))+yVal,
                            color=actCol, facecolor='none')
                        # for emc,emv in enumerate(ensMedAbs.data):
                        #     plt.annotate(
                        #         str(emc), (emv, yVal+0.03), fontsize=8, color=actCol)
                    plt.scatter(mnEnsMed, yVal, s=200, marker='|', color=actCol)
                if setDict["magBool"] == True:
                    plt.xlim([-0.1, 11])
                    plt.xticks([0, 2, 7, 10])
                else:
                    plt.xlim([-10.5,10.5])
                    plt.xticks([-10, -7, -2, 0, 2, 7, 10])
                
                plt.yticks([])
                b,t = plt.ylim()
                # plt.xlabel('Magnitude of climate speed (km/yr)')
        ax.spines[['left', 'top', 'right']].set_visible(False)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        
        savePrefix = ''
        saveStr = 'rangeplot' + '_' + md['varSve'][:11] + '_' + 'landocean' + \
            '_' + md['ensStr']
        if outDict["dpiVal"] == 'pdf':
            savename = outDict["savePath"] + savePrefix + saveStr + '.pdf'
            plt.savefig(savename, bbox_inches='tight')
        else:
            savename = outDict["savePath"] + savePrefix + saveStr + '.png'
            plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')    
                            
def plot_area_exposed(loRlzList, loDataDictList, setDict, outDict):
    ''' Make area exposed plot '''
    plt.rcParams.update({'font.size': 10})
    plt.rcParams.update({'font.family': 'Red Hat Display'})
    plt.rcParams.update({'font.weight': 'normal'})
    fig, ax = plt.subplots()
    paeTestList = list()
    # TODO: Rewrite to modularize (use a list factory)
    pae2List = list()
    pae5List = list()
    pae10List = list()
    pae30List = list()
    pae50List = list()
    pae2ListMn = list()
    pae5ListMn = list()
    pae10ListMn = list()
    pae30ListMn = list()
    pae50ListMn = list()
    for loc,rlzList in enumerate(loRlzList):
        plotDict = fpt.make_scenario_dict(rlzList, setDict)
        for pscc,pscn in enumerate(plotDict.keys()):
            plotIntRlzMn = plotDict[pscn].mean(dim=('interval','realization'))
            pae2Mn = fcv.calc_area_exposed(
                plotIntRlzMn, setDict, loDataDictList[loc], 2)
            pae5Mn = fcv.calc_area_exposed(
                plotIntRlzMn, setDict, loDataDictList[loc], 5)
            pae10Mn = fcv.calc_area_exposed(
                plotIntRlzMn, setDict, loDataDictList[loc], 10)
            pae30Mn = fcv.calc_area_exposed(
                plotIntRlzMn, setDict, loDataDictList[loc], 30)
            pae50Mn = fcv.calc_area_exposed(
                plotIntRlzMn, setDict, loDataDictList[loc], 50)
            pae2ListMn.append(pae2Mn.data)
            pae5ListMn.append(pae5Mn.data)
            pae10ListMn.append(pae10Mn.data)
            pae30ListMn.append(pae30Mn.data)
            pae50ListMn.append(pae50Mn.data)
            ic(plotIntRlzMn.scenario)
            ic(pae2ListMn, pae5ListMn, pae10ListMn, pae30ListMn, pae50ListMn)
            for itvl in plotDict[pscn].interval:
                plotIntDarr = plotDict[pscn].sel(interval=itvl)
                for rlz in plotDict[pscn].realization:
                    plotRlzDarr = plotIntDarr.sel(realization=rlz)
                    pae2 = fcv.calc_area_exposed(
                        plotRlzDarr, setDict, loDataDictList[loc], 2)
                    pae5 = fcv.calc_area_exposed(
                        plotRlzDarr, setDict, loDataDictList[loc], 5)
                    pae10 = fcv.calc_area_exposed(
                        plotRlzDarr, setDict, loDataDictList[loc], 10)
                    pae30 = fcv.calc_area_exposed(
                        plotRlzDarr, setDict, loDataDictList[loc], 30)
                    pae50 = fcv.calc_area_exposed(
                        plotRlzDarr, setDict, loDataDictList[loc], 50)
                    pae2List.append(pae2.data)
                    pae5List.append(pae5.data)
                    pae10List.append(pae10.data)
                    pae30List.append(pae30.data)
                    pae50List.append(pae50.data)
                    
                if  'CESM2-ARISE:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae2Mn], [0.5, 0.5], color='#12D0B2', linewidth=5.5)
                elif ('PreindustrialControl' in plotIntRlzMn.scenario) | ('LastMillennium' in plotIntRlzMn.scenario):
                    plt.plot([0, pae2Mn], [0.6, 0.6], color='#B8B8B8', linewidth=5.5)
                elif 'CESM2-ARISE:Control' in plotIntRlzMn.scenario:
                    plt.plot([0, pae2Mn], [0.7, 0.7], color='#F8A53D', linewidth=5.5)
                elif 'ARISE-DelayedStart:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae2Mn], [0.8, 0.8], color='#DDA2FB', linewidth=5.5)
                    
                if  'CESM2-ARISE:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae5Mn], [1.5, 1.5], color='#12D0B2', linewidth=5.5)
                elif ('PreindustrialControl' in plotIntRlzMn.scenario) | ('LastMillennium' in plotIntRlzMn.scenario):
                    plt.plot([0, pae5Mn], [1.6, 1.6], color='#B8B8B8', linewidth=5.5)
                elif 'CESM2-ARISE:Control' in plotIntRlzMn.scenario:
                    plt.plot([0, pae5Mn], [1.7, 1.7], color='#F8A53D', linewidth=5.5)
                elif 'ARISE-DelayedStart:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae5Mn], [1.8, 1.8], color='#DDA2FB', linewidth=5.5)
                    
                if  'CESM2-ARISE:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae10Mn], [2.5, 2.5], color='#12D0B2', linewidth=5.5)
                elif ('PreindustrialControl' in plotIntRlzMn.scenario) | ('LastMillennium' in plotIntRlzMn.scenario):
                    plt.plot([0, pae10Mn], [2.6, 2.6], color='#B8B8B8', linewidth=5.5)
                elif 'CESM2-ARISE:Control' in plotIntRlzMn.scenario:
                    plt.plot([0, pae10Mn], [2.7, 2.7], color='#F8A53D', linewidth=5.5)
                elif 'ARISE-DelayedStart:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae10Mn], [2.8, 2.8], color='#DDA2FB', linewidth=5.5)

                if  'CESM2-ARISE:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae30Mn], [3.5, 3.5], color='#12D0B2', linewidth=5.5)
                elif ('PreindustrialControl' in plotIntRlzMn.scenario) | ('LastMillennium' in plotIntRlzMn.scenario):
                    plt.plot([0, pae30Mn], [3.6, 3.6], color='#B8B8B8', linewidth=5.5)
                elif 'CESM2-ARISE:Control' in plotIntRlzMn.scenario:
                    plt.plot([0, pae30Mn], [3.7, 3.7], color='#F8A53D', linewidth=5.5)
                elif 'ARISE-DelayedStart:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae30Mn], [3.8, 3.8], color='#DDA2FB', linewidth=5.5)
                    
                if  'CESM2-ARISE:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae50Mn], [4.5, 4.5], color='#12D0B2', linewidth=5.5)
                elif ('PreindustrialControl' in plotIntRlzMn.scenario) | ('LastMillennium' in plotIntRlzMn.scenario):
                    plt.plot([0, pae50Mn], [4.6, 4.6], color='#B8B8B8', linewidth=5.5)
                elif 'CESM2-ARISE:Control' in plotIntRlzMn.scenario:
                    plt.plot([0, pae50Mn], [4.7, 4.7], color='#F8A53D', linewidth=5.5)
                elif 'ARISE-DelayedStart:Feedback' in plotIntRlzMn.scenario:
                    plt.plot([0, pae50Mn], [4.8, 4.8], color='#DDA2FB', linewidth=5.5)
                
                # ic(pae2List, pae5List, pae10List, pae30List, pae50List)
                pae2List = list()
                pae5List = list()
                pae10List = list()
                pae30List = list()
                pae50List = list()
                pae2ListMn = list()
                pae5ListMn = list()
                pae10ListMn = list()
                pae30ListMn = list()
                pae50ListMn = list()
        
    plt.xlim([0, 100])
    plt.ylim([0, 5])
    plt.yticks([])
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xticklabels([])
    
    md = fpd.meta_book(
            setDict, loDataDictList[loc], plotIntRlzMn) # Get metadata
    savePrefix = ''
    saveStr = 'aaexposed' + '_' + md['varSve'][:11] + '_' + setDict["landmaskFlag"] + \
        '_' + md['ensStr']
    if outDict["dpiVal"] == 'pdf':
        savename = outDict["savePath"] + savePrefix + saveStr + '.pdf'
        plt.savefig(savename, bbox_inches='tight')
    elif outDict["dpiVal"] == None:
        pass
    else:
        savename = outDict["savePath"] + savePrefix + saveStr + '.png'
        plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')    

def plot_warmrate_areaexposed(wrCsList, wrCsDictList, setDict, outDict):
    ''' Make warming rate vs area exposed (wrae) plot '''
    wrPlotDict = fpt.make_scenario_dict(wrCsList[0], setDict)
    wrDictList = wrCsDictList[0] # HARD CODED
    cspdDict = fpt.make_scenario_dict(wrCsList[1], setDict)
    cspdDictList = wrCsDictList[1] # HARD CODED
    paeThrshList = list()
    pwrList = list()
    piLine = None
    lmLine = None
    thrsh = 10 # FREE PARAMETER km/yr
    # ic(wrPlotDict.keys())
    
    plt.rcParams.update({'font.size': 10})
    plt.rcParams.update({'font.family': 'Red Hat Display'})
    plt.rcParams.update({'font.weight': 'normal'})
    fig, ax = plt.subplots()
    for pscc,pscn in enumerate(cspdDict.keys()):
        ic(pscn) # Troubleshooting
        plotIntRlzMn = cspdDict[pscn].mean(dim=('interval','realization'))
        paeThrshMn = fcv.calc_area_exposed(
            plotIntRlzMn, setDict, cspdDictList, thrsh)
        wrScnMn = wrPlotDict[pscn].mean(dim=('interval','realization'))
        # ic(wrScnMn)
        latWeights = np.cos(np.deg2rad(wrScnMn['lat']))
        darrMnWght = wrScnMn.weighted(latWeights)
        if 'ERA5' in pscn:
            wrScnGlobalMn = wrScnMn.mean(dim=['lat','lon'], skipna=True)
        else:
            wrScnGlobalMn = darrMnWght.mean(dim=['lat','lon'], skipna=True)
        fcol = fpt.colors_from_scenario(wrScnGlobalMn.scenario)
        # ic(paeThrshMn.data)
        if ('LastMillennium' not in pscn) & ('Preindustrial' not in pscn):
            plt.scatter(wrScnGlobalMn, paeThrshMn.data, s=80, color=fcol)
        for itvl in cspdDict[pscn].interval:
            for rlz in cspdDict[pscn].realization:
                cspdAct = cspdDict[pscn]
                plotRlzDarr = cspdAct.sel(interval=itvl).sel(realization=rlz)
                paeThrsh = fcv.calc_area_exposed(
                    plotRlzDarr, setDict, cspdDictList, thrsh)
                paeThrshList.append(paeThrsh)
                wrAct = wrPlotDict[pscn]
                wrScnRlz = wrAct.sel(interval=itvl).sel(realization=rlz)
                darrWghtRlz = wrScnRlz.weighted(latWeights)
                wrScnGlobal = darrWghtRlz.mean(dim=['lat','lon'], skipna=True)
                pwrList.append(wrScnGlobal)
                ic(paeThrsh.data)
                # plt.scatter( # Plot individual realizations, CAN BE MISLEADING
                    # wrScnGlobal, paeThrsh.data, color=fcol, s=20,
                    # facecolor='none')
        wrSpanRlz = np.max(np.abs(pwrList)) - np.min(np.abs(pwrList))
        pwrSpanRlz = [wrScnGlobalMn - wrSpanRlz/2, wrScnGlobalMn + wrSpanRlz/2]
        aeSpanRlz = np.max(np.abs(paeThrshList)) - np.min(np.abs(paeThrshList))
        paeSpanRlz = [
            paeThrshMn.data - aeSpanRlz / 2, paeThrshMn.data + aeSpanRlz / 2]
        if ('LastMillennium' not in pscn) & ('Preindustrial' not in pscn):
            plt.plot(
                pwrSpanRlz, [paeThrshMn.data, paeThrshMn.data], color=fcol)
            plt.plot(
                [wrScnGlobalMn.data, wrScnGlobalMn.data], paeSpanRlz, color=fcol)
        if 'Preindustrial' in pscn:
            piLine = paeSpanRlz#paeThrshMn.data + aeSpanRlz / 2
        if 'LastMillennium' in pscn:
            lmLine = paeSpanRlz#paeThrshMn.data + aeSpanRlz / 2
        paeThrshList = list()
        pwrList = list()
    xlim = [-0.1, 0.1] #0.02 land, 0.04 ocean
    plt.plot( # Zero line
        [0, 0], [0, 100], color='#242424', linewidth=0.4, linestyle='dashed',
        zorder=0.5)
    if 'CESM2-WACCM:PreindustrialControl' in setDict["plotScenarios"]:
        plt.plot( # Preindustrial boundary
            xlim, [piLine, piLine], color='#000000', linewidth=0.4, 
            linestyle='dashed', zorder=0.51)
    plt.plot( # Last Millennium boundary
        xlim, [lmLine, lmLine], color='#b8b8b8', linewidth=0.4, 
        linestyle='dashed', zorder=0.51)
    # ic(piLine)
    plt.xticks(
        [-0.1, -0.08, -0.06, -0.04, -0.02, 
            0, 0.02, 0.04, 0.06, 0.08, 0.1])
    plt.xlim(xlim)
    plt.ylim([0, 70])
    
    # plt.xlim(0, 0.1) #Supplementary Fig. 8
    # plt.ylim(10, 70) #Supplementary Fig. 8
    ax.spines[['top', 'right']].set_visible(False)
       
    md = fpd.meta_book(setDict, cspdDictList, plotIntRlzMn) # Get metadata
    savePrefix = '2_autoperiod10spc_'
    if setDict["landmaskFlag"] is None:
        setDict["landmaskFlag"] = 'none'
    saveStr = 'wrae' + '_' + str(thrsh) + 'kmyr' + '_' + \
        setDict["landmaskFlag"] + '_' + md['ensStr']
    if outDict["dpiVal"] == 'pdf':
        savename = outDict["savePath"] + savePrefix + saveStr + '.pdf'
        plt.savefig(savename, bbox_inches='tight')
    elif outDict["dpiVal"] == None:
        ic('Image generation disabled')
        pass
    else:
        savename = outDict["savePath"] + savePrefix + saveStr + '.png'
        plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')                    