# Code to double check region against GLENS output directly

def manage_area(darr, regionToPlot, areaAvgBool=True):
    ''' Manage area operations: obtain global, regional, or pointal output '''

    if regionToPlot == 'global':
        ic('global')
        locStr = 'global'
        locTitleStr = 'global'

        if areaAvgBool:
            latWeights = np.cos(np.deg2rad(darr['lat']))
            darrWght = darr.weighted(latWeights)
            darr = darrWght.mean(dim=['lat','lon'])

    elif isinstance(regionToPlot,dict):
        ic('regional')
        locStr = regionToPlot['regSaveStr']
        locTitleStr = regionToPlot['regSaveStr']

        lats = darr['lat'] #feedback and control are on same grid, fortunately
        lons = darr['lon']
        latMask = (lats>regionToPlot['regLats'][0]) & (lats<regionToPlot['regLats'][1])
        lonMask = (lons>regionToPlot['regLons'][0]) & (lons<regionToPlot['regLons'][1])
        darrBoxMask = darr[:,latMask,lonMask]

        ic(np.shape(lats))
        ic(np.shape(lons))
        LATTEST = lats[latMask]
        LONTEST = lons[lonMask]
        ic(np.shape(LATTEST))
        ic(np.shape(LONTEST))

        if areaAvgBool:
            darr = darrBoxMask.mean(dim=['lat','lon'])
        else:
            darr = darrBoxMask

    elif isinstance(regionToPlot,list):
        ic('point')
        darr = darr.sel(lat=regionToPlot[0], lon=regionToPlot[1], method="nearest")

        latStr = str(np.round_(darr.lat.data,decimals=2))
        lonStr = str(np.round_(darr.lon.data,decimals=2))
        locStr = latStr + '_' + lonStr
        locTitleStr = '(' + latStr + ',' + lonStr + ')'
    else:
        sys.exit('Invalid region! Check value for regionToPlot.')

    return darr, locStr, locTitleStr, LATTEST, LONTEST

def plot_single_basic_difference_globe(glensCntrlRlz, glensFdbckRlz, dataDict, setDict, outDict):
    ''' Plot 1 panel difference globe '''
    # Obtain levels
    glensCntrlLoi = pgf.obtain_levels(glensCntrlRlz, setDict["levOfInt"])
    glensFdbckLoi = pgf.obtain_levels(glensFdbckRlz, setDict["levOfInt"])

    # Average over years
    toiStart = dot.average_over_years(glensCntrlLoi, setDict["startIntvl"][0], setDict["startIntvl"][1]) # 2010-2019 is baseline, injection begins 2020
    toiEndCntrl = dot.average_over_years(glensCntrlLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    toiEndFdbck = dot.average_over_years(glensFdbckLoi, setDict["endIntvl"][0], setDict["endIntvl"][1])
    diffToiFdbck =  toiEndCntrl - toiEndFdbck

    # Unit conversion
    # diffToiFdbck = fcu.molmol_to_ppm(diffToiFdbck)

    # Plotting
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    cmap = cmocean.cm.balance
    cbVals = [-diffToiFdbck.quantile(0.99).data, diffToiFdbck.quantile(0.99).data]
    # cbVals = [-7, 7] #Override automatic colorbar minimum here
    md = pgf.meta_book(setDict, dataDict, labelsToPlot=None, glensCntrlLoi=glensCntrlLoi, glensFdbckRlz=glensFdbckRlz, cntrlToPlot=glensFdbckLoi)

    ## CHECKING REGIONS
    glensCntrlAoi, locStr, locTitleStr, LATTEST, LONTEST = pgf.manage_area(glensCntrlLoi, setDict["regOfInt"], setDict["areaAvgBool"])
    glensCntrlAoi = glensCntrlAoi.mean(dim='time')
    ic(glensCntrlAoi.data)
    CL = 0.
    mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
    plt.figure(figsize=(12, 2.73*2))
    ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index
    plt_tls.drawOnGlobe(ax, glensCntrlAoi.data, LATTEST, LONTEST, cmap='viridis', vmin=0, vmax=2, cbarBool=True, fastBool=True, extent='max')
    plt.title("AlaskaNorthwestCanada")

    # plt_tls.drawOnGlobe(ax, diffToiFdbck, LATTEST, LONTEST, cmap, vmin=cbVals[0], vmax=cbVals[1], cbarBool=True, fastBool=True, extent='max')
    # plt.title(md['lstDcd'] + ' ' + md['cntrlStr'] + '-' + md['fdbckStr'] + ' ' + md['levStr'] + ' ' + md['varStr'] + ' ' + 'Ens ' + str(setDict['realization']))
    # plt.title("2010-2019 Baseline - 2090-2099 SAI [50 0] ozone") #Override automatic title generation here

    savePrfx = '' #Easy modification for unique filename
    saveStr = md['varSve'] + '_' + md['levSve'] + '_' + md['lstDcd'] + '_' + md['ensStr'] + '_' + md['pid']['g1p'] + '_' + md['glbType']['fcStr']
    savename = outDict["savePath"] + savePrfx + saveStr + '.png'
    savename = outDict["savePath"] + 'testglobeAKNWC' + '.png'
    plt.savefig(savename, dpi=outDict["dpiVal"], bbox_inches='tight')
    ic(savename)
