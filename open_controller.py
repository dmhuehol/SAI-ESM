''' open_controller 
Open and plot controller data.
Written by Daniel Hueholt and Charlotte Connolly '''
import glob
import sys

from icecream import ic
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

dataDict = {
    "dataPath": '/Users/dhueholt/Documents/allarise_data/controller/',
    # "idGlensFdbck": None,  # 'feedback_*' or None
    "idArise": '*DEFAULT*',  # '*DEFAULT*' or None
    "idUkesmArise": None, #'*arise-sai-1p5*' or None
    "idDelayedStart": None, # '*DELAYED*' or None
    "idArise1p0": None, # '*LOWER-0.5*' or None
    "mask": '/Users/dhueholt/Documents/Summery_Summary/cesm_atm_mask.nc', # Landmask file location (CESM)
    "maskUkesm": '/Users/dhueholt/Documents/UKESM_data/landmask/ukesm_binary_landmask.nc' #Landmask file location (UKESM)
}
outDict = {
    "savePath": '/Users/dhueholt/Documents/allarise_fig/20230424_controller/',
    "dpiVal": 400
}
scnList = list()
rlzList = list()
for dky in dataDict.keys():
    if 'id' in dky:
        # file = 'ControlLog_b.e21.BW.f09_g17.SSP245-TSMLT-GAUSS-DEFAULT.001.txt'
        try:
            inPath = dataDict["dataPath"]+dataDict[dky]
        except:
            pass
        inGlobs = sorted(glob.glob(inPath))
        for glf in inGlobs: 
            rawFile = np.array([x.split(' ') for x in open (glf).readlines()], dtype=object)
            header = rawFile[0]
            data = rawFile[1:]
            # ic(data)
            cntrlrTemplate = np.empty(shape=(len(data), 15))
            for dc, dv in enumerate(data):
                actData = [dlc for dlc in dv[:-1] if dlc != ''] #List comprehension to remove bad '' entries
                cntrlrTemplate[dc, :] = np.asarray(actData)
            # ic(cntrlrTemplate)
            # ic(cntrlrTemplate, np.shape(cntrlrTemplate))
            # cntrlr = xr.Dataset(
            #     {"controller": (
            #         ("time", "dT0", "sumdT0", "dT1", "sumdT1", "dT2", "sumdT2",
            #         "L0", "L1N", "L1S", "L2", "30STg", "15STg", "15NTg", "30NTg"),
            #         cntrlrTemplate)},
            #     coords={
            #         "time": cntrlrTemplate[:,0],
            #         "dT0": np.array([0]),
            #         "sumdT0": np.array([0]),
            #         "dT1": np.array([0]),
            #         "sumdT1": np.array([0]),
            #         "dT2": np.array([0]),
            #         "sumdT2": np.array([0]),
            #         "L0": np.array([0]),
            #         "L1N": np.array([0]),
            #         "L1S": np.array([0]),
            #         "L2": np.array([0]),
            #         "30STg": np.array([0]),
            #         "15STg": np.array([0]),
            #         "15NTg": np.array([0]),
            #         "30NTg": np.array([0]),
            #     }
            # )
            cntrlr = xr.Dataset(
                {"controller": (
                    ("time", "info"),
                    cntrlrTemplate)},
                coords={
                    "time": cntrlrTemplate[:,0],
                    "info": np.arange(0,15),
                }
            )
            rlzList.append(cntrlr)
        scnCntrlrDs = xr.concat(rlzList, dim='realization') #Controller data for a scenario
        ic(scnCntrlrDs)
        scnList.append(scnCntrlrDs)
        rlzList = list()


fig = plt.figure
actScnDa = scnList[0]["controller"]
inject30S = actScnDa.isel(info=11)
ic(inject30S, np.shape(inject30S))