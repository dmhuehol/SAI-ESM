''' wrap_derive_binarymhw_script
Derive data using multiprocessing for the OPERATION as opposed to the file I/O.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import glob
from multiprocessing import Process, Queue
import subprocess
from cftime import DatetimeNoLeap as dtnl
import cftime
from datetime import date
import numpy as np
import xarray as xr
import marineHeatWaves as mhws

import fun_process_data as fpd
import fun_derive_data as fdd
import region_library as rlib

inPath = '/Users/dhueholt/Documents/GLENS_data/daily_SST/' #LOCAL
# inPath = '/glade/scratch/dhueholt/daily_SST/selname/regrid/mergetime/' #CASPER
inToken = ['*control*','*feedback*']
# inTokens for raw: ['*control*','*feedback*','*SSP245-TSMLT*','*BWSSP245*','*BWHIST*']
# inTokens for cdo merged: ['*control*', '*feedback*', '*SSP245-TSMLT*', '*BWSSP245*', '*BWHIST*']
outPath = '/Users/dhueholt/Documents/GLENS_data/extreme_bmhw/' #LOCAL
# outPath = '/glade/scratch/dhueholt/extreme_bmhw/' #CASPER
poolStep = 25

for scen in inToken:
    theGlob = glob.glob(inPath+scen)
    for fc,fn in enumerate(theGlob):
        ic(fn)
        mhwDefDict = {
            "defPath": '/Users/dhueholt/Documents/GLENS_data/extreme_MHW/definitionFiles/',
            "defPathCasper": '/glade/work/dhueholt/definitionFiles/',
            "defFile": 'mhwDefsFile_GLENS_global.nc',
            "defKey": 'mn_sst'
        }
        ic(mhwDefDict["defFile"]) #Helps ensure the correct definitions file is used
        inKey = 'SST'
        outKey = 'binary_mhw_pres'
        regOfInt = rlib.Globe()

        sstDset = xr.open_dataset(fn)
        sstDarr = sstDset[inKey]
        sstReg, locStr, _ = fpd.manage_area(sstDarr, regOfInt, areaAvgBool=False)
        sstRegFullTimes = sstReg.resample(time='1D').asfreq() #Add missing timesteps with NaN value
        times = sstRegFullTimes.time.data
        ordArr = fpd.make_ord_array(times)

        try:
            mhwDef = xr.open_dataset(mhwDefDict["defPath"] + mhwDefDict["defFile"])
        except:
            mhwDef = xr.open_dataset(mhwDefDict["defPathCasper"] + mhwDefDict["defFile"])
        mhwDefSst = mhwDef[mhwDefDict["defKey"]].data
        mhwDefTimes = mhwDef.time.data
        altClim = list([mhwDefTimes, mhwDefSst]) #Format required by mhws alternateClimatology feature

        lats = sstRegFullTimes.lat.data
        lons = sstRegFullTimes.lon.data
        chunkList = list()
        lonInd = np.arange(0,len(lons))
        parChunk = int(len(lonInd)/poolStep)
        chunkLonInd = np.array_split(lonInd, parChunk)
        for ltc,ltv in enumerate(lats):
            for cli in chunkLonInd:
                resQue = Queue()
                for li in cli:
                    sstTs = sstRegFullTimes.isel(lat=ltc, lon=li)
                    if __name__== '__main__': #If statement required by multiprocessing
                    # Note Pool CANNOT work here as it runs inputs in sequence!
                        shard = Process(target=fdd.derive_binary_timeseries,
                                        args=(ordArr, sstTs, altClim, resQue))
                        shard.start()
                results = [resQue.get() for li in cli]
                for r in results:
                    chunkList.append(r)
                shard.join() #Forces chunks to complete. Must be AFTER results taken from queue
                shard.close() #Free up resources from PROCESSES
                resQue.close() #Free up resources from QUEUE
            if ltc % 10 == 0: #Cheap progress bar WITHIN file
                latProgress = (ltc+1)/len(lats) * 100
                latProgMsg = str(np.round(latProgress,1)) + '% latitudes completed'
                ic(latProgMsg)

        bmhwArr = np.transpose(np.asarray(chunkList)) #Transpose to standard (time, other) format
        bmhwReshape = np.reshape(bmhwArr, (np.shape(sstRegFullTimes))) #Reshape to (time, lat, lon)
        outDset = xr.Dataset(
            {outKey: (("time","lat","lon"), bmhwReshape)},
            coords={
                "time": (('time'), times),
                "lat": lats,
                "lon": lons
            }
        )
        outDset[outKey].attrs = sstRegFullTimes.attrs
        outDset[outKey].attrs['long_name'] = 'Binary MHW presence'
        outDset[outKey].attrs['units'] = 'd/yr'

        inPcs = fn.split('/') #inFilePrect is the entire path to file
        inFn = inPcs[len(inPcs)-1] #Filename is the last part of the path
        strOut = inFn.replace(inKey, outKey)
        outDset.to_netcdf(outPath + strOut)

        filesRemaining = len(theGlob) - fc - 1
        ic(filesRemaining) #Cheap "progress bar" for OVERALL progress
