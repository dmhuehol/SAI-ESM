''' wrap_derive_data
must only have merged data to be derived form

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import glob
from multiprocessing import Process
import subprocess

import fun_derive_data as fdd

inPath = '/Users/dhueholt/Documents/GLENS_data/daily_TREFHT/mergetime/'
inToken = '*.nc'
outPath = '/Users/dhueholt/Documents/GLENS_data/clxTR/'
nProc = 2

theGlob = glob.glob(inPath+inToken)
for fc,fn in enumerate(theGlob):
    ic(fn)
    if __name__== '__main__': #If statement required by multiprocessing
        shard = Process(target=fdd.derive_annual_tropical_nights, args=(fn,outPath))
    if fc % nProc == 0 and fc != 0:
        shard.start()
        shard.join() #Forces nProc+1 processes to run to completion before beginning more
        shard.close() #Free up all associated resources
    else:
        shard.start()
    filesRemaining = len(theGlob) - fc - 1
    ic(filesRemaining) #Cheap "progress bar"
