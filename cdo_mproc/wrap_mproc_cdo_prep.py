from icecream import ic
import sys
import time

from glob import glob
from multiprocessing import Process
import numpy as np
import subprocess

inPath = '/Users/dhueholt/Documents/GLENS_data/annual_PREC/*_00*PRECC*.nc'
inList = sorted(glob(inPath))
# ic(inList)
# sys.exit('STOP')

def cdo_annual(IN_PATH, IN_TOKEN, OUT_PATH):
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/GLENS/GLENS/do_cdo_prep.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

EMEM=list([
"001",
"002",
"003",
"004",
"005",
"006",
"007",
"008",
"009",
"010",
"011",
"012",
"013",
"014",
"015",
"016",
"017",
"018",
"019",
"020",
"021"]
)
identOfInt = '.pop.h.SSH.*'
EMEM = [r + identOfInt for r in EMEM]
nProc = 5

# Shell inputs
IN_PATH = '/glade/scratch/dhueholt/monthly_IAGE/'
IN_TOKEN = ['*control*','*feedback*']
OUT_PATH = '/glade/scratch/dhueholt/annual_IAGE/'

if __name__== '__main__':
        lengthFiles = np.size(EMEM)
        for scen in IN_TOKEN:
            for rc,rv in enumerate(EMEM):
                # Instantiate a new process
                p = Process(target=cdo_annual,args=(IN_PATH, scen+rv, OUT_PATH))
                if rc % nProc == 0 and rc != 0:
                    # Run nProc number of processes at a time
                    p.start()
                    p.join()
                    p.close() #End the process
                else:
                    p.start()
                numleft = lengthFiles - rc - 1
                ic(numleft) #Displays how many files remain to be processed in current job
