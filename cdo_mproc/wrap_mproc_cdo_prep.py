from icecream import ic
import sys
import time

from glob import glob
from multiprocessing import Process
import numpy as np
import subprocess

def cdo_annual(IN_PATH, IN_TOKEN, OUT_PATH):
    subprocess.run(['sh', '/glade/u/home/dhueholt/glens_scripts/GLENS/do_cdo_prep.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def return_emem_list(inType):
    if inType == "GLENS":
        EMEM=list([
        ".001.",
        ".002.",
        ".003.",
        ".004.",
        ".005.",
        ".006.",
        ".007.",
        ".008.",
        ".009.",
        ".010.",
        ".011.",
        ".012.",
        ".013.",
        ".014.",
        ".015.",
        ".016.",
        ".017.",
        ".018.",
        ".019.",
        ".020.",
        ".021."]
        )
    elif inType == "SCIRIS":
        EMEM=list([
        ".001.",
        ".002.",
        ".003.",
        ".004.",
        ".005.",
        ".006.",
        ".007.",
        ".008.",
        ".009.",
        ".010."]
        )
    elif inType == "CMIP6":
        EMEM=list([
        "r1",
        "r2",
        "r3",
        "r4",
        "r5"]
        )
    elif inType == "historical":
        EMEM=list([
        ".001.",
        ".002.",
        ".003."]
        )
    else:
        sys.exit('Unrecognized type!')

    identOfInt = '*' #*, or something unique if folder of data has multiple variables
    EMEM = [r + identOfInt for r in EMEM]

    return EMEM

EMEM = return_emem_list('GLENS')
nProc = 5

# Shell inputs
IN_PATH = '/glade/scratch/dhueholt/ssp245/monthly_TSA/'
IN_TOKEN = ['*control*','*feedback*','*ssp245*'] #['*control*','*feedback*','*ssp245*','*CESM2-WACCM*'] #GLENS, SCIRIS, CMIP6
OUT_PATH = '/glade/scratch/dhueholt/ssp245/annual_TSA/'

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
