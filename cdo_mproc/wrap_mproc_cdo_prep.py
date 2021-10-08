'''  wrap_mproc_cdo_prep
Runs functions to process data into friendlier forms, i.e. selecting months or
calculating annual mean values. Uses multiprocessing for efficiency.

After running, output data is ready for regridding (wrap_ocean_script) or
plotting (wrap_basicplots_script, wrap_ensplots_script).

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

from glob import glob
from multiprocessing import Process
import numpy as np
import subprocess

def cdo_annualmean(IN_PATH, IN_TOKEN, OUT_PATH):
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/GLENS/GLENS/cdo_mproc/cdo_bees/do_cdo_annualmean.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_seleimnsn(IN_PATH, IN_TOKEN, OUT_PATH):
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/GLENS/GLENS/cdo_mproc/cdo_bees/do_cdo_seleimnsn.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_seleirainy(IN_PATH, IN_TOKEN, OUT_PATH):
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/GLENS/GLENS/cdo_mproc/cdo_bees/do_cdo_seleirainy.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_selfeb(IN_PATH, IN_TOKEN, OUT_PATH):
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/GLENS/GLENS/cdo_mproc/cdo_bees/do_cdo_selfeb.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_sellevel(IN_PATH, IN_TOKEN, OUT_PATH, IN_LEV=500):
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/GLENS/GLENS/cdo_mproc/cdo_bees/do_cdo_select_lev.sh', IN_PATH, IN_TOKEN, IN_LEV, OUT_PATH])
    return None

def cdo_selsept(IN_PATH, IN_TOKEN, OUT_PATH):
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/GLENS/GLENS/cdo_mproc/cdo_bees/do_cdo_selsept.sh', IN_PATH, IN_TOKEN, OUT_PATH])
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
    elif inType == "CMIP6raw":
        EMEM=list([
        ".001.",
        ".002.",
        ".003.",
        ".004.",
        ".005."]
        )
    elif inType == "historical":
        EMEM=list([
        ".001."]
        # ".002.",
        # ".003."]
        )
    else:
        sys.exit('Unrecognized type!')

    identOfInt = '*' #*, or something unique if folder of data has multiple variables
    EMEM = [r + identOfInt for r in EMEM]

    return EMEM

EMEM = return_emem_list('GLENS')
nProc = 5

# Shell inputs
IN_PATH = '/Users/dhueholt/Documents/GLENS_data/monthly_500TEMP/'
IN_TOKEN = ['*SSP245-TSMLT*','*BWSSP245*','*BWHIST*'] #['*control*','*feedback*','*SSP245-TSMLT*','*BWSSP245*','*BWHIST*'] #GLENS, ARISE, CESM2-WACCM, historical
OUT_PATH = '/Users/dhueholt/Documents/GLENS_data/annual_500TEMP/'

if __name__== '__main__':
        lengthFiles = np.size(EMEM)
        for scen in IN_TOKEN:
            for rc,rv in enumerate(EMEM):
                # Instantiate a new process
                p = Process(target=cdo_annualmean,args=(IN_PATH, scen+rv, OUT_PATH))
                if rc % nProc == 0 and rc != 0:
                    # Run nProc number of processes at a time
                    p.start()
                    p.join()
                    p.close() #End the process
                else:
                    p.start()
                numleft = lengthFiles - rc - 1
                ic(numleft) #Displays how many files remain to be processed in current job
