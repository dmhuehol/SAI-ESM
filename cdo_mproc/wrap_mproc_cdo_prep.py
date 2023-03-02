'''  wrap_mproc_cdo_prep
Runs functions to process data into friendlier forms, i.e. selecting months or
calculating annual mean values. Uses multiprocessing for efficiency.

Valid IN_TOKENS: 
    '*ssp245*': UKESM no-SAI SSP2-4.5
    '*arise-sai*': UKESM-ARISE-SAI-1.5
    '*control*': GLENS no-SAI RCP8.5
    '*feedback*': GLENS SAI
    '*SSP245-TSMLT*': CESM2-ARISE-SAI-1.5
    '*BWSSP245*': CESM2(WACCM6) no-SAI SSP2-4.5
    '*BWHIST*': CESM2(WACCM6) Historical
    '*piControl*': CESM2(WACCM6) preindustrial control

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
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-ESM/cdo_mproc/mproc_bees/do_cdo_annualmean.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    # subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_annualmean.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_annualmean_pp(IN_PATH, IN_TOKEN, OUT_PATH):
    # subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_annualmean_postproc.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_annualmean_postproc.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_mergetime(IN_PATH, IN_TOKEN, OUT_PATH):
    subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-ESM/cdo_mproc/mproc_bees/do_cdo_mergetime.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    # subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_mergetime.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_seleimnsn(IN_PATH, IN_TOKEN, OUT_PATH):
    # subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_seleimnsn.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_seleimnsn.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_seleirainy(IN_PATH, IN_TOKEN, OUT_PATH):
    # subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_seleirainy.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_seleirainy.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_selfeb(IN_PATH, IN_TOKEN, OUT_PATH):
    # subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_selfeb.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_selfeb.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_sellevel(IN_PATH, IN_TOKEN, OUT_PATH):
    # subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_select_lev.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_select_lev.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_selname(IN_PATH, IN_TOKEN, OUT_PATH):
    # subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_selname.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_selname.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_selsept(IN_PATH, IN_TOKEN, OUT_PATH):
    # subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_selsept.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_selsept.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def cdo_selmons(IN_PATH, IN_TOKEN, OUT_PATH):
    # subprocess.run(['sh', '/Users/dhueholt/Documents/GitHub/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_selmons.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    subprocess.run(['sh', '/glade/u/home/dhueholt/sai-cesm/SAI-CESM/cdo_mproc/mproc_bees/do_cdo_selmons.sh', IN_PATH, IN_TOKEN, OUT_PATH])
    return None

def return_emem_list(inType):
    if inType == "raw":
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
    elif inType == 'test':
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
        ])
    elif inType == "cdo":
        EMEM=list([
        "_001_",
        "_002_",
        "_003_",
        "_004_",
        "_005_",
        "_006_",
        "_007_",
        "_008_",
        "_009_",
        "_010_",
        "_011_",
        "_012_",
        "_013_",
        "_014_",
        "_015_",
        "_016_",
        "_017_",
        "_018_",
        "_019_",
        "_020_",
        "_021_"]
        )
    elif inType == "CMIP6":
        EMEM=list([
        "r1",
        "r2",
        "r3",
        "r4",
        "r5",
        "r6",
        "r7",
        "r8",
        "r9",
        "r10"]
        )
    else:
        sys.exit('Unrecognized type!')

    identOfInt = '*' #*, or something unique if folder of data has multiple variables
    EMEM = [r + identOfInt for r in EMEM]

    return EMEM

EMEM = return_emem_list('CMIP6')
nProc = 1

# Shell inputs
# IN_PATH = '/glade/scratch/dhueholt/monthly_ICEEXTENT/'
IN_PATH = '/Users/dhueholt/Documents/ecology_data/picontrol_monthly_tas/'
IN_TOKEN = ['*piControl*',] # See docstring for valid tokens
# OUT_PATH = '/glade/scratch/dhueholt/sept_ICEEXTENT/'
OUT_PATH = '/Users/dhueholt/Documents/ecology_data/picontrol_annual_tas/'

if __name__== '__main__':
        lengthFiles = np.size(EMEM)
        for scen in IN_TOKEN:
            for rc,rv in enumerate(EMEM):
                # Instantiate a new process
                p = Process(target=cdo_annualmean, args=(IN_PATH, scen+rv, OUT_PATH))
                if rc % nProc == 0 and rc != 0:
                    # Run nProc number of processes at a time
                    p.start()
                    p.join()
                    p.close() #End the process
                else:
                    p.start()
                numleft = lengthFiles - rc - 1
                ic(numleft) #Displays how many files remain to be processed in current job
