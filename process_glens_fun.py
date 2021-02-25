import xarray as xr
import numpy as np
import sys

# baselineFname = baselineStr[len(baselinePath):len(baselineStr)] #only filename
def extract_meta_from_fname(fname):
# Extract useful metadata from GLENS output filename
    # baselineFname = fname[len(fname):len(fname)] #only filename
    fnamePcs = fname.split('.') #period is delimiter in filename
    glensMeta = {"unkn1": fnamePcs[0],
        "unkn2": fnamePcs[1],
        "unkn3": fnamePcs[2],
        "unkn4": fnamePcs[3],
        "runType": fnamePcs[4],
        "ensMemNum": int(fnamePcs[5]),
        "unkn5": fnamePcs[6],
        "unkn6": fnamePcs[7],
        "dataVar": fnamePcs[8],
        "rawDates": fnamePcs[9],
        "strtYr": int(fnamePcs[9][0:4]),
        "strtMnth": int(fnamePcs[9][4:6]),
        "endYr": int(fnamePcs[9][7:11]),
        "endMnth": int(fnamePcs[9][11:13]),
        "extension": fnamePcs[10]
    }

    return glensMeta

def discover_data_var(glensDsetCntrl):
    fileKeys = list(glensDsetCntrl.keys())
    notDataKeys = ['time_bnds', 'date', 'datesec']
    dataKey = "empty"
    for cKey in fileKeys:
        if cKey in notDataKeys:
            print('not this one: ' + cKey)
        else:
            dataKey = cKey
            print('Data key is: ' + dataKey)
            if dataKey == "empty":
                sys.exit("Data discovery failed! Ending run now...")

    return dataKey

def find_matching_year_bounds(glensCntrl,glensFdbck):
    cntrlYrs = glensCntrl['time'].dt.year.data
    fdbckYrs = glensFdbck['time'].dt.year.data
    bothYrs,cntrlInd,fdbckInd = np.intersect1d(cntrlYrs, fdbckYrs, return_indices=True)
    firstYrInBoth = bothYrs[0]
    lastYrInBoth = bothYrs[len(bothYrs)-1]
    bndDct = {
        "strtYrMtch": firstYrInBoth,
        "endYrMtch": lastYrInBoth,
        "cntrlStrtMtch": cntrlInd[0],
        "fdbckStrtMtch": fdbckInd[0],
        "cntrlEndMtch": cntrlInd[len(cntrlInd)-1],
        "fdbckEndMtch": fdbckInd[len(fdbckInd)-1],
        "mtchYrs": bothYrs}

    return bndDct
