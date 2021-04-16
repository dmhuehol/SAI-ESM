IN_CARD='/Users/dhueholt/Documents/GLENS_data/controlForCDO/test.nc'
MERGE_FILE='/Users/dhueholt/Documents/GLENS_data/controlForCDO/testMergeShell.nc'
SHIFT_FILE='/Users/dhueholt/Documents/GLENS_data/controlForCDO/shiftMergeShell.nc'

cdo mergetime $IN_CARD $MERGE_FILE
cdo shifttime,'-1days' $MERGE_FILE $SHIFT_FILE
