# Empty s cript reserved for testing code
import cmasher
import matplotlib as mpl

import fun_special_plot as fsp
import fun_plot_tools as fpt

rbd = {
    "nRlz": 21,
    "muteThr": 15
}
discColors = cmasher.take_cmap_colors(cmasher.cm.ghostlight,
                                      N=rbd["nRlz"], cmap_range=(0,1))
muteList = fpt.mute_by_numbers(rbd["muteThr"])
disc = mpl.colors.ListedColormap(muteList)
cbarDict = {
    "cmap": disc,
    "range": [1,21],
    "direction": 'horizontal',
    "label": 'percent'
}
savePath = '/Users/dhueholt/Documents/GLENS_fig/20230321_checkThresh/'
saveName = 'glens121'
fsp.save_colorbar(cbarDict, savePath, saveName, dpiVal=400)