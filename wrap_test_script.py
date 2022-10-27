import cmasher
import cmocean
import fun_plot_tools as fpt
import fun_special_plot as fsp
import matplotlib as mpl
import seaborn

# muteThresh = 13
# muteList = fpt.mute_by_numbers(muteThresh)
# muteList = fpt.mute_by_numbers_arise(muteThresh)
# disc = mpl.colors.ListedColormap(muteList)
precipPal = seaborn.diverging_palette(58, 162, s=100, l=45, as_cmap=True)

cbarDict = {
    "cmap": precipPal,
    "range": [-25,25],
    "direction": 'horizontal',
    "label": ''
}
savePath = '/Users/dhueholt/Documents/GLENS_fig/20221020_itcz/'
saveName = 'precipColorbar'

fsp.save_colorbar(cbarDict, savePath, saveName, dpiVal=400)
