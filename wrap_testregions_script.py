''' wrap_testregions_script
Runs region testing function to plot input regions on a map.

Written by Daniel Hueholt | June 2021
Graduate Research Assistant at Colorado State University
'''

from icecream import ic
import sys

import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import process_glens_fun as pgf
import region_library as rlib

boxIPCC = (rlib.AlaskaNorthwestCanada(), rlib.CentralAsia(), rlib.CanadaGreenlandIceland(), rlib.CentralNorthAmerica(),
           rlib.EastAfrica(), rlib.EastAsia(), rlib.EastNorthAmerica(), rlib.NorthAsia(), rlib.NorthAustralia(), rlib.NortheastBrazil(),
           rlib.SouthAustraliaNewZealand(), rlib.SoutheastAsia(), rlib.TibetanPlateau(), rlib.WestAsia(), rlib.WestNorthAmerica(),
           rlib.Antarctica(), rlib.Arctic(), rlib.PacificIslandsRegion2(), rlib.PacificIslandsRegion3(), rlib.SouthernTropicalPacific(),
           rlib.WestIndianOcean())
boxMerCrossIPCC = (rlib.SouthEuropeMediterranean(), rlib.SouthernAfrica(), rlib.Sahara(), rlib.WestAfrica())
polyIPCC = (rlib.Amazon(), rlib.CentralAmericaMexico(), rlib.SmallIslandsRegionsCaribbean(), rlib.SouthAsia(),
            rlib.SoutheasternSouthAmerica(), rlib.WestCoastSouthAmerica())
polyMerCrossIPCC = (rlib.CentralEurope(),rlib.NorthEurope())

# The same 9 colors repeated a bunch of times
colors = ('viridis','magma','Purples_r', 'Greens_r','Greys_r','Oranges_r','spring','winter','cool',
          'viridis','magma','Purples_r', 'Greens_r','Greys_r','Oranges_r','spring','winter','cool',
          'viridis','magma','Purples_r', 'Greens_r','Greys_r','Oranges_r','spring','winter','cool',
          'viridis','magma','Purples_r', 'Greens_r','Greys_r','Oranges_r','spring','winter','cool',
          'viridis','magma','Purples_r', 'Greens_r','Greys_r','Oranges_r','spring','winter','cool')
CL = 0.
mapProj = cartopy.crs.EqualEarth(central_longitude = CL)
fig = plt.figure(figsize=(12, 2.73*2))
ax = plt.subplot(1, 1, 1, projection=mapProj) #nrow ncol index

for rc,reg in enumerate(boxIPCC + boxMerCrossIPCC + polyIPCC + polyMerCrossIPCC):
    ic(reg)
    fig,ax = rlib.test_region(reg,colors[rc],fig,ax)
plt.title("WG1-AR5 IPCC regions")
plt.savefig('/Users/dhueholt/Documents/GLENS_fig/20210624_finalRegionsNewQ/map_IPCCWG1AR5.png', dpi=400)
