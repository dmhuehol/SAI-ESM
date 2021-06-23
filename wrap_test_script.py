import region_library as rlib
from icecream import ic

# boxIPCC = (rlib.AlaskaNorthwestCanada(), rlib.CentralAsia(), rlib.CanadaGreenlandIceland(), rlib.EastAfrica(), rlib.EastAsia(), rlib.EastNorthAmerica(), rlib.NorthAsia(), rlib.NorthAustralia(), rlib.NortheastBrazil(), rlib.SouthAustraliaNewZealand(), rlib.SoutheastAsia(), rlib.TibetanPlateau(), rlib.WestAsia(), rlib.WestNorthAmerica(), rlib.Antarctica(), rlib.Arctic(), rlib.PacificIslandsRegion2(), rlib.PacificIslandsRegion3(), rlib.SouthernTropicalPacific(), rlib.WestIndianOcean())
merCrossIPCC = (rlib.SouthEuropeMediterranean(), rlib.SouthernAfrica(), rlib.Sahara(), rlib.WestAfrica())
for reg in merCrossIPCC:
    ic(reg)
    rlib.test_region(reg)
