import region_library as rlib
from icecream import ic

# boxIPCC = (rlib.AlaskaNorthwestCanada(), rlib.CentralAsia(), rlib.CanadaGreenlandIceland(), rlib.CentralNorthAmerica(), rlib.EastAfrica(), rlib.EastAsia(), rlib.EastNorthAmerica(), rlib.NorthAsia(), rlib.NorthAustralia(), rlib.NortheastBrazil(), rlib.SouthAustraliaNewZealand(), rlib.SoutheastAsia(), rlib.TibetanPlateau(), rlib.WestAsia(), rlib.WestNorthAmerica(), rlib.Antarctica(), rlib.Arctic(), rlib.PacificIslandsRegion2(), rlib.PacificIslandsRegion3(), rlib.SouthernTropicalPacific(), rlib.WestIndianOcean())
# merCrossIPCC = (rlib.SouthEuropeMediterranean(), rlib.SouthernAfrica(), rlib.Sahara(), rlib.WestAfrica())
# polyIPCC = (rlib.Amazon(), rlib.CentralAmericaMexico(), rlib.SmallIslandsRegionsCaribbean(), rlib.SouthAsia(), rlib.SoutheasternSouthAmerica(), rlib.WestCoastSouthAmerica())
polyCrossIPCC = (rlib.CentralEurope(), rlib.NorthEurope())

for reg in polyCrossIPCC:
    ic(reg)
    rlib.test_region(reg)
