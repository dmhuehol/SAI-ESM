import region_library as rlib

boxIPCC = (rlib.AlaskaNorthwestCanada(), rlib.CentralAsia(), rlib.CanadaGreenlandIceland(), rlib.EastAfrica(), rlib.EastAsia(), rlib.EastNorthAmerica(), rlib.NorthAsia(), rlib.NorthAustralia(), rlib.NortheastBrazil(), rlib.SouthAustraliaNewZealand(), rlib.SoutheastAsia(), rlib.TibetanPlateau(), rlib.WestAsia(), rlib.WestNorthAmerica(), rlib.Antarctica(), rlib.Arctic(), rlib.PacificIslandsRegion2(), rlib.PacificIslandsRegion3(), rlib.SouthernTropicalPacific(), rlib.WestIndianOcean())

for reg in boxIPCC:
    rlib.test_region(reg)
