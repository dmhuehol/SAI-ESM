''' generate_periods
Generate periods during the Last Millennium prior to 1850 and avoiding
volcanic perturbations. Volcano dates taken from the eVolv2k dataset,
which is used in the CESM2(WACCM6ma) Last Millennium simulation.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys

from icecream import ic
import numpy as np
import xarray as xr

def periods(th=10, vol_spc=5, start_yr=850, end_yr=1850):
    ''' Generate and return dictionary with periods 
    The volcano list is generated automatically,
    but the periods need hand tuning!
    th: threshold for serious volcano in Tg
    vol_spc: years to consider contaminated by volcano
    start_yr: start year
    end_yr: end year
    '''
    vol_path = "/Users/dhueholt/Documents/ecology_data/"
    vol_data = "eVolv2k_v2_ds_1.nc"
    ds_vol = xr.open_dataset(vol_path + vol_data)
    da_erupt_yrs = ds_vol.yearCE
    np_ssi = ds_vol.ssi.data
    np_ssi_gth = np_ssi > th
    da_erupt_yrs_gth = da_erupt_yrs[np_ssi_gth]
    erupt_yrs_gth_lm = da_erupt_yrs_gth[da_erupt_yrs_gth >= 850].data
    l_vol_contam = list()
    for eygl in erupt_yrs_gth_lm:
        l_vol_contam.append(np.arange(eygl, eygl + vol_spc + 1)) #+1 for inclusivity
    np_vol_contam = np.ravel(np.array(l_vol_contam))
    # ic(sorted(np_vol_contam))

    # I would prefer this to be automatic but xkcd.com/1205/ applies here for now
    clist = [ #th=10 spc=5
        [850, 869], [870, 889], [890, 909], [910, 929],
        [945, 964], [965, 984], [985, 1004], [1005, 1024],
        [1025, 1044], [1045, 1064], [1065, 1084], [1085, 1104],
        [1114, 1133], [1134, 1153], [1188, 1207], [1208, 1227], 
        [1236, 1255], [1292, 1311], [1312, 1331], [1351, 1370], 
        [1371, 1390], [1391, 1410], [1411, 1430], [1431, 1450],
        [1464, 1483], [1484, 1503], [1504, 1523], [1524, 1543], 
        [1544, 1563], [1564, 1583], [1606, 1625], [1646, 1665], 
        [1666, 1685], [1701, 1720], [1721, 1740], [1741, 1760], 
        [1761, 1780], [1789, 1808]
    ]
    # clist = [ #th=5 spc=5
    #     [850, 869], [870, 889], [945, 964], [982, 1001],
    #     [1002, 1021], [1034, 1053], [1054, 1073], [1074, 1093],
    #     [1114, 1133], [1134, 1153], [1197, 1216], [1292, 1311],
    #     [1312, 1331], [1351, 1370], [1371, 1390], [1391, 1410],
    #     [1411, 1430], [1431, 1450], [1483, 1502], [1503, 1522],
    #     [1523, 1542], [1543, 1562], [1563, 1582], [1606, 1625], 
    #     [1646, 1665], [1701, 1720], [1721, 1740], [1741, 1760],
    #     [1789, 1808]
    # ]
    # clist = [ #Randomly-chosen 20-year periods obeying th=5 spc=5 non-overlap
    #     [851, 870], [983, 1012], [1083, 1102], [1149, 1170],
    #     [1208, 1227], [1306, 1325], [1366, 1385], [1397, 1416],
    #     [1520, 1539], [1674, 1693]
    # ]

    # ind = [2, 3, 7, 13, 22, 27, 28, 30, 32, 37] #10 random inds chosen with random.org
    # ind = [0, 3, 7, 8, 9, 14, 15, 20, 21, 25] # 10 random inds chosen with random.org
    # ind = np.arange(0, 10)
    # ind = np.arange(0, len(clist)) # Supplementary Fig. 9a
    
    per_dict = {
        "eruptions": erupt_yrs_gth_lm,
        "contaminated_yrs": np_vol_contam,
        "per_all": clist,
        # Next line for all but Supplementary Fig. 9a
        # "per_ens": [clist[ei] for ei in ind] # Use for more flexible index selection
        # Below for Fig. 3a; Supplementary Fig. 9b
        "per_ens":[ # Chosen ad hoc early in development to avoid large volcanoes
            [851, 870], [871, 890], [891, 910], [911, 930], [945, 964],
            [971, 990], [991, 1010], [1011, 1030], [1031, 1050], [1051, 1070],]
            
    }

    return per_dict
