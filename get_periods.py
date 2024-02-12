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

def periods(th=10, vol_spc=5, span=20, spcr=20, start_yr=850, end_yr=1850):
    ''' Generate and return dictionary with periods 
    th: threshold for serious volcano in Tg
    vol_spc: years to consider contaminated by volcano
    span: length of output time period
    spacer: time between periods; spcr > span for non-overlapping
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
    
    clist = list()
    for c in np.arange(start_yr, end_yr, spcr):
        per = np.arange(c, c+span, 1)
        per_bounds = np.array([c, c + span - 1]) # -1 for inclusivity
        check_vol = np.intersect1d(np_vol_contam, per)
        if check_vol.size > 0:
            continue # Intersecting members exist
        else:
            clist.append(per_bounds)
    
    rnd = [3,4,7,8,9,16,17,20,21,22] # 10 random inds chosen with random.org (compare to ARISE)
    per_dict = {
        "eruptions": erupt_yrs_gth_lm,
        "contaminated_yrs": np_vol_contam,
        "per_all": clist,
        "per_ens": [clist[ri] for ri in rnd]
            
    }
    
    return per_dict

# def periods():
#     per_dict = {
#         "volcanoes": [
#             1835, 1831, 1822, 1815, 1783, 1766, 1755, 1739, 1721, 1707, 1673,
#             1667, 1640, 1600, 1595, 1585, 1510, 1477, 1257, 946, 939],
#         "discontinuities": [
#             '1110', '1170', '1230', '1460', '1695',
#         ],
#         "per_all_5yr_gth10": None,
        
#         "per_ens": [
#             # [851, 870], [871, 890], [891, 910], [911, 930], [957, 976],
#             # [977, 996], [997, 1016], [1017, 1036], [1037, 1056], [1051, 1070],],
#             [851, 870], [871, 890], [891, 910], [911, 930], [951, 970],
#             [971, 990], [991, 1010], [1011, 1030], [1031, 1050], [1051, 1070],],
#         # "per_all_10yr": [
#         #     [851, 870], [871, 890], [891, 910], [911, 930], [957, 976],
#         #     [977, 996], [997, 1016], [1017, 1036], [1037, 1056], [1051, 1070],
#         #     [1071, 1090], [1127, 1146], [1131, 1150], [1181, 1200], [1201, 1220],
#         #     [1268, 1287], [1288, 1307], [1308, 1327], [1328, 1347], [1348, 1367],
#         #     [1368, 1387], [1388, 1407], [1408, 1427], [1428, 1447],
#         #     [1488, 1507], [1521, 1540], [1541, 1560], [1561, 1580], [1611, 1630],
#         # ]
#         "per_all_10yr": [
#             [851, 870], [871, 890], [891, 910], [911, 930], [951, 970],
#             [971, 990], [991, 1010], [1011, 1030], [1031, 1050], [1051, 1070],
#             [1065, 1084], [1122, 1141], [1125, 1144], [1175, 1194], [1195, 1214],
#             [1263, 1282], [1283, 1302], [1303, 1322], [1323, 1342], [1343, 1362],
#             [1363, 1382], [1383, 1402], [1403, 1422], [1423, 1442],
#             [1483, 1502], [1516, 1535], [1536, 1555], [1556, 1575], [1606, 1625],
#         ]
#     }
#     return per_dict
