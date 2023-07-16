# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 21:58:12 2022

@author: Bohao Li
"""

import pandas as pd
import numpy as np
import netCDF4 as nc
from datetime import datetime
from osgeo import gdal, osr
import os
from tqdm import tqdm
import glob
import matplotlib.pyplot as plt
import xarray as xr
import time
from concurrent.futures import ProcessPoolExecutor
# os.environ['PROJ_LIB'] = r'C:\ProgramData\Anaconda3\envs\geoplot\Lib\site-packages\pyproj\proj_dir\share\proj'
# os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
"""
The HI output is scaled up by a scaler of 10 and converted to an integer.
"""


def get_arr_from_nc(path: str, var: str):
    """
    :param path: path
    :param var: variable
    :return:array
    """
    data = nc.Dataset(path)
    array = np.asarray(data.variables[var]).astype(np.float32)[:, ::-1, :]
    return array


def coordinate_convert(array):
    array_temp = array.copy()
    array_temp[:, :, :720] = array[:, :, 720:1440]
    array_temp[:, :, 720: 1440] = array[:, :, :720]
    return array_temp


def unit_convert_pr(array):
    """
    Kg/m3/s -> mm/day * 10
    :param array:
    :return:
    """

    array[array == 1e20] = -1
    array = (array * 86400)
    array[array < 0] = -9999
    # array = array.astype(np.int16)
    return array


def unit_convert_tas(array):
    """
    K -> ℃
    :param array:
    :return:
    """

    array[array == 1e20] = -10000
    array = array - 273.15
    array[array < -9999] = -9999
    # array = array.astype(np.int16)
    return array


def unit_convert_C2F(array):
    """
    ℃ ->℉
    :param array:
    :return:
    """

    array[array == 1e20] = -10000
    array = array * (9/5) + 32
    array[array < -9999] = -9999
    # array = array.astype(np.int16)
    return array


def write_2d_tif(extreme_index_name, output_folder, mode, senario, year, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
        senario + "_" + str(year) + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1,
                            gdal.GDT_Int16, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (0, 0.25, 0, 90, 0, -0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # Get geographic coordinate system
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # Define the output coordinate system as WGS84
    out_tif.SetProjection(srs.ExportToWkt())  # Creating projection information for new layers

    # writing the results
    out_tif.GetRasterBand(1).WriteArray(array)
    out_tif.GetRasterBand(1).SetNoDataValue(-9999)
    out_tif.FlushCache()
    # print(f'output successfully')
    del out_tif


def write_3d_tif(extreme_index_name, output_folder, mode, senario, year, array):
    bands = array.shape[0]
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
        senario + "_" + str(year) + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600,
                            bands, gdal.GDT_Int16, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (0, 0.25, 0, 90, 0, -0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # Get geographic coordinate system
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # Define the output coordinate system as WGS84
    out_tif.SetProjection(srs.ExportToWkt())  # Creating projection information for new layers

    for b in range(bands):

        # 数据写出
        out_tif.GetRasterBand(b+1).WriteArray(array[b])
        out_tif.GetRasterBand(b+1).SetNoDataValue(-9999)
    out_tif.FlushCache()
    # print(f'output successfully')
    del out_tif


class IndexExtreme:
    def __init__(self):
        self.days_count = 365
        self.year = 1950
        pass

    def pr_by_years(self, array):
        return array.sum(axis=0)

    def pr_daily_extreme(self, array):
        return array.max(axis=0)

    def pr_50mm_count(self, array):
        return (array >= 50).sum(axis=0)

    def get_90th_percen(self, array):
        """

        Parameters
        ----------
        array : 3d array
            DESCRIPTION.
            input 3d array with date range in 1950-2014

        Returns
        -------
        TYPE 2d np.array
            DESCRIPTION.
            90% percentile array
        """
        return np.percentile(array, 90, axis=0)


    def heat_index(self, tasmax_folder, hurs_folder, mode, senario, year):
        """
        Calculation of HI on a year-by-year basis Calculation methodology is the NWS2011 framework, reference website https://ehp.niehs.nih.gov/doi/full/10.1289/ehp.1206273
        In order to reduce storage space, the final result *10 is stored as type int16
        Parameters
        ----------
        tasmax_folder : TYPE
            DESCRIPTION.
        hurs_folder : TYPE
            DESCRIPTION.
        mode : TYPE
            DESCRIPTION.
        senario : TYPE
            DESCRIPTION.
        year : TYPE
            DESCRIPTION.
        out_folder : TYPE
            DESCRIPTION.

        Returns heat_index
        -------
        None.

        """

        forcing = os.listdir(hurs_folder + mode)[0].split("_")[-3]
        gngr = os.listdir(hurs_folder + mode)[0].split("_")[-2]
        tas = unit_convert_C2F(unit_convert_tas(get_arr_from_nc(
            tasmax_folder + f"{mode}\\{senario}\\tasmax_day_{mode}_{senario}_{forcing}_{gngr}_{year}.nc", "tasmax")))
        hurs = get_arr_from_nc(
            hurs_folder + f"{mode}\\hurs_day_{mode}_{senario}_{forcing}_{gngr}_{year}.nc", "hurs")

        result = np.full(tas.shape, -9999, dtype=np.int16)
        result[tas <= 40] = tas[tas <= 40] * 10

        A = -10.3 + (1.1 * tas) + (0.047 * hurs)
        result[(A < 79) * (tas > 40)] = A[(A < 79) * (tas > 40)] * 10

        B = -42.379 + (2.04901523 * tas) + (10.14333127 * hurs) - (0.22475541 * tas * hurs) - (6.83783 * 0.001
                                                                                               * tas * tas) - (5.481717 * 0.01 * hurs * hurs) + (1.22874 * 0.001 * tas * tas * hurs) + (8.5282 * 0.0001
                                                                                                                                                                                        * tas * hurs * hurs) - (1.99 * 0.000001 * tas * tas * hurs * hurs)

        result[(hurs <= 13) * (tas >= 80) * (tas <= 112) * (A >= 79) * (tas > 40)] = (B - (((13-hurs)/4) *
                                                                              np.sqrt((17-np.abs(tas-95))/17)))[(hurs <= 13) * (tas >= 80) * (tas <= 112) * (A >= 79) * (tas > 40)] * 10

        result[(hurs > 85) * (tas >= 80) * (tas <= 87) * (A >= 79) * (tas > 40) ] = (B + 0.02 *
                                                                       ((hurs-85) * (87-tas)))[(hurs > 85) * (tas >= 80) * (tas <= 87) * (A >= 79) * (tas > 40)] * 10

        result[(result==-9999) * (tas!=-9999)] = B[(result==-9999) * (tas!=-9999)] * 10

        result[tas == -9999] = -9999
        return result


def mission(m):

    index = IndexExtreme()
    senarios = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"]
    hurs = "E:\\hurs\\"
    tasmax = "F:\\NEX-GDDP6\\tasmax\\"
    out_folder = "F:\\HI\\"
    df = pd.DataFrame(columns=["link"])
    count=0
    for s in tqdm(senarios):
        if s == "historical":
            for y in range(1950, 2015):
                # try:
                    out_tif_name = "HI" + '_' + m + '_' + \
                        s + "_" + str(y) + '.tif'
                    if out_tif_name not in os.listdir(out_folder):
                        
                        hi = index.heat_index(tasmax, hurs, m, s, y)
                        write_3d_tif("HI", out_folder, m, s, y, hi)
                    else:
                        print(f"data {out_tif_name} already processed")
                # except:
                #     forcing = os.listdir(hurs + m)[0].split("_")[-3]
                #     gngr = os.listdir(hurs + m)[0].split("_")[-2]
                #     df.loc[count, "link"] = f"https://nex-gddp-cmip6.s3-us-west-2.amazonaws.com/NEX-GDDP-CMIP6/{m}/{s}/{forcing}/hurs/hurs_day_{m}_{s}_{forcing}_{gngr}_{y}.nc"
                #     count += 1
                #     print(f"{y} has encountered a problem..")
        else:
            for y in range(2015, 2101):
                try:
                    out_tif_name = "HI" + '_' + m + '_' + \
                        s + "_" + str(y) + '.tif'
                    if out_tif_name not in os.listdir(out_folder):
                        hi = index.heat_index(tasmax, hurs, m, s, y)
                        write_3d_tif("HI", out_folder, m, s, y, hi)
                    else:
                        print(f"data {out_tif_name} already processed")
                except:
                    forcing = os.listdir(hurs + m)[0].split("_")[-3]
                    gngr = os.listdir(hurs + m)[0].split("_")[-2]
                    df.loc[count, "link"] = f"https://nex-gddp-cmip6.s3-us-west-2.amazonaws.com/NEX-GDDP-CMIP6/{m}/{s}/{forcing}/hurs/hurs_day_{m}_{s}_{forcing}_{gngr}_{y}.nc"
                    count += 1
                    print(f"{y} has encountered a problem..")
    df.to_csv(f"E:\\hurs_sup\\{m}_link.csv")

def main():
    NEX_folder = "F:\\NEX-GDDP6\\tasmax"

    # mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
    #              "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
    #              "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
    #              "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
    #              "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
    #              "NorESM2-MM", "UKESM1-0-LL"]
    mode_list = ["ACCESS-CM2"]
    

    # with ProcessPoolExecutor(max_workers=2) as pool:
    #     pool.map(mission, mode_list)
    
    for m in mode_list:
        mission(m)

if __name__ == "__main__":

    main()
