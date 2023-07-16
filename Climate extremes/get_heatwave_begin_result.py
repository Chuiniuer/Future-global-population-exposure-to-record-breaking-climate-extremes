# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 19:19:25 2023

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


def get_arr_from_nc(path: str, var: str):
    """
    :param path: path
    :param var: variable
    :return: array
    """
    data = nc.Dataset(path)
    array = np.asarray(data.variables[var]).astype(np.float32)
    return array

def unit_convert_tas(array):
    """
    K -> ℃
    :param array:
    :return:
    """

    array[array==1e20] = -10000
    array = array - 273.15
    array[array<-9999] = -9999
    # array = array.astype(np.int16)
    return array

def tas_F2C(array):
    """
    F -> ℃

    Parameters
    ----------
    array : TYPE
        DESCRIPTION.

    Returns
    -------

    """
    array = array.atype(np.float32)
    array = array / 10
    array = (array - 32) / 1.8
    array[array<-100] = -9999
    return array

def write_2d_tif(extreme_index_name, output_folder, mode, scenario, year, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
                   scenario + "_" + str(year) + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Byte, options=["COMPRESS=LZW"])

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

def write_3d_tif(extreme_index_name, output_folder, mode, scenario, year, array):
    bands = array.shape[0]
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
        scenario + "_" + str(year) + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600,
                            bands, gdal.GDT_Byte, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (0, 0.25, 0, 90, 0, -0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # Get geographic coordinate system
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # Define the output coordinate system as WGS84
    out_tif.SetProjection(srs.ExportToWkt())  # Creating projection information for new layers

    for b in range(bands):

        # writing the results
        out_tif.GetRasterBand(b+1).WriteArray(array[b])
        out_tif.GetRasterBand(b+1).SetNoDataValue(2)
    out_tif.FlushCache()
    # print(f'output successfully')
    del out_tif

def mission(m):
    df = pd.DataFrame(columns=["scenario", "year"])
    count = 0
    input_folder = "G:\\weather extreme\\Revision\\heatwave\\heatwave_begin\\"
    output_folder = "F:\\Heatwave\\heatwave_result\\"
    log_folder = "F:\\Heatwave\\heatwavebegin_log\\"
    scenarios = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"]
    for s in tqdm(scenarios):
        if s == "historical":
            years = range(1950, 2015)
        else:
            years = range(2015, 2101)
        
        for y in years:
            try:
                if "HeatWave" + '_' + m + '_' + \
                    s + "_" + str(y) + '.tif' in os.listdir(output_folder):
                        continue
                if y == 1950 or y == 2015:
                    result = gdal.Open(input_folder + f"HeatWave_{m}_{s}_{y}.tif").ReadAsArray()
                    write_3d_tif("HeatWave", output_folder, m, s, y, result)
                else:
                    hw = gdal.Open(input_folder + f"HeatWave_{m}_{s}_{y}.tif").ReadAsArray().astype(np.float32)
                    is_begin = gdal.Open(input_folder + f"NextYearBeginHeat_{m}_{s}_{y}.tif").ReadAsArray().astype(np.float32)
                    result = hw.copy()
                    result[0, :, :] = hw[0, :, :] - is_begin
                    # result[result==-1] = 0
                    write_3d_tif("HeatWave", output_folder, m, s, y, result)
            except:
                print(f"{m}_{s}_{y} read error")
                df.loc[count, "scenario"] = s
                df.loc[count, "year"] = y
                count += 1
    df.to_csv(log_folder + f"{m}_read_error_files.csv")
def main():

    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    # mode_list = ["CanESM5", "CNRM-CM6-1"]
    # for m in mode_list:
    #     mission(m)
    with ProcessPoolExecutor(max_workers=8) as pool:
        pool.map(mission, mode_list)

if __name__ == "__main__":
    main()