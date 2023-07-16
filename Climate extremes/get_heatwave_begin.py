# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 18:18:39 2023

Note that HI is stored in degrees Fahrenheit
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
    array = array.astype(np.float32)
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

class IndexExtreme:
    def __init__(self):
        self.days_count = 365
        self.year = 1950
        pass


    def get_heatwave_begin(self, time_series):

        #Here next_year_begin_heat is used to characterise whether or not the end-of-year heatwave bridges to the beginning of the second year, with 1 meaning yes and 0 meaning no, and if it's 1, the first heatwave of the second year needs to be removed
        next_year_begin_heat = 0
        if time_series[0] == 2:
            return time_series[:self.days_count], 0
        # Variables for counting the number of heat wave events and for counting whether or not a heat wave counts as a heat wave
        heat_wave_begin = 0
        heat_days = 0
        count = 0
        result = np.zeros((self.days_count,), np.int8)
        # Determine if there was a heat wave at the beginning of the year, default no
        for t in range(len(time_series)):

            if t != self.days_count-1:
                if time_series[t] == 1:
                    if count == 0:
                        heat_wave_begin = t
                        count += 1
                    heat_days += 1
                else:
                    if heat_days >= 3:
                        result[heat_wave_begin] = 1

                    count = 0
                    heat_days = 0
            else:
                if self.year == 2014 or self.year == 2100:
                    result[-1] = 0
                    break
                else:
                    if time_series[t] == 0:
                        if heat_days >= 3:
                            result[heat_wave_begin] = 1
                            count = 0
                        heat_days = 0
                        break
                    else:
                        if count == 0:
                            heat_wave_begin = t
                            count += 1
                        heat_days += 1
                        for tn in range(t+1, t+3):
                            if time_series[tn] == 1:
                                heat_days += 1
                            else:
                                heat_days = 0
                        if heat_days >= 3:
                            result[heat_wave_begin] = 1

                        if time_series[t+1] + time_series[t+2] + time_series[t+3] == 3:
                            next_year_begin_heat = 1
                        else:
                            next_year_begin_heat = 0
                        break
        # Actually, here it means that an event is only kicked out if the last day is high temperature and the following three days are consecutively high temperature
        return result, next_year_begin_heat

    def get_heat_wave(self, input_folder, mode, scenario, year, out_folder, th):
        # Note that here th is a two-dimensional array
        input_folder = input_folder + mode + "\\" + scenario + "\\"
        df = pd.DataFrame(columns=["data"])
        num = 0
        # try:
        if "HeatWave" + '_' + mode + '_' + \
            scenario + "_" + str(year) + '.tif' not in os.listdir(out_folder):

            force = os.listdir(input_folder)[0].split("_")[-3]
            gngr = os.listdir(input_folder)[0].split("_")[-2]
            path = input_folder + f"tasmax_day_{mode}_{scenario}_{force}_{gngr}_{year}.nc"
            tmp_array = unit_convert_tas(get_arr_from_nc(path, "tasmax"))[:, ::-1, :]
            th1 = np.repeat(th[np.newaxis, :, :], tmp_array.shape[0], axis=0)
            array = (tmp_array >= th1).astype(np.int8)
            array[tmp_array==-9999]=2
            del tmp_array
            self.year = year
            self.days_count = array.shape[0]

            if year != 2014 and year != 2100:
                next_year = input_folder + f"tasmax_day_{mode}_{scenario}_{force}_{gngr}_{year+1}.nc"
                next_tmp_array = unit_convert_tas(get_arr_from_nc(next_year, "tasmax"))[:, ::-1, :]
                if next_tmp_array.shape[0] != array.shape[0]:
                    th1 = np.repeat(th[np.newaxis, :, :], next_tmp_array.shape[0], axis=0)
                next_array = (next_tmp_array >= th1).astype(np.int8)
                next_array[next_tmp_array==-9999] = 2
                del next_tmp_array
                array = np.concatenate((array, next_array), axis=0)

            next_year_begin_heat = np.zeros((array.shape[1], array.shape[2]), dtype=np.int8)
            temp_result = np.zeros((self.days_count, array.shape[1], array.shape[2]), dtype=np.int8)
            for i in tqdm(range(array.shape[1])):
                for j in range(array.shape[2]):
                    temp_result[:, i, j], next_year_begin_heat[i, j] = self.get_heatwave_begin(array[:, i, j])

            write_3d_tif("HeatWave", out_folder, mode, scenario, year, temp_result)
            write_2d_tif("NextYearBeginHeat", out_folder, mode, scenario, year+1, next_year_begin_heat)
        # except:
        #     df.loc[num, "data"] = f"HW_{mode}_{scenario}_{year}"
        #     num += 1
        #     df.to_csv(r"G:\HW_check.csv")
def mission(m):
    index = IndexExtreme()
    input_folder = "/work/home/bjsfdxlbh/Bohao/population exposure/tasmax/"
    out_folder = "/work/home/bjsfdxlbh/Bohao/population exposure/heatwave_output/heatwave_begin/"
    threshold_folder = "/work/home/bjsfdxlbh/Bohao/population exposure/threshold/threshold/"
    threshold_path = threshold_folder + f"tasmax90p_{m}.tif"
    threshold = gdal.Open(threshold_path).ReadAsArray()[::-1, :]
    
    threshold[threshold<=35] = 35
    
    for s in tqdm(["historical", "ssp126", "ssp245", "ssp370", "ssp585"]):

        if s=="historical":
            for y in range(1950, 2015):
                try:
                    index.get_heat_wave(input_folder, m, s, y, out_folder, threshold)
                except:
                    with open("G:\\weather extreme\\log_0131.txt", "a") as f:
                        f.write(f"{m}_{s}_{y}_has encountered a problem!\n")
        else:
            for y in range(2015, 2101):
                try:
                    index.get_heat_wave(input_folder, m, s, y, out_folder, threshold)
                except:
                    with open("G:\\weather extreme\\log_0131.txt", "a") as f:
                        f.write(f"{m}_{s}_{y}_has encountered a problem!\n")

def main():
    # 模式列表
    # mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
    #               "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
    #               "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
    #               "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
    #               "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
    #               "NorESM2-MM", "UKESM1-0-LL"]
    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1"]
    for m in mode_list:
        mission(m)
    # with ProcessPoolExecutor(max_workers=4) as pool:
    #     pool.map(mission, mode_list)

if __name__ == "__main__":
    main()
