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

def write_2d_tif(extreme_index_name, output_folder, mode, senario, year, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
                   senario + "_" + str(year) + '.tif'
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

def write_3d_tif(extreme_index_name, output_folder, mode, senario, year, array):
    bands = array.shape[0]
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
        senario + "_" + str(year) + '.tif'
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

        # 数据写出
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
    # 降雨指标

    def get_humid_heatwave_begin(self, time_series):
        
        # Here next_year_begin_heat is used to characterise whether or not the end-of-year heatwave bridges to the beginning of the second year, with 1 meaning yes and 0 meaning no, and if it's 1, the first heatwave of the second year needs to be removed
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
        return result, next_year_begin_heat
    
    def get_heat_wave(self, input_folder, mode, senario, year, out_folder, th=40.6):
        #Note that the HI scaler here is 10
        df = pd.DataFrame(columns=["data"])
        num = 0
        try:
            if "HumidHeatWave" + '_' + mode + '_' + \
                senario + "_" + str(year) + '.tif' not in os.listdir(out_folder):
            
                path = input_folder + f"HI_{mode}_{senario}_{year}.tif"
                tmp_array = tas_F2C(gdal.Open(path).ReadAsArray())
                array = tmp_array.copy()
                for c in range(array.shape[0]):
                    array[c] = (array[c] >= th).astype(np.int8)
                # array = (tmp_array >= th).astype(np.int8)
                array[tmp_array==-9999]=2
                del tmp_array
                self.year = year
                self.days_count = array.shape[0]
                
                if year != 2014 and year != 2100:
                    next_year = input_folder + f"HI_{mode}_{senario}_{year+1}.tif"
                    next_tmp_array = tas_F2C(gdal.Open(next_year).ReadAsArray())
                    next_array = next_tmp_array.copy()
                    for c in range(next_array.shape[0]):
                        next_array[c] = (next_array[c] >= th).astype(np.int8)
                    # next_array = (next_tmp_array >= th).astype(np.int8)
                    next_array[next_tmp_array==-9999] = 2
                    del next_tmp_array
                    array = np.concatenate((array, next_array), axis=0)
                    next
                
                next_year_begin_heat = np.zeros((array.shape[1], array.shape[2]), dtype=np.int8)
                temp_result = np.zeros((self.days_count, array.shape[1], array.shape[2]), dtype=np.int8)
                for i in tqdm(range(array.shape[1])):
                    for j in range(array.shape[2]):
                        temp_result[:, i, j], next_year_begin_heat[i, j] = self.get_humid_heatwave_begin(array[:, i, j])
                
                write_3d_tif("HumidHeatWave", out_folder, mode, senario, year, temp_result)
                write_2d_tif("NextYearBeginHeat", out_folder, mode, senario, year+1, next_year_begin_heat)
        except:
            df.loc[num, "data"] = f"HI_{mode}_{senario}_{year}"
            num += 1
            df.to_csv(r"F:\HI_check.csv")
def mission(m):
    index = IndexExtreme()
    input_folder = "F:\\HI\\"
    th_folder = "F:\\hi_th\\"
    out_folder = "G:\\weather extreme\\Revision\\heatwave\\humid_heatwave_begin\\"
    for s in tqdm(["historical", "ssp126", "ssp245", "ssp370", "ssp585"]):
        
        th = th_folder + f"HI90p_{m}.tif"
        th = gdal.Open(th).ReadAsArray()[::-1, :]
        th[th<=40.6] = 40.6
        if s=="historical":
            for y in range(1950, 2015):
                index.get_heat_wave(input_folder, m, s, y, out_folder, th)
        else:
            for y in range(2015, 2101):
                index.get_heat_wave(input_folder, m, s, y, out_folder, th)

def main():
    # 模式列表
    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8"]
    # mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
    #               "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
    #               "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
    #               "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
    #               "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
    #               "NorESM2-MM", "UKESM1-0-LL"]
    # for m in mode_list:
    #     mission(m)
    with ProcessPoolExecutor(max_workers=4) as pool:
        pool.map(mission, mode_list)

if __name__ == "__main__":
    main()
