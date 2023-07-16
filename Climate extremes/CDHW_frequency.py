# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 22:29:13 2023

@author: Bohao Li
"""
# Getting the frequency of CHDW
import pandas as pd
import numpy as np
import netCDF4 as nc 
from datetime import datetime
from osgeo import gdal, osr
import os
from tqdm import tqdm
import glob
import time
import datetime
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor
def make_dir(path):
    isExists = os.path.exists(path)
    
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

def unit_convert(array):
    """
    Kg/m3/s -> mm/day
    :param array:
    :return:
    """

    array[array==1e20] = -1
    array = (array * 86400)
    array[array<0] = -9999
    # array = array.astype(np.int16)
    return array

def write_2d_tif(extreme_index_name, output_folder, mode, scenario, year, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
                   scenario + "_" + str(year) + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Int16, options=["COMPRESS=LZW"])

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
    
def coordinate_convert(array):
    array_temp = array.copy()
    array_temp[:, :, :720] = array[:, :, 720:1440]
    array_temp[:, :, 720: 1440] = array[:, :, :720]
    return array_temp

def get_arr_from_nc(path: str, var: str):
    """
    :param path: path
    :param var: variable
    :return: array
    """
    data = nc.Dataset(path)
    array = np.asarray(data.variables[var]).astype(np.float32)
    return array

def mission(m):
    count = 0
    df = pd.DataFrame(columns=["file", "scenario", "year"])
    NEX_folder = "F:\\NEX-GDDP6\\tasmax"
    scenarios = ["ssp126", "ssp245", "ssp370", "ssp585"] 
    output_folder = "F:\\CDHWFQ\\"
    # folder = NEX_folder + "\\" + m + "\\"
    
    drought_folder = "G:\\weather extreme\\drought_events\\"
    heatwave_folder = "F:\\Heatwave\\heatwave_result\\"
    
    # data = os.listdir(folder)
    
    # Corresponding to the timing of the high temperature document based on the timing of the drought month document
    for s in scenarios:
        
        try:
            drought_events = gdal.Open(drought_folder+f"droughtEvents_{m}_{s}.tif").ReadAsArray()
        
        except:
            print(f"{m}_{s}")
            continue
        folder = NEX_folder + "\\" + m + "\\" + s + "\\"
        force = os.listdir(folder)[0].split("_")[-3]
        gngr = os.listdir(folder)[0].split("_")[-2]
        
        for y in range(1950, 2101):
            if "CDHWFreq" + '_' + m + '_' + \
                           s + "_" + str(y) + '.tif' in os.listdir(output_folder):
                               continue
            
            if y >= 2015:
                folder = NEX_folder + "\\" + m + "\\" + s + "\\"
                ref_path = folder + f"tasmax_day_{m}_{s}_{force}_{gngr}_{y}.nc"
                heatwave_begin_array = gdal.Open(heatwave_folder + f"HeatWave_{m}_{s}_{y}.tif").ReadAsArray()
            else:
                folder = NEX_folder + "\\" + m + "\\" + "historical" + "\\"
                ref_path = folder + f"tasmax_day_{m}_historical_{force}_{gngr}_{y}.nc"
                heatwave_begin_array = gdal.Open(heatwave_folder + f"HeatWave_{m}_historical_{y}.tif").ReadAsArray()
            times = nc.Dataset(ref_path).variables["time"]
            times = nc.num2date(times,times.units, times.calendar)
            
            result = np.zeros((heatwave_begin_array.shape[1], heatwave_begin_array.shape[2]), dtype=np.int16)
            
            print(len(times), heatwave_begin_array.shape, s, m, y)
            try:
                for i in tqdm(range(drought_events.shape[1])):
                    for j in range(drought_events.shape[2]):
                        
                        freq = 0
                        
                        if drought_events[3, i, j] == 2:
                            result[i, j] = -9999
                        else:
                            for k in range(12):
                                if drought_events[(y-1950)*12+k, i, j] == 1:
                                    for t in range(len(times)):
                                        month = int(str(times[t]).split("-")[1])
                                        
                                        try:
                                            if (month==(k+1)) and (heatwave_begin_array[t, i, j]==1):
                                                freq += 1
                                            elif month > (k+1):
                                                break
                                        except:
                                            if t >= 300:
                                                break
                                            else:
                                                raise ValueError
                            result[i, j] = freq
            except:
                
                print(f"{m}_{s}")
                break
            
            
            write_2d_tif("CDHWFreq", output_folder, m, s, y, result)

if __name__ == "__main__":
    mode_list = ["GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    
    # for m in mode_list:
    #     mission(m)
    
    with ProcessPoolExecutor(max_workers=8) as pool:
        pool.map(mission, mode_list)





