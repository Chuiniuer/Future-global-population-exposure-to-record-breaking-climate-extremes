# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:49:55 2022

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
import time
import datetime
from concurrent.futures import ProcessPoolExecutor


"""
Find the cumulative monthly rainfall and average temperature and output the data for SPEI calculation.
"""
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

def write_2d_tif(extreme_index_name, output_folder, mode, senario, year, month, array):
    
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
                   senario + "_" + str(year) + "_" + str(month) + '.tif'
        

    driver = gdal.GetDriverByName('GTiff')
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

    # Set the display range of the image
    geotransform = (0, 0.25, 0, 90, 0, -0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # Get geographic coordinate system
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # Define the output coordinate system as WGS84
    out_tif.SetProjection(srs.ExportToWkt())  # Creating projection information for new layers

    # Writing the results
    out_tif.GetRasterBand(1).WriteArray(array)
    out_tif.GetRasterBand(1).SetNoDataValue(np.nan)
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

def main():
   
    # 模式列表
    # mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
    #               "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
    #               "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
    #               "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
    #               "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
    #               "NorESM2-MM", "UKESM1-0-LL"]
    mode_list = [
                  "MIROC6"]
    for m in mode_list:
        mission(m)
    # with ProcessPoolExecutor(max_workers=4) as pool:
    #     pool.map(mission, mode_list)
def mission(m):
    NEX_folder = "F:\\CMIP6_GDDP\\"
    scenarios = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"] 
    output_folder = "G:\\weather extreme\compound spei heatwave"
    make_dir(output_folder + "\\" + m)
    output_folder = output_folder + "\\" + m
    folder = NEX_folder + "\\" + m + "\\"
    # data = os.listdir(folder)
    data = ["G:\\weather extreme\\compound spei heatwave\\pr_day_MIROC6_ssp126_r1i1p1f1_gn_2085.nc"]
    for d in data:
        s = d.split("_")[3]
        #Determine if a file has been generated
        # year = d.split("_")[-1].split(".")[0]
        # count = 0
        # for mon in range(1, 12+1):
        #     if mon < 10:
        #         if "monthly_tas" + '_' + m + '_' + s+ "_" + str(year) + "_" + "0" + str(mon) + '.tif' not in os.listdir(output_folder):
        #             count += 1
        #             break
        #         else:
        #             continue
        #     else:
        #         if "monthly_tas" + '_' + m + '_' + s+ "_" + str(year) + "_" + str(mon) + '.tif' not in os.listdir(output_folder):
        #             count += 1
        #             break
        #         else:
        #             continue
        # if count == 0:
        #     continue
        # print(count)
        # file = folder + "\\" + d 
        file = d
        times = nc.Dataset(file).variables["time"]
        times = nc.num2date(times,times.units, times.calendar)
        nodata = float(nc.Dataset(file).variables["pr"]._FillValue)
        # we calculate the mean temperature
        pr = get_arr_from_nc(file, "pr")
        pr[pr==nodata] = np.nan
        for t in range(len(times)):
            year = str(times[t]).split("-")[0]
            month = str(times[t]).split("-")[1]
            if t == 0:
                t1 = t
                m1 = month
                continue
            elif m1 == "12":
                out_arr = pr[t-1:, ::-1, :].sum(axis=0)
                
                write_2d_tif("monthly_pr", output_folder, m, s, year, m1, out_arr)
                t1 = t
                m1 = month
            else:
                if m1 == month:
                    continue
                else:
                    out_arr = pr[t1:t, ::-1, :].sum(axis=0)
                    write_2d_tif("monthly_pr", output_folder, m, s, year, m1, out_arr)
                    t1 = t
                    m1 = month
                        
if __name__ == "__main__":
    main()
    
                            
                            
                            
                        
                    
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        