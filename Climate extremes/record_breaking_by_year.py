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

os.environ['PROJ_LIB'] = r'C:\ProgramData\Anaconda3\envs\geoplot\Lib\site-packages\pyproj\proj_dir\share\proj'

def get_arr_from_tif(path):
    array = gdal.Open(path).ReadAsArray()
    return array

def write_2d_tif(extreme_index_name, output_folder, mode, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

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
    
def write_result_tif(extreme_index_name, output_folder, mode, year, senario, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + "_" + senario + "_" + str(year) + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

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
    
class IndexExtreme:
    def __init__(self):
        pass

    def pr_by_years(self, array):
        return array.sum(axis=0)

    def pr_daily_extreme(self, array):
        return array.max(axis=0)

    def pr_50mm_count(self, array):
        return (array>=50).sum(axis=0)
    
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
    
    def get_maximum_indices(self, array):
        return array.max(axis=0)
    

if __name__ == "__main__":
    input_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\EH\\"
    threshold = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_his_record_extreme_threshold\\"
    indices = ["WD"]
    senarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    mode_list = []
    files = [input_folder + d for d in os.listdir(input_folder) if d[-4:] == ".tif"]
    for s in files:
        s1 = s.split("\\")[-1]
        mode_list.append(s1.split("_")[1])
    mode_list = list(set(mode_list))
    
    out_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_record_breaking_result"
    
    for m in mode_list:
        for s in senarios:
            for i in indices:
                try:
                    for y in range(2015, 2101):
                        array = get_arr_from_tif(input_folder + f"EH_{m}_{s}_{y}.tif")
                        th = get_arr_from_tif(threshold + f"{i}_{m}.tif")
                        result = (array > th).astype(np.float32)
                        result[array==-9999] = -9999
                        write_result_tif(i, out_folder, m, y, s, result)
    
                except:
                    print(f"{m}, {s}, {i}, {y} has encountered a problem.")
    