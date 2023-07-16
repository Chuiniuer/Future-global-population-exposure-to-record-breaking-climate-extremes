# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 08:21:09 2022
Record-breaking probability for statistical indicators
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
# os.environ['PROJ_LIB'] = r'C:\ProgramData\Anaconda3\envs\geoplot\Lib\site-packages\pyproj\proj_dir\share\proj'
# os.environ['KMP_DUPLICATE_LIB_OK']='True'

def get_arr_from_tif(path):
    array = gdal.Open(path).ReadAsArray()
    return array

def get_arr_from_nc(path: str, var: str):
    """
    :param path: path
    :param var: variable
    :return: array
    """
    data = nc.Dataset(path)
    array = np.asarray(data.variables[var]).astype(np.float32)
    return array

def coordinate_convert_2d(array):
    array_temp = array.copy()
    array_temp[:, :720] = array[:, 720:1440]
    array_temp[:, 720: 1440] = array[:, :720]
    return array_temp

def unit_convert_pr(array):
    """
    Kg/m3/s -> mm/day * 10
    :param array:
    :return:
    """

    array[array==1e20] = -1
    array = (array * 86400)
    array[array<0] = -9999
    # array = array.astype(np.int16)
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
def write_result(extreme_index_name, output_folder, mode, period, s, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + "_" + period + "_" + s + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (-180, 0.25, 0, 90, 0, -0.25)  # ?????
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
    
    

if __name__ == "__main__":
    input_folder = "E:\\weather extreme\\SFH\\SFH\\"
    # mode_list = []
    # files = [input_folder + d for d in os.listdir(input_folder) if d[-4:] == ".tif"]
    # for s in files:
    #     s1 = s.split("\\")[-1]
    #     mode_list.append(s1.split("_")[1])
    # mode_list = list(set(mode_list))
    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    record_breaking_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_record_breaking_result\\"
    files = os.listdir(record_breaking_folder)
    out_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_record_breaking_prob\\"
    senarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    mid_range = range(2031, 2061)
    lat_range = range(2071, 2101)
    indices = ["WD"]
    
    ref_mask = gdal.Open(r"G:\weather extreme\Revision\heat(tasmax)\threshold\tasmax90p_ACCESS-CM2.tif").ReadAsArray()[::-1, :]
    
    # Traversing document statistics
    for m in tqdm(mode_list):
        for s in senarios:
            for i in indices:
                # Medium-term future
                path_list = []
                for y in mid_range:
                    if f"{i}_{m}_{s}_{y}.tif" in files:
                        path_list.append(record_breaking_folder + f"{i}_{m}_{s}_{y}.tif")
                    # Traversing files integrated together to calculate probabilities
                for j in range(len(path_list)):
                    if j == 0:
                        arr = np.expand_dims(get_arr_from_tif(path_list[j]), axis=0)
                    else:
                        arr = np.concatenate((arr, np.expand_dims(get_arr_from_tif(path_list[j]), axis=0)), axis=0)
                arr = arr.astype(np.float32)
                result = arr.sum(axis=0)/arr.shape[0]
                result[arr[0]==-9999] = -9999
                result[ref_mask==-9999] = -9999
                result = coordinate_convert_2d(result)
                write_result(i, out_folder, m, "mid", s, result)
                
                # distant future
                path_list = []
                for y in lat_range:
                    if f"{i}_{m}_{s}_{y}.tif" in files:
                        path_list.append(record_breaking_folder + f"{i}_{m}_{s}_{y}.tif")
                    #Traversing files integrated together to calculate probabilities
                for j in range(len(path_list)):
                    if j == 0:
                        arr = np.expand_dims(get_arr_from_tif(path_list[j]), axis=0)
                    else:
                        arr = np.concatenate((arr, np.expand_dims(get_arr_from_tif(path_list[j]), axis=0)), axis=0)
                arr = arr.astype(np.float32)
                result = arr.sum(axis=0)/arr.shape[0]
                result[arr[0]==-9999] = -9999
                result[ref_mask==-9999] = -9999
                result = coordinate_convert_2d(result)
                write_result(i, out_folder, m, "lat", s, result)