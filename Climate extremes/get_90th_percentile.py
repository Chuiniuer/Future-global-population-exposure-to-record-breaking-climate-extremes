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

os.environ['KMP_DUPLICATE_LIB_OK']='True'
"""

"""

def get_arr_from_nc(path: str, var: str):
    """
    :param path: path
    :param var: variable
    :return: array
    """
    data = nc.Dataset(path)
    array = np.asarray(data.variables[var]).astype(np.float32)
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

    array[array==1e20] = -1
    array = (array * 86400)
    array[array<0] = -9999
    # array = array.astype(np.int16)
    return array

def unit_convert_tas(array):
    """
    K -> ℃ scale: 100
    :param array:
    :return:
    """

    array[array==1e20] = -1
    array = (array - 273.15)*100
    array[array<0] = -1
    # array = array.astype(np.int16)
    return array.astype(np.int16)

def write_2d_tif(extreme_index_name, output_folder, mode, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (0, 0.25, 0, -60, 0, 0.25)  # ?????
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
    #Calculation of quartiles in 8*8 chunks -·

    #Traversing different models, different historical years in different scenarios
    NEX_folder = "G:\\NEX-GDDP6\\tasmax"

    mode_list = os.listdir(NEX_folder)
    
    senario = "historical"
    year_range = range(1950, 2015)
    
    Index = IndexExtreme()
    
    for m in tqdm(mode_list[3:]):
        gn_gr = os.listdir(NEX_folder + "\\" + m + "\\" + \
                               senario)[0].split("_")[5]
        stress = os.listdir(NEX_folder + "\\" + m + "\\" + \
                               senario)[0].split("_")[4]
        for y in year_range:
            if y == 1950:
                arr = unit_convert_tas(get_arr_from_nc(NEX_folder + "\\" + m + "\\" + \
                                       senario + "\\" + f"tasmax_day_{m}_{senario}_{stress}_{gn_gr}_{y}.nc", "tasmax"))
            else:
                arr = np.concatenate((arr, unit_convert_tas(get_arr_from_nc(NEX_folder + "\\" + m + "\\" + \
                                   senario + "\\" + f"tasmax_day_{m}_{senario}_{stress}_{gn_gr}_{y}.nc", "tasmax"))), axis=0)
                perc = Index.get_90th_percen(arr).astype(np.float32)/100
        write_2d_tif("tasmax90p", "E:\\weather extreme\\tas\\threshold", m, perc)