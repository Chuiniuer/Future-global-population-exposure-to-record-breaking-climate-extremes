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
from concurrent.futures import ProcessPoolExecutor
"""
Calculation of extreme heat days
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
    K -> ℃
    :param array:
    :return:
    """

    array[array==1e20] = -10000
    array = array - 273.15
    array[array<-9999] = -9999
    # array = array.astype(np.int16)
    return array

def write_2d_tif(extreme_index_name, output_folder, mode, senario, year, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
                   senario + "_" + year + '.tif'
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
    
def main(m):
    th_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\threshold\\"
    out_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\EH\\"
    for s in ["historical", "ssp126", "ssp245", "ssp370", "ssp585"]:
        folder = NEX_folder + "\\" + m + "\\" + s + "\\"
        files = [folder + d for d in os.listdir(folder) if d[-3:]=='.nc']
        
        th = th_folder + f"tasmax90p_{m}.tif"
        th = gdal.Open(th).ReadAsArray()[::-1, :]
        
        #new threshold
        th[th<=35] = 35
        
        for f in files:
            
            y = f.split('_')[-1].split('.')[0]
            if "EH" + '_' + m + '_' + \
                        s + "_" + y + '.tif' in os.listdir(out_folder):
                            print("already existed")
                            continue

            array = unit_convert_tas(gdal.Open(f).ReadAsArray())
            plt.imshow(array[0])
            temp_arr = array.copy()
            for ax in range(array.shape[0]):
                temp_arr[ax] = (array[ax] > th).astype(np.float32)
            result = temp_arr.sum(axis=0)
            result[th==-9999] = -9999
            
            write_2d_tif("EH", out_folder, m, s, y, result)
    
if __name__ == "__main__":

    NEX_folder = "F:\\NEX-GDDP6\\tasmax"

    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CNRM-ESM2-1",
                 "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                 "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                 "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                 "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                 "NorESM2-MM", "UKESM1-0-LL"]
    for m in mode_list:
        main(m)
    # with ProcessPoolExecutor(max_workers=8) as pool:
    #     pool.map(main, mode_list)

    
    
    

        
                
                        
            
    
    
    
    