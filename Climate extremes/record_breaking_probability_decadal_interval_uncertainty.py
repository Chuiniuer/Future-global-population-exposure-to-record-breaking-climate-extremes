# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 11:55:46 2022
#record_breaking_uncertainty
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

def write_result(out_path, array):
    driver = gdal.GetDriverByName('GTiff')
    
    out_tif = driver.Create(out_path, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

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
    
def get_arr_from_tif(path):
    array = gdal.Open(path).ReadAsArray()
    return array

def write_2d_tif(extreme_index_name, output_folder, mode, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (0.125, 0.25, 0, -59.875, 0, 0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # Get geographic coordinate system
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326) # Define the output coordinate system as WGS84
    out_tif.SetProjection(srs.ExportToWkt())   # Creating projection information for new layers

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

    def get_10th_percen(self, array):
        return np.percentile(array, 10, axis=0)
    def get_median(self, array):
        return np.median(array, axis=0)

if __name__ == "__main__":
    NEX_folder = "G:\\NEX-GDDP6\\tasmax\\"

    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    
    input_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_record_breaking_proba_10years\\"
    indices = ["WD"]
    senarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    periods = range(2020, 2100, 10)
    index = IndexExtreme()
    out_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_record_prob_10years_uncertain\\"
    for i in indices:
        for s in senarios:
            for p in periods:
                for m in range(len(mode_list)):
                    file = input_folder + f"{i}_{mode_list[m]}_{p}_{s}.tif"
                    if m==0:
                        arr = np.expand_dims(get_arr_from_tif(file), axis=0)
                    else:
                        arr = np.concatenate((arr, np.expand_dims(get_arr_from_tif(file), axis=0)), axis=0)
                #median
                nodatamask = (arr==-9999).sum(axis=0)
                nodatamask[nodatamask>=1] = 1
                result = index.get_median(arr)
                result[nodatamask==1] = -9999
                write_result(out_folder + '\\' + i + "_" + str(p) + "_" + s + '_median.tif', result)
                #90th
                result = index.get_90th_percen(arr)
                result[nodatamask==1] = -9999
                write_result(out_folder + '\\' + i + "_" + str(p) + "_" + s + '_90th.tif', result)
                #10th
                result = index.get_10th_percen(arr)
                result[nodatamask==1] = -9999
                write_result(out_folder + '\\' + i + "_" + str(p) + "_" + s + '_10th.tif', result)
    
    
    
    
