
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
from concurrent.futures import ProcessPoolExecutor
# os.environ['PROJ_LIB'] = r'C:\ProgramData\Anaconda3\envs\geoplot\Lib\site-packages\pyproj\proj_dir\share\proj'
# os.environ['KMP_DUPLICATE_LIB_OK']='True'
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
    K -> ℃
    :param array:
    :return:
    """

    array[array==1e20] = -1000
    array = array - 273.15
    array[array<-1000] = -9999
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
    
    #极端高温指标
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
    
def mission(m):
    NEX_folder = "G:\\HI"
    senario = "historical"
    year_range = range(1950, 2015)
    
    log = pd.DataFrame(columns=["log"])
    count = 0

    Index = IndexExtreme()
    result = np.zeros((600, 1440), dtype=np.float32)
    for c in range(8):
        for l in range(8):
            for y in year_range:
                try:
                    ds = gdal.Open(NEX_folder + "\\" +f"HI_{m}_{senario}_{y}.tif")
                    
                    if y == 1950:
                        arr = gdal.Open(NEX_folder + "\\" +\
                                             f"HI_{m}_{senario}_{y}.tif").ReadAsArray()[:, c*75:(c+1)*75, l*180:(l+1)*180]
                    else:
                        arr = np.concatenate((arr, gdal.Open(NEX_folder + "\\" + \
                                             f"HI_{m}_{senario}_{y}.tif").ReadAsArray()[:, c*75:(c+1)*75, l*180:(l+1)*180]), axis=0)
                except:
                    log.loc[count, "log"] = f"{m}_{senario}"
                    count += 1
                    print(f"{y} has encountered a problem.")
            arr = tas_F2C(arr)
            perc = Index.get_90th_percen(arr)
            result[c*75:(c+1)*75, l*180:(l+1)*180] = perc
    result = result[::-1, :]
    write_2d_tif("HI90p", "E:\\weather extreme\\hi_th", m, result)
    log.to_csv(f"F:\\hi_th_log_{m}.csv")
def main():
    #Calculation of quartiles in 8*8 chunks··
    NEX_folder = "F:\\tasmin\\"

    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    for m in mode_list:
        mission(m)        
    # with ProcessPoolExecutor(max_workers=4) as pool:
    #     pool.map(mission, mode_list)
if __name__ == "__main__":
    main()
    
        
    
    
    
    
    
    
    
    

