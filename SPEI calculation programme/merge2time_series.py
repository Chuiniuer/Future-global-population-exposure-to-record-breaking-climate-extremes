# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 21:42:22 2022

@author: 当年令狐冲华山论剑屎壳郎脚下一个粪球
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
from concurrent.futures import ProcessPoolExecutor
os.environ['PROJ_LIB'] = r'C:\ProgramData\Anaconda3\envs\geoplot\Lib\site-packages\pyproj\proj_dir\share\proj'
os.environ['KMP_DUPLICATE_LIB_OK']='True'

def get_arr_from_tiff(path: str):
    """
    :param path: path
    :return:array
    """
    ds = gdal.Open(path)
    array = ds.ReadAsArray()
    return array

def unit_convert_tas(array):
    """
    K -> ℃
    :param array:
    :return:
    """

    # array[array==1e20] = -10000
    array = array - 273.15
    # array[array<-9999] = -9999
    # array = array.astype(np.int16)
    return array

def write_nc(path: str, arr: np.ndarray, main_var: str, unit: str, nodata: float):
    newfile = nc.Dataset(path, 'w', format='NETCDF4')

    #define dimensions
    lon = newfile.createDimension("longitude", size=1440)
    lat = newfile.createDimension("latitude", size=600)
    times = newfile.createDimension("time", size=None)

    # define variables for storing data
    lon = newfile.createVariable("lon", np.float32, dimensions="longitude")
    lat = newfile.createVariable("lat", np.float32, dimensions="latitude")
    time = newfile.createVariable("times", "S19", dimensions="time")
    crucial_var = newfile.createVariable(main_var, np.float32, dimensions=("time", "latitude", "longitude"))
    date_range = pd.date_range(datetime(1950, 1, 15), datetime(2100, 12, 31), freq="1m")

    # add data to variables
    lon[:] = np.arange(0.125, 360.125, 0.25)
    lat[:] = np.arange(89.875, -60.125, -0.25)
    crucial_var[:, :, :] = arr
    for i in range(arr.shape[0]):
        time[i] = date_range[i].strftime("%Y-%m-%d")
    print(date_range[0].strftime("%Y-%m-%d"))
    print(time)
    # add attributes
    # add global attributes
    newfile.title = f"daily {main_var} data"
    newfile.start_time = time[i]
    newfile.times = time.shape[0]
    newfile.history = "Created" + datetime(2022, 12, 18).strftime("%Y-%m-%d")

    # add local attributes to variable
    lon.description = "longitude, range from 0 to 360"
    lon.units = "degrees"

    lat.description = "latitude, south is negative"
    lat.units = "degrees north"

    time.description = "time, unlimited dimension"
    time.units = "time since {0:s}".format(time[0])

    crucial_var.description = f"We've used original unit here without scaler."
    crucial_var.units = unit
    crucial_var.missing_value = nodata
    # close file
    newfile.close()

def write_2d_tif(extreme_index_name, output_folder, mode, senario, year, month, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
                   senario + "_" + str(year) + "_" + str(month) + '.tif'
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
   

    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    # mode_list = ["UKESM1-0-LL"]
    
    # for m in mode_list:
    #     merge2time_series(m)
    with ProcessPoolExecutor(max_workers=2) as pool:
        pool.map(merge2time_series, mode_list)
        
def merge2time_series(m):
    scenarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    input_folder = "G:\\monthly_tas\\" + m + '\\'
    out_folder = "G:\\monthly_tas_aggregated\\"
    for s in tqdm(scenarios):
        # Iterate over the years to be consolidated into a new array.
        for y in range(1950, 2100+1):
            for m1 in range(1, 12+1):
                # Determination of documents by date
                if y < 2015:
                    if m1 < 10:
                        path = input_folder + f"monthly_tas_{m}_historical_{y}_0{m1}.tif"
                    else:
                        path = input_folder + f"monthly_tas_{m}_historical_{y}_{m1}.tif"
                else:
                    if m1 < 10:
                        path = input_folder + f"monthly_tas_{m}_{s}_{y}_0{m1}.tif"
                    else:
                        path = input_folder + f"monthly_tas_{m}_{s}_{y}_{m1}.tif"
                if y == 1950 and m1 == 1:
                    arr_aggregated = np.expand_dims(unit_convert_tas(get_arr_from_tiff(path)), axis=0)
                else:
                    arr_aggregated = np.concatenate((arr_aggregated, np.expand_dims(unit_convert_tas(get_arr_from_tiff(path)), axis=0)), axis=0)
        arr_aggregated[np.isnan(arr_aggregated)] = 1e20
        print(arr_aggregated.shape)
        write_nc(out_folder+f"monthly_tas_{m}_{s}.nc", arr_aggregated, "tas", "Celsius", 1e20)
if __name__ == "__main__":
    main()