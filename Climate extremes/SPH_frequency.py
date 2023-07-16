# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:47:24 2023

@author: Bohao Li
"""

import pandas as pd
import numpy as np
import netCDF4 as nc
from datetime import datetime
from osgeo import gdal, osr
import os
# from tqdm import tqdm
import glob
# import matplotlib.pyplot as plt
import xarray as xr
import time
from concurrent.futures import ProcessPoolExecutor

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


def get_sequential_flood_heatwave(flood_event_series, humid_heatwave_series):
    """
    Get the frequency of compound flood heat events
    The pane tries 4 first

    Parameters
    ----------
    flood_event_series : TYPE np.ndarray
        DESCRIPTION.
        Binary sequence of flood events, 2 is an outlier, 0 is a non-flood, 1 is a flood
    humid_heatwave_series : TYPE np.ndarray
        DESCRIPTION.
        Sequence of hot and humid heat wave events, 2 is an outlier, 0 is a non-high temperature heat wave start node, 1 is a high temperature heat wave start node

    Returns 
    -------
    Returns the frequency of the composite event, of type np.ndarray

    """
    
    if flood_event_series[0] == 2:
        return -9999
    
    freq = 0
    time_window = 6
    for t in range(len(flood_event_series)):
        if flood_event_series[t] == 1:
            for i in range(t, t+time_window+1):
                if humid_heatwave_series[t] == 1:
                    freq += 1
                    break
    return freq

def mission(m):
    
    Humid_heat_folder = "/work/home/bjsfdxlbh/Bohao/population_exposure/heatwave_output/humid_heatwave_result/"
    flood_event_folder = "/work/home/bjsfdxlbh/Bohao/population_exposure/heatwave_output/flood_event/"
    output_folder = "/work/home/bjsfdxlbh/Bohao/population_exposure/heatwave_output/SPHFQ/"
    
    scenarios = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"]
    for s in scenarios:
        
        if s == "historical":
            years = range(1950, 2015)
        else:
            years = range(2015, 2101)
        for y in years:
            
            if "SPHFreq" + '_' + m + '_' + \
                           s + "_" + str(y) + '.tif' in os.listdir(output_folder):
                               continue
            Humid_heat_array = gdal.Open(Humid_heat_folder + f"HumidHeatWave_{m}_{s}_{y}.tif").ReadAsArray()
            flood_event_array = gdal.Open(flood_event_folder + f"FloodEvent_{m}_{s}_{y}.tif").ReadAsArray()
            
            result = np.zeros((flood_event_array.shape[1], flood_event_array.shape[2]), dtype=np.int16)
            
            for i in range(flood_event_array.shape[1]):
                for j in range(flood_event_array.shape[2]):
                    result[i][j] = get_sequential_flood_heatwave(flood_event_array[:, i, j], Humid_heat_array[:, i, j])
            write_2d_tif("SPHFreq", output_folder, m, s, y, result)
            
if __name__ == "__main__":
    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    
    # for m in mode_list:
    #     mission(m)
    
    with ProcessPoolExecutor() as pool:
        pool.map(mission, mode_list)
    