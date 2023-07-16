# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 20:03:04 2023

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

# Getting heavy precipitation events
# The WAP input scale is 100

def write_3d_tif(extreme_index_name, output_folder, mode, senario, year, array):
    bands = array.shape[0]
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
        senario + "_" + str(year) + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600,
                            bands, gdal.GDT_Byte, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (0, 0.25, 0, 90, 0, -0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # Getting the geographic coordinate system
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


def mission(m):

    senarios = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"]
    input_folder = "E:\\WAP\\"
    out_folder = "F:\\flood_event2\\"
    threshold_folder = "E:\\WAP_threshold\\"
    log_folder = "F:\\flood_event_log_folder\\"
    WAP_threshold = gdal.Open(threshold_folder + f"WAP95p_{m}.tif").ReadAsArray()
    
    df = pd.DataFrame(columns=["scenario", "year"])
    count = 0
    for s in tqdm(senarios):
        if s == "historical":
            for y in range(1950, 2015):
                if "FloodEvent" + '_' + m + '_' + \
                    s + "_" + str(y) + '.tif' in os.listdir(out_folder):
                        continue
                try:
                    WAP = gdal.Open(input_folder + f"WAP_{m}_{s}_{y}.tif").ReadAsArray().astype(np.float32)*0.01
                    th = np.repeat(WAP_threshold[np.newaxis, :, :], WAP.shape[0], axis=0)
                    result = (WAP > th).astype(np.int8)
                    result[th==-9999] = 2
                    write_3d_tif("FloodEvent", out_folder, m, s, y, result)
                except:
                    df.loc[count, "scenarios"] = s
                    df.loc[count, "year"] = y
                    count += 1
                    print(f"{y} has encountered a problem..")
        else:
            for y in range(2015, 2101):
                try:
                    if "FloodEvent" + '_' + m + '_' + \
                        s + "_" + str(y) + '.tif' in os.listdir(out_folder):
                            continue
                    WAP = gdal.Open(input_folder + f"WAP_{m}_{s}_{y}.tif").ReadAsArray().astype(np.float32)*0.01
                    th = np.repeat(WAP_threshold[np.newaxis, :, :], WAP.shape[0], axis=0)
                    result = (WAP > th).astype(np.int8)
                    result[th==-9999] = 2
                    write_3d_tif("FloodEvent", out_folder, m, s, y, result)
                except:
                    df.loc[count, "scenarios"] = s
                    df.loc[count, "year"] = y
                    count += 1
                    print(f"{y} has encountered a problem..")
    df.to_csv(log_folder + f"{m}.csv")

if __name__ == "__main__":
    
    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    for m in mode_list:
        mission(m)        
    # with ProcessPoolExecutor() as pool:
    #     pool.map(mission, mode_list)