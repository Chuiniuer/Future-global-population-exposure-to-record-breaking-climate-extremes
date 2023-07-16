# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 16:39:16 2023

@author: HP
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
import xarray as xr
import time
from concurrent.futures import ProcessPoolExecutor


def get_arr_from_nc(path, var):
    data = nc.Dataset(path)
    array = np.asarray(data.variables[var]).astype(np.float32)
    return array

def write_3d_tif(extreme_index_name, output_folder, mode, senario, array):
    bands = array.shape[0]
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
        senario + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600,
                            bands, gdal.GDT_Byte, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (0, 0.25, 0, 90, 0, -0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # Get geographic coordinate system
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

if __name__ == "__main__":
    input_folder = "F:\\SPEI calculated\\"
    result = get_arr_from_nc(input_folder + "SPEI_EC-Earth3_ssp585_3.nc", "SPEI")
    nodata_mask = np.isnan(result)
    result[result>-1] = 0
    result[result!=0] = 1
    result[nodata_mask] = 2
    write_3d_tif("droughtEvents", "F:\\Drought_event_comp\\", 'EC-Earth3', "ssp585", result)
    
    result = get_arr_from_nc(input_folder + "SPEI_MIROC6_ssp126_3.nc", "SPEI")
    nodata_mask = np.isnan(result)
    result[result>-1] = 0
    result[result!=0] = 1
    result[nodata_mask] = 2
    write_3d_tif("droughtEvents", "F:\\Drought_event_comp\\", 'MIROC6', "ssp126", result)
    
    