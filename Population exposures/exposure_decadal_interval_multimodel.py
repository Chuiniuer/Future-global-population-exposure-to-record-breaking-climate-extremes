# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 11:10:01 2023

@author: Bohao Li
"""

from osgeo import gdal, osr
import numpy as np
import os
# os.environ['PROJ_LIB'] = r'C:\ProgramData\Anaconda3\envs\geoplot\Lib\site-packages\pyproj\proj_dir\share\proj'

def write_result(out_path, array):
    driver = gdal.GetDriverByName('GTiff')
    
    out_tif = driver.Create(out_path, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

    # Setting the image display range
    geotransform = (-180, 0.25, 0, 90, 0, -0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # Getting the geographic coordinate system
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # Getting the geographic coordinate system
    out_tif.SetProjection(srs.ExportToWkt())  # Creating projection information for new layers

    # writing the results
    out_tif.GetRasterBand(1).WriteArray(array)
    out_tif.GetRasterBand(1).SetNoDataValue(np.nan)
    out_tif.FlushCache()
    # print(f'output successfully')
    del out_tif

def get_arr_from_tif(path):
    array = gdal.Open(path).ReadAsArray()
    return array

if __name__ == "__main__":
    input_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_record_breaking_proba_10years\\"
    output_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\exposure_span10\\multi-model result\\"
    pop_path = "G:\\weather extreme\\Population data resampled\\"
    senarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    periods = range(2020, 2100, 10)
    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    indices=["WD"]
    
    for i in indices:
        for s in senarios:
            for p in periods:
                for u in mode_list:
                    pop_folder = pop_path + f"{s}" + "\\"
                    file = input_folder + f"{i}_{u}_{p}_{s}.tif"
                    extreme_ds = gdal.Open(file)
                    extreme_array = extreme_ds.ReadAsArray()
                    extreme_array[extreme_array==-9999] = np.nan

                    if s == "ssp126":
                        s_temp = "SSP1"
                    elif s == "ssp245":
                        s_temp = "SSP2"
                    elif s == "ssp370":
                        s_temp = "SSP3"
                    else:
                        s_temp = "SSP5"
                    pop_data = get_arr_from_tif(pop_folder + f"{s_temp}_{p}.tif")
                    result = pop_data * extreme_array
                    write_result(output_folder + f"{i}_{s}_{p}_{u}_exposure.tif", result)