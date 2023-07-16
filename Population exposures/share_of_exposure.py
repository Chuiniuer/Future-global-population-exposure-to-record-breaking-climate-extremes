# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:46:03 2023

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
    srs.ImportFromEPSG(4326)  # Define the output coordinate system as WGS84
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
    input_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_record_prob_10years_uncertain\\"
    pop_output_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\share_of_pop\\"
    prob_output_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\share_of_prob\\"
    
    pop_path = "G:\\weather extreme\\Population data resampled\\"
    senarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    periods = range(2020, 2100, 10)
    uncertainty = ["10th", "median", "90th"]
    indices=["WD"]
    
    for i in indices:
        for s in senarios:
            for p in periods:
                for u in uncertainty:
                    
                    if i == "R50":
                        i_copy = "preceed50mm_pr"
                    elif i == "RX1D":
                        i_copy = "max_daily_pr"
                    elif i == "PRCPTOT":
                        i_copy = "accumulate_pr"
                    elif i == "HeatWave":
                        i_copy = "HeatWave"
                    elif i == "Tn90p":
                        i_copy = "Tn90p"
                    elif i == "Tx90p":
                        i_copy = "Tx90p"
                    elif i == "SFH":
                        i_copy = "SFH"
                    elif i == "CHDW":
                        i_copy = "CHDW"
                    else:
                        i_copy = i
                    pop_folder = pop_path + f"{s}" + "\\"
                    file = input_folder + f"{i_copy}_{p}_{s}_{u}.tif"
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
                    pop_2020 = get_arr_from_tif(pop_folder + f"{s_temp}_2020.tif")
                    nodata_mask = (pop_data==-9999)
                    pop_result = (pop_data-pop_2020) * extreme_array
                    pop_result[nodata_mask] = -9999
                    write_result(pop_output_folder + f"{i}_{s}_{p}_{u}_popShareExp.tif", pop_result)
                    
                    prob_result = pop_2020 * extreme_array
                    prob_result[nodata_mask] = -9999
                    write_result(prob_output_folder + f"{i}_{s}_{p}_{u}_probShareExp.tif", prob_result)
