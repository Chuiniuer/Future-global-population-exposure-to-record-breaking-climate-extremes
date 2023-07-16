# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 20:42:36 2023
population exposure=probability * avg(population)
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
    input_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\WD_record_breaking_uncertain\\"
    output_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\exposure\\"
    pop_path = "G:\\weather extreme\\Population data resampled\\"
    senarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    periods = ["mid", "lat"]
    uncertainty = ["10th", "median", "90th"]
    indices=["WD"]
    
    for i in indices:
        for s in senarios:
            for p in periods:
                for u in uncertainty:
                    pop_folder = pop_path + f"{s}" + "\\"
                    file = input_folder + f"{i}_{p}_{s}_{u}.tif"
                    extreme_ds = gdal.Open(file)
                    extreme_array = extreme_ds.ReadAsArray()
                    extreme_array[extreme_array==-9999] = np.nan
                    if p == "mid":
                        date = ["2040", "2050", "2060"]
                    else:
                        date = ["2080", "2090", "2100"]
                    for d in date:
                        if s == "ssp126":
                            s_temp = "SSP1"
                        elif s == "ssp245":
                            s_temp = "SSP2"
                        elif s == "ssp370":
                            s_temp = "SSP3"
                        else:
                            s_temp = "SSP5"
                        pop_data = pop_folder + f"{s_temp}_{d}.tif"
                        if date[0] in ["2040", "2080"]:

                            arr = np.expand_dims(get_arr_from_tif(pop_data), axis=0)
                        else:
                            arr = np.concatenate((arr, np.expand_dims(get_arr_from_tif(pop_data), axis=0)), axis=0)

                    result = arr.mean(axis=0) * extreme_array
                    write_result(output_folder + f"{i}_{s}_{p}_{u}_exposure.tif", result)