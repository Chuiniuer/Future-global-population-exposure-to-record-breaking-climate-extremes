# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 19:14:58 2023

Resampling all the population data
@author: Bohao Li
"""

from osgeo import gdal
import numpy as np
import os

def mkdir(path):
    folder = os.path.exists(path)
    
    if not folder:
        os.makedirs(path)
    
    
if __name__ == "__main__":
    for s in ["SSP1", "SSP2", "SSP3", "SSP5"]:
        input_folder = f"E:\\weather extreme\\人口数据\\{s}\\Total\\GeoTIFF\\"
        output_temp = f"E:\\weather extreme\\Population data resampled temp\\"
        output_folder = f"E:\\weather extreme\\Population data resampled\\{s}\\"
        mkdir(output_folder)
        for y in range(2010, 2110, 10):
            gdal.Warp(output_temp + f"{s}_{y}.tif", input_folder + f"{s}_{y}.tif", outputBounds=(-180, -60, 180, 90), xRes=0.125, yRes=-0.125, resampleAlg="bilinear", srcNodata=2.1474836e+09, outputType=gdal.GDT_Float32, dstNodata=np.nan)
            gdal.Warp(output_folder + f"{s}_{y}.tif", output_temp + f"{s}_{y}.tif", resampleAlg="sum", xRes=0.25, yRes=0.25, outputBounds=(-180, -60, 180, 90), workingType=gdal.GDT_Float32, dstNodata=np.nan, srcNodata=np.nan)