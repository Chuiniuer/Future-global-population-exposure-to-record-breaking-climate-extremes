# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 12:54:28 2023

Note: the scale of the HI data is 10

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

def wirte_2d_tif(extreme_index_name, out_tif_name, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Int16, options=["COMPRESS=LZW"])

    # 设置影像的显示范围
    geotransform = (0, 0.25, 0, 90, 0, -0.25)  # ?????
    out_tif.SetGeoTransform(geotransform)

    # 获取地理坐标系
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # 定义输出的坐标系统为WGS84
    out_tif.SetProjection(srs.ExportToWkt())  # 给新建图层创建投影信息

    # 数据写出
    out_tif.GetRasterBand(1).WriteArray(array)
    out_tif.GetRasterBand(1).SetNoDataValue(-9999)
    out_tif.FlushCache()
    # print(f'output successfully')
    del out_tif
    
class IndexExtreme:
    def __init__(self):
        pass

    def HImax_high_count(self, array, th):
        
        result = array.copy()
        for s in tqdm(range(array.shape[0])):
            result[s] = array[s] >= th
        return result.sum(axis=0)    
    
def main(d):
    out_folder = "F:\\HI_heat\\"
    mode = d.split("_")[1]
    if d.split("\\")[-1] not in os.listdir(out_folder):
        
        
        array = tas_F2C(gdal.Open(d).ReadAsArray())
        
        
        
        HI_th_folder = "F:\\hi_th\\"
        
        th = gdal.Open(HI_th_folder + f"HI90p_{mode}.tif").ReadAsArray()[::-1, :]
        
        th[th<=40.6]=40.6
        index = IndexExtreme()
        result = index.HImax_high_count(array, th)
        wirte_2d_tif("HI406", out_folder + d.split("\\")[-1], result)
    
    else:
        pass
    
if __name__ == "__main__":
    
    
    
    input_folder = "F:\\HI\\"
    
    data = [input_folder + f for f in os.listdir(input_folder) if f[-4:]==".tif"]
    
    
    for d in data:
        main(d)
    
    # with ProcessPoolExecutor(max_workers=8) as pool:
    #     pool.map(main, data)

    
        
        
        
    
    
