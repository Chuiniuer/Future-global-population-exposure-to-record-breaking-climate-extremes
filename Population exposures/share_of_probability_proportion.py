# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 18:57:37 2023

@author: Bohao Li
"""

from osgeo import gdal
import xarray as xr          
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import statsmodels.api as sm
from scipy import stats
from tqdm import tqdm
import os

def saving_new_img(out_path, array, gt, proj, xsize, ysize, nodata):
    """
    :param out_path: path
    :param array: Array of output files
    :param gt: affine transform function (math.)
    :param proj: coordinate reference system
    :param xsize: number of pixels in x axis
    :param ysize: number of pixels in y axis
    :return:
    """
    #重新生成新的文件
    driver = gdal.GetDriverByName("GTiff")
    driver.Register()
    out_ds = driver.Create(out_path, xsize=xsize, ysize=ysize, bands=1, eType=gdal.GDT_Float32, options=["COMPRESS=LZW"])
    out_ds.SetGeoTransform(gt)
    out_ds.SetProjection(proj)
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(array)
    out_band.SetNoDataValue(nodata)
    out_band.FlushCache()

    out_band = None
    out_ds = None
    

if __name__ == "__main__":
    # indices = ["PRCPTOT", "RX1D", "R50", "Tx90p", "HW", "Tn90p", "SFH", "CHDW"]
    indices = ["WD"]
    scenarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    uncertainties = ["median"]
    
    prob_trend_folder = 'G:\\weather extreme\\Revision\\heat(tasmax)\\share_of_prob_trend\\slope\\'
    total_trend_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\exposure_trend\\slope\\"
    total_pvalue_folder  = "G:\\weather extreme\\Revision\\heat(tasmax)\\exposure_trend\\significance\\"
    
    out_folder = "G:\\weather extreme\\Revision\\heat(tasmax)\\share_prob_proportion\\"
    
    for i in indices:
        for s in scenarios:
            for u in uncertainties:
                
                if i == "HW":
                    i1 = "HeatWave"
                else:
                    i1 = i
                
                prob = gdal.Open(prob_trend_folder + f"{i1}_{s}_{u}_slope.tif").ReadAsArray()
                total_trend = gdal.Open(total_trend_folder + f"{i}_{s}_{u}_slope.tif").ReadAsArray()
                pvalue = gdal.Open(total_pvalue_folder + f"{i}_{s}_{u}_pvalue.tif").ReadAsArray()
                
                prob[prob==-9999] = np.nan
                total_trend[total_trend==-9999] = np.nan
                
                total_trend[pvalue > 0.05] = np.nan
                
                total_trend[total_trend <= 0] = np.nan
                
                prob[pvalue > 0.05] = np.nan
                
                
                result = prob / total_trend
                
                result[result<0] = 0
                result[result>1] = 1
                
                gt = gdal.Open(prob_trend_folder + f"{i1}_{s}_{u}_slope.tif").GetGeoTransform()
                proj = gdal.Open(prob_trend_folder + f"{i1}_{s}_{u}_slope.tif").GetProjection()
                
                saving_new_img(out_folder + f"{i}_{s}_{u}_sharePro.tif", result, gt, proj, result.shape[1], result.shape[0], -9999)
                
                