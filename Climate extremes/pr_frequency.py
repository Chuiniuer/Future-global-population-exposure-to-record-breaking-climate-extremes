import pandas as pd
import numpy as np
import netCDF4 as nc
from datetime import datetime
from osgeo import gdal, osr
import os
from tqdm import tqdm
import glob

os.environ['KMP_DUPLICATE_LIB_OK']='True'
"""
Year-by-year reading of multi-modal data to calculate annual rainfall-related indicators
"""

def get_arr_from_nc(path: str, var: str):
    """
    :param path: path
    :param var: variable
    :return: array
    """
    data = nc.Dataset(path)
    array = np.asarray(data.variables[var])
    return array

def coordinate_convert(array):
    array_temp = array.copy()
    array_temp[:, :, :720] = array[:, :, 720:1440]
    array_temp[:, :, 720: 1440] = array[:, :, :720]
    return array_temp

def unit_convert(array):
    """
    Kg/m3/s -> mm/day
    :param array:
    :return:
    """

    array[array==1e20] = -1
    array = (array * 86400)
    array[array<0] = -9999
    # array = array.astype(np.int16)
    return array

def wirte_2d_tif(extreme_index_name, output_folder, mode, senario, year, array):
    driver = gdal.GetDriverByName('GTiff')
    out_tif_name = output_folder + '\\' + extreme_index_name + '_' + mode + '_' + \
                   senario + "_" + year + '.tif'
    out_tif = driver.Create(out_tif_name, 1440, 600, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])

    # getting the image display range
    geotransform = (0, 0.25, 0, -60, 0, 0.25)  # ?????
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
class IndexExtreme:
    def __init__(self):
        pass

    def pr_by_years(self, array):
        return array.sum(axis=0)

    def pr_daily_extreme(self, array):
        return array.max(axis=0)

    def pr_50mm_count(self, array):
        return (array>=50).sum(axis=0)


if __name__ == "__main__":
    NEX_folder = "F:\\CMIP6_GDDP\\"

    mode_list = []
    files = [NEX_folder + d for d in os.listdir(NEX_folder) if d[-3:] == ".nc"]
    for s in files:
        s1 = s.split("\\")[-1]
        mode_list.append(s1.split("_")[2])
    mode_list = list(set(mode_list))


    Index = IndexExtreme()
    
    down_folder = "E:\\weather extreme\\pr1"
    down_files = [down_folder + d for d in os.listdir(down_folder) if d[-4:] == ".tif"]
    for m in tqdm(mode_list):
        #Get all the rainfall files for the corresponding model
        m_paths = []
        for f in files:
            if f.split("\\")[-1].split("_")[2] == m:
                m_paths.append(f)
        #Iteration of documents to calculate annual indicators
        for f in m_paths:
            try:
                senario = f.split("\\")[-1].split("_")[3]
                y = f.split("\\")[-1].split("_")[6].split(".")[0]
                arr = unit_convert(get_arr_from_nc(f, "pr"))
  
                arr_acc = Index.pr_by_years(arr)
                arr_acc[arr_acc<0] = -9999
                
                #Already generated files are skipped directly
                if down_folder + "accumulate_pr" + '_' + m + '_' + \
                    senario + y + '.tif' not in down_files:
                        wirte_2d_tif("accumulate_pr", down_folder, m, senario, y, arr_acc)
                
                if down_folder + "preceed50mm_pr" + '_' + m + '_' + \
                    senario + y + '.tif' not in down_files:
                    arr_50 = Index.pr_50mm_count(arr)
                    arr_50[arr_acc<0] = -9999
                    wirte_2d_tif("preceed50mm_pr", down_folder, m, senario, y, arr_50)
                if down_folder + "max_daily_pr" + '_' + m + '_' + \
                    senario + y + '.tif' not in down_files:
                    arr_daily = Index.pr_daily_extreme(arr)
                    arr_daily[arr_daily<0] = -9999
                    wirte_2d_tif("max_daily_pr", down_folder, m, senario, y, arr_daily)
            except:
                print(f"{m}, {y} has encountered a problem")

