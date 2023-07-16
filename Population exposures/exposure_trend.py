from osgeo import gdal
import xarray as xr          
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import statsmodels.api as sm
from scipy import stats
from tqdm import tqdm
import os
# os.environ['PROJ_LIB'] = r'C:\ProgramData\Anaconda3\envs\geoplot\Lib\site-packages\pyproj\proj_dir\share\proj'
# Define the function that fits the linear equation and obtains the trend and p-value
count = 0
def lm_trend(x):
    if -9999 in list(x):
        return (-9999, -9999)

    if list(x) == [0, 0, 0, 0, 0, 0]:  
        return (-9999, -9999)
        
    years = np.arange(1,9)
    years = sm.add_constant(years)
    model = sm.OLS(x,years)
    result = model.fit()
    #print(result.summary())
    slope = result.params[1]
    p = result.pvalues[1]
    return(slope , p )
def saving_new_img(out_path, array, gt, proj, xsize, ysize, nodata):
    """
    :param out_path: 输出文件路径
    :param array: 输出文件数组
    :param gt: 仿射变换函数
    :param proj: 坐标信息
    :param xsize: x像元数
    :param ysize: y像元数
    :return:
    """
    #Re-generation of new documents
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
    
    # indices = ["R50", "RX1D", "PRCPTOT", "HW", "Tn90p", "Tx90p", "SFH", "CHDW"]
    indices = ["WD"]
    scenarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    periods = range(2020, 2100, 10)
    uncertainties = ["10th", "median", "90th"]
    
    for i in indices:
        for s in tqdm(scenarios):
            for u in uncertainties:
                
                if i == "R50":
                    i_copy = "R50"
                elif i == "RX1D":
                    i_copy = "max_daily_pr"
                elif i == "PRCPTOT":
                    i_copy = "accumulate_pr"
                elif i == "HW":
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
                
                
                input_folder = "G:/weather extreme/Revision/heat(tasmax)/exposure_span10/uncertainty/"
                ds = gdal.Open(input_folder + f"{i_copy}_{s}_2020_{u}_exposure.tif")
                slope_folder = "G:/weather extreme/Revision/heat(tasmax)/exposure_trend/slope/"
                p_folder = "G:/weather extreme/Revision/heat(tasmax)/exposure_trend/significance/"
                
                # try:
                #     gdal.Open(slope_folder+f"{i}_{s}_{u}_slope.tif").ReadAsArray()
                #     gdal.Open(p_folder+f"{i}_{s}_{u}_pvalue.tif").ReadAsArray()
                # except:
                
            
                xsize = ds.RasterXSize
                ysize = ds.RasterYSize
                # Size of each block after chunking (except for the last block)
                gt = ds.GetGeoTransform()
                proj = ds.GetProjection()
                data = np.array(
                                    [gdal.Open(input_folder + f"{i_copy}_{s}_{y}_{u}_exposure.tif").ReadAsArray()
                                     for y in periods])
                #         slope = np.zeros((data.shape[1], data.shape[2]), dtype=np.float32)
                #         p = np.zeros((data.shape[1], data.shape[2]), dtype=np.float32)
                #         for i in tqdm(range(data.shape[1])):
                #             for j in tqdm(range(data.shape[2])):
                #                 slope[i, j], p[i, j] = lm_trend(data[:, i, j])
                data[np.isnan(data)] = -9999
                
                data = xr.DataArray(data,dims = ['year','lat','lon'])
                    # 使用xarray的apply_ufunc方法计算
                data_lm_trend = xr.apply_ufunc(
                    lm_trend,
                    data,
                    input_core_dims = [['year']],
                    dask="parallelized",
                    output_core_dims = [[],[]],
                    output_dtypes=['float32', "float32"],
                    vectorize=True)
                saving_new_img(slope_folder+f"{i}_{s}_{u}_slope.tif", np.array(data_lm_trend[0]), gt, proj, xsize,
                                           ysize, -9999)
                saving_new_img(p_folder+f"{i}_{s}_{u}_pvalue.tif", np.array(data_lm_trend[1]), gt, proj, xsize,
                                           ysize, -9999)
