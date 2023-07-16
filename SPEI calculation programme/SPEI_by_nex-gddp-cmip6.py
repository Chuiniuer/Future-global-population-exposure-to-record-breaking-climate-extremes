# -*- coding: utf-8 -*-
"""
Created on Tue Dec  28 15:55:47 2022

@author: Bohao Li
"""
import numpy as np
import pandas as pd
from osgeo import gdal
from climate_indices import indices
from climate_indices import compute
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import datetime
import netCDF4 as nc

def cal_pet(tas_time_series, lat, data_start_year):
    """
    

    Parameters
    ----------
    tas_time_series : TYPE
        DESCRIPTION.
    lat_time_series : TYPE
        DESCRIPTION.
    data_start_year : TYPE
        DESCRIPTION.

    Returns
    -------
    pet_data : TYPE
        Return pet time series.

    """
    pet_data = indices.pet(temperature_celsius=tas_time_series,
                           latitude_degrees=lat,
                           data_start_year=data_start_year
        )
    return pet_data

def cal_spei(prcp_time_series, pet_time_series, scale, data_start_year, calibration_year_initial, calibration_year_final):
    '''

    Parameters
    ----------
    prcp_time_series : TYPE
        DESCRIPTION.
    pet_time_series : TYPE
        DESCRIPTION.
    scale : TYPE
        DESCRIPTION.
    data_start_year : TYPE
        DESCRIPTION.
    calibration_year_initial : TYPE
        DESCRIPTION.
    calibration_year_final : TYPE
        DESCRIPTION.

    Returns
    -------
    spei : TYPE
        Return time series

    '''
    spei = indices.spei(precips_mm=prcp_time_series,
                        pet_mm=pet_time_series,
                        scale=scale,
                        distribution=indices.Distribution.gamma,
                        periodicity=compute.Periodicity.monthly,
                        data_start_year=data_start_year,
                        calibration_year_initial=calibration_year_initial,
                        calibration_year_final=calibration_year_final
        )
    # spei[np.isnan(spei)] = -99
    return spei

def get_arr_from_nc(path: str, var: str):
    """
    :param path: path
    :param var: variable
    :return: array
    """
    data = nc.Dataset(path)
    array = np.asarray(data.variables[var])
    return array

def save_tiff(path, array, gt, proj, xsize, ysize, bands, nodata, eType=gdal.GDT_Float32):
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    outds = driver.Create(path, xsize=xsize, ysize=ysize, bands=bands, eType=eType)
    #gdalconst-variables can see the type
    outds.SetGeoTransform(gt)
    outds.SetProjection(proj)
    for i in range(bands):
        outband = outds.GetRasterBand(i+1)
        outband.WriteArray(array[i, :, :])
    outband.SetNoDataValue(nodata)
    outband.FlushCache()
    
    outband = None
    outds = None

def write_nc(path: str, arr: np.ndarray, main_var: str, unit: str, nodata: float):
    newfile = nc.Dataset(path, 'w', format='NETCDF4')

    #define dimensions
    lon = newfile.createDimension("longitude", size=1440)
    lat = newfile.createDimension("latitude", size=600)
    times = newfile.createDimension("time", size=None)

    # define variables for storing data
    lon = newfile.createVariable("lon", np.float32, dimensions="longitude")
    lat = newfile.createVariable("lat", np.float32, dimensions="latitude")
    time = newfile.createVariable("times", "S19", dimensions="time")
    crucial_var = newfile.createVariable(main_var, np.float32, dimensions=("time", "latitude", "longitude"))
    date_range = pd.date_range(datetime(1950, 1, 15), datetime(2100, 12, 31), freq="1m")

    # add data to variables
    lon[:] = np.arange(0.125, 360.125, 0.25)
    lat[:] = np.arange(89.875, -60.125, -0.25)
    crucial_var[:, :, :] = arr
    for i in range(arr.shape[0]):
        time[i] = date_range[i].strftime("%Y-%m-%d")
    print(date_range[0].strftime("%Y-%m-%d"))
    print(time)
    # add attributes
    # add global attributes
    newfile.title = f"monthly {main_var} data"
    newfile.start_time = time[i]
    newfile.times = time.shape[0]
    newfile.history = "Created" + datetime(2022, 12, 28).strftime("%Y-%m-%d")

    # add local attributes to variable
    lon.description = "longitude, range from 0 to 360"
    lon.units = "degrees"

    lat.description = "latitude, south is negative"
    lat.units = "degrees north"

    time.description = "time, unlimited dimension"
    time.units = "time since {0:s}".format(time[0])

    crucial_var.description = f"the time scaler is 3 month."
    crucial_var.units = unit
    crucial_var.missing_value = nodata
    # close file
    newfile.close()

def main():
    mode_list = ["ACCESS-CM2", "ACCESS-ESM1-5", "CNRM-ESM2-1",
                  "CMCC-ESM2", "CanESM5", "CNRM-CM6-1", "EC-Earth3",
                  "EC-Earth3-Veg-LR", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8",
                  "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L",
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
                  "NorESM2-MM", "UKESM1-0-LL"]
    senarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    
    input_tas_folder = "G:\\monthly_tas_aggregated\\"
    input_pr_folder = "H:\\weather extreme\\compound spei heatwave\\monthly_pr_aggregated\\"
    
    for m in mode_list:
        for se in tqdm(senarios):
            tempa_data = input_tas_folder + f"monthly_tas_{m}_{se}.nc"
            prcp_data = input_pr_folder + f"monthly_pr_{m}_{se}.nc"
            lat_data = r"H:\weather extreme\compound spei heatwave\latmask\latmask.tif"
            out_folder = "G:\\SPEI calculated\\"
            
            tempa_array = get_arr_from_nc(tempa_data, "tas")
            tempa_array[tempa_array==1e20] = np.nan
            
            pre_array = get_arr_from_nc(prcp_data, "pr")
            pre_array[pre_array==1e20] = np.nan
            
            lat_ds = gdal.Open(lat_data)
            lat_array = lat_ds.ReadAsArray()
            
            gt = lat_ds.GetGeoTransform()
            proj = lat_ds.GetProjection()
            
            pet = tempa_array.copy()
            
            for i in tqdm(range(tempa_array.shape[1])):
                for j in range(tempa_array.shape[2]):
                    pet[:, i, j] = cal_pet(tempa_array[:, i, j], lat_array[i, j], 1950)
                    
            print("___SPEI caculating___")
            #scale of spei
            for s in tqdm([3]):
                spei = tempa_array.copy()
                for i in range(tempa_array.shape[1]):
                    for j in range(tempa_array.shape[2]):
                        spei[:, i, j] = cal_spei(pre_array[:, i, j], pet[:, i, j], s, 1950, 1950, 2100)
                write_nc(out_folder + f"SPEI_{m}_{se}_{s}.nc", s, spei)
if __name__ == "__main__":
    main()