from tqdm import tqdm
from pathlib import Path
import numpy as np
import xarray as xr
import time
import xarray as xr
import xskillscore as xs

def timing_decorator(func):
    '''@timing_decorator ontop of function you want to time'''
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time_min = execution_time / 60
        print(f"\n{func.__name__} took: {execution_time_min:.1f}min / {execution_time:.1f}s to run ")
        return result    
    return wrapper

class TWS:
    def preprocess(saveFolder, simDirectory, graceLWEFile):
        
        def mergeSimFiles(clones, timeMin, timeMax, latGrid, lonGrid): 
                       
            count = 0
            for clone, files in tqdm(clones.items(), desc=f"Reading clones:"):
                count = count + 1
                maskArray =xr.open_mfdataset(files).sel(time=slice(timeMin, timeMax)).load()
                maskArray = maskArray.interp(lat=latGrid, lon=lonGrid)
                
                if count == 1: 
                    ds_final=maskArray
                else:
                    ds_final = ds_final.combine_first(maskArray)
            ds_final = ds_final.rename({'total_thickness_of_water_storage': 'lwe_thickness'})
            return ds_final        
        
        simDirectory = Path(simDirectory)
        saveFolder = Path(saveFolder)
        saveFolder.mkdir(parents=True, exist_ok=True)
    
        with xr.open_dataset(graceLWEFile, chunks={})[['lwe_thickness']] as grace_ds:
            timeMin, timeMax = np.datetime_as_string(grace_ds.time.min().values,unit='D'), \
                                    np.datetime_as_string(grace_ds.time.max().values, unit='D')
            latGrid = grace_ds.lat.values
            lonGrid = grace_ds.lon.values
            
        simFiles = simDirectory.glob(f'**/*netcdf/totalWaterStorageThickness_monthAvg_output.nc')
        clones = {}
        for nc_file in simFiles:
            parent_dir = str(nc_file.parents[2])
            
            if parent_dir in clones:
                clones[parent_dir].append(str(nc_file))
            else:
                clones[parent_dir] = [str(nc_file)]
        
        ds = mergeSimFiles(clones, timeMin, timeMax, latGrid=latGrid, lonGrid=lonGrid)
        ds = ds * 100
        ds = ds.reindex_like(grace_ds, method='nearest')
        ds.to_netcdf(saveFolder / 'sim_lwe_thicknessTemp.nc', mode='w')
        ds = xr.open_dataset(saveFolder / 'sim_lwe_thicknessTemp.nc', chunks='auto')['lwe_thickness']
        ds = ds - ds.mean(dim='time')
        ds.to_netcdf(saveFolder / 'sim_lwe_thickness.nc', mode='w')
        
        ds = xr.open_dataset(graceLWEFile, chunks='auto')['lwe_thickness']
        ds = ds - ds.mean(dim='time')
        ds.to_netcdf(saveFolder / 'obs_lwe_thickness.nc', mode='w')

    def scoreTWS(saveFolder, graceLWEFile):        
        saveFolder = Path(saveFolder)
        saveFolder.mkdir(parents=True, exist_ok=True)
        
        ds = xr.open_dataset(saveFolder / 'sim_lwe_thickness.nc')['lwe_thickness']
        grace = xr.open_dataset(saveFolder / 'obs_lwe_thickness.nc')['lwe_thickness']
        grace = grace.sel(time=slice(ds.time[0], ds.time[-1]))
        grace.attrs = {}
        
        ds_r2 = xs.r2(grace, ds, dim='time', skipna=False)
        ds_r2.to_dataset(name='r2').to_netcdf(saveFolder / 'r2.nc', mode='w')
        
        ds_rmse = xs.rmse(grace, ds, dim='time', skipna=True)
        ds_std = ds.std(dim='time', skipna=True)
        ds_rmse.to_dataset(name='rmse').to_netcdf(saveFolder / 'rmse.nc', mode='w')
        ds_rrmse = ds_rmse / ds_std
        ds_rrmse.to_dataset(name='rrmse').to_netcdf(saveFolder / 'rrmse.nc', mode='w')
        