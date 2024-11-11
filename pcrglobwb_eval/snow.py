from tqdm import tqdm
from pathlib import Path
import xarray as xr
from functools import partial
import pandas as pd
import shutil
import re
import numpy as np
import time
from datetime import datetime
import pendulum as pend
import gc

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

class snow:
    @timing_decorator
    def preProcessSnowObs(saveFolder, snowCoverDataDirectory, gridFile):

        saveFolder = Path(saveFolder)
        saveFolder.mkdir(exist_ok=True, parents=True)
        snowCoverDataDirectory = Path(snowCoverDataDirectory)
        obsFiles = sorted(snowCoverDataDirectory.glob(f"**/*.nc"))
        obsFiles = obsFiles
        _chunks = {'time': 1, 'lat': 10000, 'lon': 10000}
        
        def create_zarr_store(file_paths, latGrid, lonGrid):
            def fileNames_to_dateIndex(file_paths):
                date_strings = [re.search(r"\d{8}", str(file_path)).group(0) for file_path in file_paths]
                dates = [datetime.strptime(date_string, "%Y%m%d") for date_string in date_strings]
                date_index = pd.DatetimeIndex(dates)                
                return date_index
            time_index = fileNames_to_dateIndex(file_paths)
            ds = xr.open_dataset(file_paths[0], chunks = _chunks)[['scfg']].rename({'scfg' : 'snow_cover'})
            ds = ds.reindex(time=time_index, lat=latGrid, lon=lonGrid, method='nearest')
            ds['lat'], ds['lon'] = ds.lat.astype(np.float32), ds.lon.astype(np.float32)
            ds.attrs = {}
            ds.encoding = {}
            ds['snow_cover'].attrs = {}
            ds['snow_cover'].encoding = {}
            ds = ds.chunk(_chunks)
            ds.to_zarr(saveFolder / 'snowCover30sec.zarr', compute=False, mode='w')

        grid = xr.open_dataset(gridFile, chunks = _chunks)
        latGrid= grid.lat.values
        lonGrid= grid.lon.values
        del grid
        create_zarr_store(obsFiles, latGrid, lonGrid)
        count= 0
        for file in tqdm(obsFiles, desc="Saving Obs Array", disable=False):
            count = count + 1
            with xr.open_dataset(file, cache=False, chunks = _chunks)['scfg'] as ds:
                ds = ds.rename('snow_cover')
                ds = xr.where(ds < 200, ds, np.nan)
                ds = xr.where(ds > 0, 1, ds)
                ds = ds.sortby('lat', ascending=False)
                ds = ds.reindex(lat=latGrid, lon=lonGrid, method='nearest')
                ds['lat'], ds['lon'] = ds.lat.astype(np.float32), ds.lon.astype(np.float32)
                ds = ds.chunk(_chunks)
                ds.drop_vars(['lat', 'lon']).to_zarr(saveFolder / f'snowCover30sec.zarr', region={"time": slice(count-1, count)})

    @timing_decorator
    def scoreSnowCover(saveFolder, simDirectory, snowCoverDataDirectory, partialRun=False, regionSubset=None):
        
        saveFolder = Path(saveFolder)
        if partialRun:
            saveFolder = saveFolder / f'partial/{regionSubset[0]}_{regionSubset[1]}_{regionSubset[2]}_{regionSubset[3]}'
        simDirectory = Path(simDirectory)
        snowCoverDataDirectory = Path(snowCoverDataDirectory)
        
        saveFolder.mkdir(exist_ok=True, parents=True)
        _chunks = {'time': 1, 'lat': -1, 'lon': -1}
        
        def readSim(simDirectory, partialRun):
            sweFiles = sorted(simDirectory.glob(f"**/*.nc"))
            data = xr.open_mfdataset(sweFiles, decode_cf=True,
                                concat_dim='time', engine="h5netcdf", 
                                combine='nested', combine_attrs='drop',
                                join='override', parallel=True,
                                mask_and_scale=True, 
                                chunks=_chunks)[['snow_water_equivalent']]
            data = data.rename({'snow_water_equivalent': 'snow_cover'})
            data['snow_cover'] = xr.where(data.snow_cover > 0, 1., data.snow_cover)
            data = data['snow_cover']
            mask = xr.open_dataset(sweFiles[0], decode_cf=True, engine="h5netcdf", mask_and_scale=True).isel(time=0).drop_vars('time')['snow_water_equivalent']
            mask.attrs = {}
            mask = xr.where(mask.notnull(), 1, np.nan)
            if partialRun:
                data = data.isel(lon=slice(regionSubset[0],regionSubset[1]), lat=slice(regionSubset[2],regionSubset[3]))
            return data, mask
        
        sim_ds, mask = readSim(simDirectory, partialRun)
        sim_ds = sim_ds.sel(time=slice('1985-01-01', '2019-12-31'))
        # sim_ds = sim_ds.sel(time=slice('2010-01-01', '2010-01-06'))
        obs_ds = xr.open_zarr(snowCoverDataDirectory, chunks=_chunks)['snow_cover'].sel(lon=slice(sim_ds.lon.values[0], sim_ds.lon.values[-1]), lat=slice(sim_ds.lat.values[0], sim_ds.lat.values[-1]))
        mask = mask.reindex(lat=obs_ds.lat, lon=obs_ds.lon, method='nearest').compute()
        
        sim_start = pend.from_format(sim_ds.time.values[0].astype(str)[:10], 'YYYY-MM-DD')
        sim_end = pend.from_format(sim_ds.time.values[-1].astype(str)[:10], 'YYYY-MM-DD')
        obs_start = pend.from_format(obs_ds.time.values[0].astype(str)[:10], 'YYYY-MM-DD')
        obs_end = pend.from_format(obs_ds.time.values[-1].astype(str)[:10], 'YYYY-MM-DD')
        if sim_start > obs_start:
            latest_start = sim_start
        else:
            latest_start = obs_start
        if sim_end < obs_end:
            earliest_end = sim_end
        else:
            earliest_end = obs_end
        
        latest_start,earliest_end = latest_start.to_date_string(), earliest_end.to_date_string()
        obs_ds = obs_ds.sel(time=slice(latest_start, earliest_end))
        sim_ds = sim_ds.sel(time=obs_ds.time)
        
        count = 0
        for t in tqdm(range(0, sim_ds.time.size)):
            count = count + 1
            sim_sub = sim_ds.isel(time=t).drop_vars('time').compute()
            obs_sub = obs_ds.isel(time=t).drop_vars('time').compute()
            sim_sub = sim_sub.reindex(lat=obs_sub.lat, lon=obs_sub.lon, method='nearest')
            sim_sub = sim_sub.where(obs_sub.notnull())
            obs_sub = obs_sub.where(sim_sub.notnull())
            
            sim_T_obs_T = sim_sub.where((sim_sub == 1) & (obs_sub == 1), np.nan)
            sim_T_obs_F = sim_sub.where((sim_sub == 1) & (obs_sub == 0), np.nan)
            
            sim_F_obs_T = sim_sub.where((sim_sub == 0) & (obs_sub == 1), np.nan)
            sim_F_obs_T = xr.where(sim_F_obs_T==0, 1, np.nan)
            
            sim_F_obs_F = sim_sub.where((sim_sub == 0) & (obs_sub == 0), np.nan)
            sim_F_obs_F = xr.where(sim_F_obs_F==0, 1, np.nan)
            
            if count == 1:
                obs_sub_final = obs_sub.fillna(0)
                sim_sub_final = sim_sub.fillna(0)
                sim_T_obs_T_final =sim_T_obs_T.fillna(0)
                sim_T_obs_F_final =sim_T_obs_F.fillna(0)
                sim_F_obs_T_final =sim_F_obs_T.fillna(0)
                sim_F_obs_F_final =sim_F_obs_F.fillna(0)
                
            else:
                obs_sub_final += obs_sub.fillna(0)
                sim_sub_final += sim_sub.fillna(0)
                sim_T_obs_T_final += sim_T_obs_T.fillna(0)
                sim_T_obs_F_final += sim_T_obs_F.fillna(0)
                sim_F_obs_T_final += sim_F_obs_T.fillna(0)
                sim_F_obs_F_final += sim_F_obs_F.fillna(0)
        
        # Calculate each metric as an xarray DataArray
        FAR_den = (sim_T_obs_T_final + sim_T_obs_F_final)
        FAR_den= xr.where(FAR_den==0, 1E-20, FAR_den)
        FAR  = sim_T_obs_F_final / FAR_den
 
        bias_den = (sim_T_obs_T_final + sim_T_obs_F_final + sim_F_obs_T_final + sim_F_obs_F_final)
        bias_den = xr.where(bias_den==0, 1E-20, bias_den)
        bias = (sim_T_obs_T_final + sim_T_obs_F_final) / bias_den


        pod_den = (sim_T_obs_T_final + sim_F_obs_T_final)
        pod_den = xr.where(pod_den==0, 1E-20, pod_den)
        pod  = sim_T_obs_T_final / pod_den
 
        SR_den = (sim_T_obs_T_final + sim_T_obs_F_final)
        SR_den = xr.where(SR_den==0, 1E-20, SR_den)
        SR   = sim_T_obs_T_final / SR_den

        # Combine the DataArrays into a single Dataset
        ds = xr.Dataset({'FAR' : FAR,
                         'bias': bias,
                         'pod' : pod,
                         'SR'  : SR,
                         'obs' : obs_sub_final,
                         'sim' : sim_sub_final,
                         'sim_T_obs_T' : sim_T_obs_T_final,
                         'sim_T_obs_F' : sim_T_obs_F_final,
                         'sim_F_obs_T' : sim_F_obs_T_final,
                         'sim_F_obs_F' : sim_F_obs_F_final})
        
        ds = ds.where(mask==1).astype(np.float32)
        ds.to_zarr(saveFolder / 'metrics.zarr', mode='w')
    
        data = np.array([sim_T_obs_T_final.sum().values,
                         sim_T_obs_F_final.sum().values, 
                         sim_F_obs_T_final.sum().values, 
                         sim_F_obs_F_final.sum().values])
        df = pd.DataFrame(data=[data], columns=['hits', 'false_alarms', 'misses', 'corr_negative'], index=[0])
        df['total'] = df['hits'] + df['false_alarms'] + df['misses'] + df['corr_negative']
        print(df)
        df.to_parquet(saveFolder / 'contingencyTable.parquet')