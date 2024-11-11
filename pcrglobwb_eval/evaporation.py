from tqdm import tqdm
from pathlib import Path
import numpy as np
import xarray as xr
import time
import xarray as xr
import hydroeval as he
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
pd.set_option('display.max_columns', None)


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

class evaporation:
    def saveTimeSeries(fluxRawDataDirectory, saveFolder, simDirectory):
        """Convert grdc data downloaded from data portal into workable format           
        """   
        saveFolder = Path(saveFolder)
        simDirectory = Path(simDirectory)        
        fluxRawDataDirectory = Path(fluxRawDataDirectory)

        fluxMetaFile = fluxRawDataDirectory / 'meta.csv'
        fluxFiles = sorted((fluxRawDataDirectory / 'FLUXNET/data/extracted').glob('*.csv'))
        
        evaporationFiles = sorted(simDirectory.glob(f"**/*netcdf/totalEvaporation_dailyTot_output.nc"))
        
        def create_search_criteria(evaporationFiles, no_years=2):
            with xr.open_mfdataset(evaporationFiles, chunks={}, concat_dim='time', engine="netcdf4", combine='nested') as ds:
                time_min, time_max = ds.time.dt.year.min().values, ds.time.dt.year.max().values
                lat_min, lat_max = ds.lat.min().values.item(), ds.lat.max().values.item()
                lon_min, lon_max = ds.lon.min().values.item(), ds.lon.max().values.item()
                return {'time_min': time_min, 'time_max': time_max, 
                        'lat_min': lat_min,'lat_max': lat_max,
                        'lon_min': lon_min,'lon_max': lon_max,
                        'no_years': no_years}
            
        def findFluxSites(fluxMetaFile, search_criteria):
            def checkOverlap(search_criteria, result, yearsOverlap=3):
                #Check if the time ranges overlap
                if result['ts_end'] < search_criteria['time_min'] or result['ts_start'] > search_criteria['time_max']:
                    return False
                
                # Check if the spatial ranges overlap
                if search_criteria['lat_min'] < result['latitude'] < search_criteria['lat_max']:
                    if search_criteria['lon_min'] < result['longitude'] < search_criteria['lon_max']:
                        # Check if the difference between the start and end dates is at least 3 years
                        if (result['ts_end'] - search_criteria['time_min']) >= 3:
                            return True
                        else:
                            return False
                    else:
                        return False    
                else:
                    return False 

            meta_df = pd.read_csv(fluxMetaFile)
            df = pd.DataFrame(columns=['latitude', 'longitude', 'ts_start', 'ts_end'])
            for index, row in meta_df.iterrows():
                result = {'id': row['id'],
                    'latitude': row['lat'], 'longitude': row['lon'],
                    'ts_start': row['start-year'], 'ts_end': row['end-year']}
                if checkOverlap(search_criteria, result) == True:
                    df = pd.concat([df, pd.DataFrame(
                        {'id': [result['id']],
                         'latitude': [result['latitude']], 
                         'longitude':[result['longitude']],
                         'ts_start': [result['ts_start']], 
                         'ts_end': [result['ts_end']]})], ignore_index=True)
            return df
        
        def getIndices(evaporationFiles, metaInfo):
            lats = metaInfo['latitude'].values
            lons = metaInfo['longitude'].values
            with xr.open_mfdataset(evaporationFiles, combine="nested",concat_dim=["time"], decode_coords=False, 
                            combine_attrs="drop", decode_times=False, engine='h5netcdf',  decode_cf=False, parallel=True, chunks={}) as sim_ds:
                lon_indices = []
                for lon in lons:
                    lon_diff = np.abs(sim_ds.coords["lon"].values - lon)
                    lon_index = np.where(lon_diff == np.min(lon_diff))[0][0]
                    lon_indices.append(lon_index)
                lat_indices = []
                for lat in lats:
                    lat_diff = np.abs(sim_ds.coords["lat"].values - lat)
                    lat_index = np.where(lat_diff == np.min(lat_diff))[0][0]
                    lat_indices.append(lat_index)
                if lat_indices == [] or lon_indices == []:
                    return None, None
                lat_indices = xr.DataArray(lat_indices, dims='station', coords={"station": metaInfo['id'].values}).astype('int32')
                lon_indices = xr.DataArray(lon_indices, dims='station', coords={"station": metaInfo['id'].values}).astype('int32')
                
                
                ds = sim_ds.isel(time=0).load().isel(lon=lon_indices, lat=lat_indices)#.load()
                ds = ds.where(ds < 1e+20)
                ds = ds.dropna(dim='station', how='all')
                stations = ds.coords["station"].values
                lat_indices = lat_indices.sel(station=stations)
                lon_indices = lon_indices.sel(station=stations)
                return lat_indices, lon_indices, stations

        def read_netcdfs(evaporationFiles, lon_indices, lat_indices, step_size=1):
            data = xr.open_mfdataset(evaporationFiles, decode_cf=True,
                                   concat_dim='time', engine="h5netcdf", 
                                   combine='nested', combine_attrs='drop',
                                   join='override', parallel=True,
                                   mask_and_scale=False, 
                                   chunks={})['total_evaporation']
            timeSteps = [i for i in range(0, data.time.size+1, step_size)]
            for timeStep in tqdm(timeSteps, desc="Reading Points", disable=False):
                data_sub = data.isel(time=slice(timeStep,timeStep+step_size)).load().isel(lon=lon_indices, lat=lat_indices)
                if timeStep == 0: dsList = [data_sub]
                else: dsList.append(data_sub)
            selectedData = xr.concat(dsList, 'time')
            return selectedData

        def create_gpd(ds):
            df = pd.DataFrame({'lat': ds.lat.values, 'lon': ds.lon.values, 'station': ds.station.values})
            geometry = [Point(xy) for xy in zip(df.lon, df.lat)]
            gdf = gpd.GeoDataFrame(df, geometry=geometry)
            gdf = gdf.drop(['lat', 'lon'], axis=1)
            return gdf.to_file(saveFolder / f"points/points.shp")

        search_criteria = create_search_criteria(evaporationFiles)
        metaInfo = findFluxSites(fluxMetaFile, search_criteria)
        if len(metaInfo) != 0:
            lat_indices, lon_indices, stations = getIndices(evaporationFiles, metaInfo)
            sim_ds = read_netcdfs(evaporationFiles, lon_indices, lat_indices)
            if len(stations) != 0:
                (saveFolder / "timeseries").mkdir(parents=True, exist_ok=True)	
                (saveFolder / "points").mkdir(parents=True, exist_ok=True)	
                (saveFolder / "scores").mkdir(parents=True, exist_ok=True)	
                
                create_gpd(sim_ds)
                
                for ID in tqdm(stations, desc="Processing FLUXNET stations", disable=False):
                    filePath = [s for s in fluxFiles if ID in str(s)][0]
                    df = pd.read_csv(filePath)
                    df = df[['TIMESTAMP', 'LE_F_MDS']].rename(columns={'LE_F_MDS': 'obs_evap'})
                    metaInfoSubset = metaInfo[metaInfo['id'] == ID][['id', 'latitude', 'longitude']]
                    metaInfoSubset = pd.DataFrame(np.repeat(metaInfoSubset, df.shape[0], axis=0))
                    metaInfoSubset.columns = ['id', 'latitude', 'longitude']
                    df = pd.concat([metaInfoSubset, df], axis=1)
                    df['TIMESTAMP'] = pd.to_datetime(df['TIMESTAMP'].astype(str), format='%Y%m%d')
                    df = df.set_index('TIMESTAMP')
                    df = df.reindex(sim_ds.time.values)
                    df['obs_evap'][df['obs_evap']  < 0] = np.nan
                    non_nan_count = df['obs_evap'].count()
                    if non_nan_count > 1095:
                        df['sim_evap']=sim_ds.sel(station=ID).values * 1000
                        df['obs_evap'] = ((df['obs_evap'] * 86400) / 2257) / 1000 
                        first_non_missing, last_non_missing = df['obs_evap'].first_valid_index(), df['obs_evap'].last_valid_index()
                        df = df.loc[first_non_missing:last_non_missing]
                        df.to_parquet(saveFolder / f"timeseries/{ID}.parquet")

    def scoreTimeSeries(saveFolder):
        saveFolder = Path(saveFolder)
        timeSeriesFiles = sorted(saveFolder.glob(f"timeseries/*.parquet"))
        if len( timeSeriesFiles) != 0:
            df = pd.DataFrame(columns=['station', 'KGE', 'KGE_r', 'KGE_alpha', 'KGE_beta'])
            for file in tqdm(timeSeriesFiles, desc="Scoring Stations Daily", disable=False):
                station=file.stem
                file = pd.read_parquet(file)
                kge, r, alpha, beta = he.evaluator(he.kge, file["sim_evap"], file["obs_evap"])
                if df.empty:
                    df = pd.DataFrame({'station': [station], 
                                       'KGE': kge,
                                       'KGE_r': r,
                                       'KGE_alpha': alpha,
                                       'KGE_beta': beta})
                else:
                    df = pd.concat([df, pd.DataFrame({'station': [station], 
                                                      'KGE': kge,
                                                      'KGE_r': r,
                                                      'KGE_alpha': alpha,
                                                      'KGE_beta': beta})], ignore_index=True)
            scores = gpd.read_file(saveFolder / "points/points.shp")
            scores = scores.merge(df, on="station")
            scores.to_file(saveFolder / f"scores/day_scores.shp")
            
            df = pd.DataFrame(columns=['station', 'KGE', 'KGE_r', 'KGE_alpha', 'KGE_beta'])
            for file in tqdm(timeSeriesFiles, desc="Scoring Stations Monthly", disable=False):
                station=file.stem
                file = pd.read_parquet(file)
                file = file.resample('M').agg({'id': 'first',
                                                'latitude': 'first',
                                                'longitude': 'first',
                                                'obs_evap': 'mean',
                                                'sim_evap': 'mean'})
                kge, r, alpha, beta = he.evaluator(he.kge, file["sim_evap"], file["obs_evap"])
                if df.empty:
                    df = pd.DataFrame({'station': [station], 
                                       'KGE': kge,
                                       'KGE_r': r,
                                       'KGE_alpha': alpha,
                                       'KGE_beta': beta})
                else:
                    df = pd.concat([df, pd.DataFrame({'station': [station], 
                                                      'KGE': kge,
                                                      'KGE_r': r,
                                                      'KGE_alpha': alpha,
                                                      'KGE_beta': beta})], ignore_index=True)
            scores = gpd.read_file(saveFolder / "points/points.shp")
            scores = scores.merge(df, on="station")
            scores.to_file(saveFolder / "scores/month_scores.shp")